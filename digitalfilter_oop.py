import numpy as np
import numpy_custom as npc
import prepare_traverse as preptrav
from scipy import interpolate
from scipy.signal import convolve
import matplotlib.pyplot as plt
import math
import os
import time
import nplotlib as nplt
import h5py

class Filter:

    def __init__(self, filename='none', res=0.1, outdir='cwd', dt_split=1, profile='hyperbolic-tangent',
                 fwidth=2, verbose=True, pdfr=np.sqrt(3), keep_all=False,
                 res_y=10, res_z=10, oy=0, oz=0, turb_prof='top-hat',
                 tu_intensity=0, lengthscale=3, y_range=(0, 0), z_range=(0, 0),
                 inner_diameter_norm=0, prof2d=False, prof1d=False, mdot=0, bulk_velocity=0,
                 density=1, plot_flag=False, polar=False, hdf5 = False, temp_fluct = 0, fluct_comps=('u', 'v', 'w'),
                 lengthscale_x=0, lengthscale_y=0, lengthscale_z=0):

        ################################################################################################################
        ####################################### Setting Parameters #####################################################
        ################################################################################################################

        # time step split for interpolation between time steps: dt_split = number of interpolation time steps between
        # calculated time steps
        self.dt_split = dt_split
        # activating plot of fields
        self.plot_flag = plot_flag
        # activates Digital Filter in polar coordinates. Ture = polar, False = cartesian
        self.polar = polar
        # activates hdf5 export
        self.hdf5 = hdf5
        # standard deviation of temperature fluctuations
        self.temp_fluct = temp_fluct
        # filter width
        self.fwidth = fwidth
        # resolution
        self.res = res
        # bulk velocity in x direction
        self.u_bulk = bulk_velocity
        # Variables to create fluctuations for
        self.fluct_comps = fluct_comps
        # save profiles
        self.verbose = verbose
        if outdir is not 'cwd':
            self.outdir = outdir
        else:
            self.outdir = os.getcwd()
        # create output folder
        if self.verbose or self.plot_flag[0]:
            if not os.path.exists(self.outdir + '/PODFS'):
                os.makedirs(self.outdir + '/PODFS')

        # initialising parameters for statistics
        self.mean_diff = np.zeros((0, len(fluct_comps) + 1))
        self.max_diff = np.zeros((0, len(fluct_comps) + 1))

        ################################################################################################################
        ######################################## Preparation of input field ############################################
        ################################################################################################################

        # switches between 2D field data as input (prof2d = True) and 1D data or no input data (prof2d = False)
        # reads in profiles or creates profiles from scratch
        if prof2d:
            self.profile = Profil2D(filename, res, mdot, bulk_velocity, density, polar, fluct_comps=fluct_comps,
                                    temp_fluct=temp_fluct, lengthscale_x=lengthscale_x, lengthscale_y=lengthscale_y,
                                    lengthscale_z=lengthscale_z)
        else:
            self.profile = Profil1DProf(res, profile=profile, turb_prof=turb_prof, u_bulk=bulk_velocity,
                                        tu_intensity=tu_intensity, filename=filename,
                                        res_y=res_y, res_z=res_z, oy=oy, oz=oz,
                                        lengthscale=lengthscale, y_range=y_range, z_range=z_range,
                                        inner_diameter_norm=inner_diameter_norm, prof1d=prof1d)

        # scaling of profiles, interpolation on output grid, calculation of RST
        self.profile.prepare_for_filter()

        # plot of prepared input fields
        if self.plot_flag[0] and self.plot_flag[1]:
            self.profile.plot_input_profile(outdir=self.outdir)

        ################################################################################################################
        ######################################## Preparation of filter #################################################
        ################################################################################################################

        # definition of filter width in 3 directions
        self.nfx = self.profile.lnx * self.fwidth
        self.nfy = self.profile.lny * self.fwidth
        self.nfz = self.profile.lnz * self.fwidth

        ######################################## Filter Matrix #########################################################

        # initialising field of filter matrix, format: self.filter_coeffs[y_index][z_index]
        self.filter_coeffs = [[[] for i in range(self.nfy.shape[0])] for j in range(self.nfy.shape[1])]

        # fill filter matrix with filter coefficients
        self.fill_coeff_field()

        ######################################## Random Fields #########################################################

        # initialisation of random fields for all components in fluct_comps
        self.rand_fields = [RandomField(name, self.nfx, self.nfy, self.nfz, pdfr, polar) for name in self.fluct_comps]

        ######################################## Correlated Fields #####################################################

        # initialisation of correlated fields for all components in random fields
        self.corr_fields = [CorrField(el.name, el, self.filter_coeffs, self.profile.outside_flag, keep_all)
                            for el in self.rand_fields]

        ######################################## Output Fields #########################################################

        # initialisation of output fields for all components in random fields
        self.out_fields = [OutField(el.name, el, self.dt_split, keep_all) for el in self.rand_fields]

        # define and initialise additional variables fpr HYDRA export
        add_vars = ['awhirl', 'apitch', 'ptotal']
        self.out_fields_2 = [OutField(el, self.rand_fields[0], self.dt_split, keep_all) for el in add_vars]
        self.fluct_comps_2 = ['awhirl', 'apitch', 'ptotal']

        # getting filter timestep from resolution and velocity
        self.dt_filter = self.res / \
                         np.mean(self.profile.dataarrays[self.profile.get_datanames().index('u')].data[
                                     np.logical_not(self.profile.outside_flag)])



    def fill_coeff_field(self):
        # At the moment same coefficient matrix for all grid points ==> maybe make different for lenght scale distribution
        matrix = FilterMatrix(self.nfx[0, 0], self.nfy[0, 0], self.nfz[0, 0],
                              self.profile.lnx[0, 0], self.profile.lny[0, 0], self.profile.lnz[0, 0])
        for i in range(len(self.filter_coeffs)):
            for j in range(len(self.filter_coeffs[i])):
                self.filter_coeffs[i][j] = matrix

    ####################################################################################################################
    ################################################# Filter Process ###################################################
    ####################################################################################################################

    def filtering(self, timesteps, periodic=True, hdf5 = False, PODFS = False, conv_stat=False, time_summary=False):
        if conv_stat:
            fig = plt.figure()

        # calculate time step for output fields
        if self.dt_split > 1:
            timesteps = math.ceil(timesteps / self.dt_split)

        # create HDF5 file
        if self.verbose and hdf5:
            self.create_hdf5_file(timesteps)

        # allocate array for PODFS inlet
        if PODFS:
            if hdf5:
                coord_comps = ['phi', 'r']
                coord = [self.profile.dataarrays[self.profile.get_datanames().index(ident)].data for ident in coord_comps]
                coord = [el.reshape(el.size) for el in coord]
                coord[0] = coord[0] * 180 / np.pi
            else:
                coord_comps = ['x', 'y', 'z']
                coord = [self.profile.dataarrays[self.profile.get_datanames().index(ident)].data for ident in
                         coord_comps]
                coord = [el.reshape(el.size) for el in coord]
            self.prepare_POD = Prepare_POD(hdf5, timesteps, self.out_fields[0].res_y, self.out_fields[0].res_z, self.fluct_comps, self.fluct_comps_2, self.dt_filter, coord, coord_comps)

        ########################################### loop over timesteps ################################################
        starttime = time.time()
        for i in range(0, timesteps):
            print('Processing Inlet {:d}/{:d}'.format(i + 1, timesteps))
            t0 = time.time()
            # creation of spatial correlations
            for el in self.corr_fields:
                el.add_field_conv()
            t1 = time.time()
            # adapt profile to Reynolds stresses
            self.profile.rst_transform.adapt_2_profile(self.corr_fields,
                                                       self.out_fields,
                                                       self.profile,
                                                       self.out_fields_2)
            t2 = time.time()

            # save profiles
            if self.verbose:
                self.save_field(i, hdf5)

            # prepare data for PODFS
            if PODFS:
                Prepare_POD.prepare_for_POD(self.prepare_POD, self.out_fields, self.out_fields_2, i)

            t3 = time.time()

            # prepare random fields for next time step
            for el in self.rand_fields:
                el.next_field(timesteps, i, periodic)
            t4 = time.time()

            # convergence statistics
            if conv_stat:
                self.diff_stats()
                self.plot_diff_stats(fig)
                if self.verbose and (i + 1) % 10 == 0:
                    fig2 = self.plot_k_field(i + 1)
                    fig2.savefig(self.outdir + '/PODFS/k_field_{:04d}.png'.format(i + 1),
                                 dpi=200)
                    plt.close(fig2)

            # Statistics during filtering
            if (i+1) % 50 == 0:
                self.save_stats(i + 1, periodic)

            # approximation of remaining time
            endtime = time.time()
            avg_time = (endtime - starttime) / (i + 1)
            remain_time = (timesteps - (i + 1)) * avg_time
            print('Verbleibende Zeit: {:02d}h {:02d}m {:02d}s'.format(int(remain_time // 3600),
                                                                      int((remain_time % 3600) // 60),
                                                                      int(remain_time % 60)))

            if time_summary:
                print('Time Summary: {: 8.3f} s\n'
                      'Spatial Correlations: {:4.2f}\n'
                      'Adapt to RST: {:4.2f}\n'
                      'Saving: {:4.2f}\n'
                      'New Random: {:4.2f}\n'
                      'ConvStats: {:4.2f}'.format(endtime - t0,
                                                  (t1 - t0) / (endtime - t0),
                                                  (t2 - t1) / (endtime - t0),
                                                  (t3 - t2) / (endtime - t0),
                                                  (t4 - t3) / (endtime - t0),
                                                  (endtime - t4) / (endtime - t0)))

            ############################################ plot fields ###################################################

            # plot random fields
            if self.plot_flag[0] and self.plot_flag[2]:
                step = str("{:04d}".format(i + 1))
                if self.polar:
                    phi = np.concatenate((self.profile.dataarrays[self.profile.get_datanames().index('phi')].data[:,:91]-22.5*np.pi/180, self.profile.dataarrays[self.profile.get_datanames().index('phi')].data[:,:91], self.profile.dataarrays[self.profile.get_datanames().index('phi')].data[:,:91]+22.5*np.pi/180), axis=1)
                    r = np.concatenate((self.profile.dataarrays[self.profile.get_datanames().index('r')].data[:,:91], self.profile.dataarrays[self.profile.get_datanames().index('r')].data[:,:91], self.profile.dataarrays[self.profile.get_datanames().index('r')].data[:,:91]), axis=1)
                    y = r * -np.sin(phi)
                    z = r * np.cos(phi)

                    nplt.contourf2('rand' + step, '-', y,
                                   z,
                                   np.concatenate((self.rand_fields[0].field[0][0:self.profile.res_z, 0:self.profile.res_y-2], self.rand_fields[0].field[0][0:self.profile.res_z, 0:self.profile.res_y-2], self.rand_fields[0].field[0][0:self.profile.res_z, 0:self.profile.res_y-2]), axis=1) , 100, 'y', 'z',
                                   'r', self.outdir + '/rand' + step,
                                   figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                    nplt.close('rand' + step)

                else:
                    nplt.contourf2('rand' + step, '-', self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                                   self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                                   self.rand_fields[0].field[0][0:self.profile.res_z, 0:self.profile.res_y], 100, 'y', 'z',
                                   'r', self.outdir + '/rand' + step,
                                   figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                    nplt.close('rand' + step)

            # plot correlated fields
            if self.plot_flag[0] and self.plot_flag[3]:
                step = str("{:04d}".format(i + 1))
                if self.polar:

                    if step == str("{:04d}".format( 1)):
                        u_min = np.min(self.corr_fields[0].field[0][:,:91])
                        u_max = np.max(self.corr_fields[0].field[0][:,:91])
                        v_tan_min = np.min(self.corr_fields[1].field[0][:,:91])
                        v_tan_max = np.max(self.corr_fields[1].field[0][:,:91])
                        v_rad_min = np.min(self.corr_fields[2].field[0][:,:91])
                        v_rad_max = np.max(self.corr_fields[2].field[0][:,:91])

                    nplt.contourf(u_min, u_max,
                                  'filtered_u' + step, y,
                                  z,
                                  np.concatenate((self.corr_fields[0].field[0][:,:91], self.corr_fields[0].field[0][:,:91], self.corr_fields[0].field[0][:,:91]), axis=1), 100, 'y', 'z',
                                  'u [m/s]', self.outdir + '/filtered_u' + step,
                                  figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                    nplt.close('filtered_u' + step)

                    nplt.contourf(v_tan_min, v_tan_max,
                                  'filtered_v_tan' + step, y,
                                  z,
                                  np.concatenate((self.corr_fields[1].field[0][:,:91], self.corr_fields[1].field[0][:,:91], self.corr_fields[1].field[0][:,:91]), axis=1), 100, 'y', 'z',
                                  'v_tan [m/s]', self.outdir + '/filtered_v_tan' + step,
                                  figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                    nplt.close('filtered_v_tan' + step)

                    nplt.contourf(v_rad_min, v_rad_max,
                                  'filtered_v_rad' + step, y,
                                  z,
                                  np.concatenate((self.corr_fields[2].field[0][:,:91], self.corr_fields[2].field[0][:,:91], self.corr_fields[2].field[0][:,:91]), axis=1), 100, 'y', 'z',
                                  'v_rad [m/s]', self.outdir + '/filtered_v_rad' + step,
                                  figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                    nplt.close('filtered_v_rad' + step)

                    if 'Tt' in self.fluct_comps:
                        if step == str("{:04d}".format(1)):
                            Tt_min = np.min(self.corr_fields[3].field[0][:,:91])
                            Tt_max = np.max(self.corr_fields[3].field[0][:,:91])

                        nplt.contourf(Tt_min, Tt_max,
                                  'filtered_Tt' + step, y,
                                  z,
                                  np.concatenate((self.corr_fields[3].field[0][:,:91], self.corr_fields[3].field[0][:,:91], self.corr_fields[3].field[0][:,:91]), axis=1), 100, 'y', 'z',
                                  'Tt [K]', self.outdir + '/filtered_Tt' + step,
                                  figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                        nplt.close('filtered_Tt' + step)

                else:

                    nplt.contourf(np.min(self.corr_fields[0].field[0]), np.max(self.corr_fields[0].field[0]),
                                  'filtered_u' + step, self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                                  self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                                  self.corr_fields[0].field[0], 100, 'y', 'z',
                                  'u [m/s]', self.outdir + '/filtered_u' + step,
                                  figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                    nplt.close('filtered_u' + step)

                    nplt.contourf(np.min(self.corr_fields[1].field[0]), np.max(self.corr_fields[1].field[0]),
                                  'filtered_v' + step, self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                                  self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                                  self.corr_fields[1].field[0], 100, 'y', 'z',
                                  'v [m/s]', self.outdir + '/filtered_v' + step,
                                  figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                    nplt.close('filtered_v' + step)

                    nplt.contourf(np.min(self.corr_fields[2].field[0]), np.max(self.corr_fields[2].field[0]),
                                  'filtered_w' + step, self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                                  self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                                  self.corr_fields[2].field[0], 100, 'y', 'z',
                                  'w [m/s]', self.outdir + '/filtered_w' + step,
                                  figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                    nplt.close('filtered_w' + step)

                    if 'Tt' in self.fluct_comps:
                        nplt.contourf(np.min(self.corr_fields[3].field[0]), np.max(self.corr_fields[3].field[0]),
                                  'filtered_Tt' + step, self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                                  self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                                  self.corr_fields[3].field[0], 100, 'y', 'z',
                                  'Tt [K]', self.outdir + '/filtered_Tt' + step,
                                  figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                        nplt.close('filtered_Tt' + step)

            # plot adapted fields
            if self.plot_flag[0] and self.plot_flag[4]:
                step = str("{:04d}".format(i + 1))
                if self.polar:
                    if step == str("{:04d}".format(1)):
                        u_min2 = np.min(self.out_fields[0].field[-1][:,:91])
                        u_max2 = np.max(self.out_fields[0].field[-1][:,:91])
                        v_tan_min2 = np.min(self.out_fields[1].field[0][:,:91])
                        v_tan_max2 = np.max(self.out_fields[1].field[0][:,:91])
                        v_rad_min2 = np.min(self.out_fields[2].field[0][:,:91])
                        v_rad_max2 = np.max(self.out_fields[2].field[0][:,:91])
                        whirl_min = np.min(self.out_fields_2[0].field[-1][:,:91])
                        whirl_max = np.max(self.out_fields_2[0].field[-1][:,:91])
                        pitch_min = np.min(self.out_fields_2[1].field[-1][:,:91])
                        pitch_max = np.max(self.out_fields_2[1].field[-1][:,:91])
                        ptotal_min = np.min(self.out_fields_2[2].field[-1][:,:91])
                        ptotal_max = np.max(self.out_fields_2[2].field[-1][:,:91])

                    nplt.contourf(u_min2, u_max2,
                                  'adapted_u' + step, y,
                                  z,
                                  np.concatenate((self.out_fields[0].field[-1][:,:91], self.out_fields[0].field[-1][:,:91],self.out_fields[0].field[-1][:,:91]),axis=1) , 100, 'y', 'z',
                                  'u [m/s]', self.outdir + '/adapted_u' + step,
                                  figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                    nplt.close('adapted_u' + step)

                    nplt.contourf(v_tan_min2, v_tan_max2,
                                  'adapted_v_tan' + step, y,
                                  z,
                                  np.concatenate((self.out_fields[1].field[-1][:,:91], self.out_fields[1].field[-1][:,:91],self.out_fields[1].field[-1][:,:91]),axis=1) , 100, 'y', 'z',
                                  'v_tan [m/s]', self.outdir + '/adapted_v_tan' + step,
                                  figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                    nplt.close('adapted_v_tan' + step)

                    nplt.contourf(v_rad_min2, v_rad_max2,
                                  'adapted_v_rad' + step, y,
                                  z,
                                  np.concatenate((self.out_fields[2].field[-1][:,:91], self.out_fields[2].field[-1][:,:91],self.out_fields[2].field[-1][:,:91]),axis=1) , 100, 'y', 'z',
                                  'v_rad [m/s]', self.outdir + '/adapted_v_rad' + step,
                                  figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                    nplt.close('adapted_v_rad' + step)

                    nplt.contourf(whirl_min, whirl_max,
                                  'adapted_whirl' + step, y,
                                  z,
                                  np.concatenate((self.out_fields_2[0].field[-1][:,:91], self.out_fields_2[0].field[-1][:,:91],self.out_fields_2[0].field[-1][:,:91]),axis=1) , 100, 'y', 'z',
                                  'whirl [°]', self.outdir + '/adapted_whirl' + step,
                                  figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                    nplt.close('adapted_whirl' + step)

                    nplt.contourf(pitch_min, pitch_max,
                                  'adapted_pitch' + step, y,
                                  z,
                                  np.concatenate((self.out_fields_2[1].field[-1][:,:91], self.out_fields_2[1].field[-1][:,:91],self.out_fields_2[1].field[-1][:,:91]),axis=1) , 100, 'y', 'z',
                                  'pitch [°]', self.outdir + '/adapted_pitch' + step,
                                  figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                    nplt.close('adapted_pitch' + step)

                    nplt.contourf(ptotal_min, ptotal_max,
                                  'adapted_ptotal' + step, y,
                                  z,
                                  np.concatenate((self.out_fields_2[2].field[-1][:,:91], self.out_fields_2[2].field[-1][:,:91],self.out_fields_2[2].field[-1][:,:91]),axis=1) , 100, 'y', 'z',
                                  'pt [Pa]', self.outdir + '/adapted_ptotal' + step,
                                  figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                    nplt.close('adapted_ptotal' + step)

                    if 'Tt' in self.fluct_comps:
                        if step == str("{:04d}".format(1)):
                            Tt_min2 = np.min(self.out_fields[3].field[-1][:,:91])
                            Tt_max2 = np.max(self.out_fields[3].field[-1][:,:91])
                        nplt.contourf(Tt_min2, Tt_max2,
                                  'adapted_Tt' + step, y,
                                  z,
                                  np.concatenate((self.out_fields[3].field[-1][:,:91], self.out_fields[3].field[-1][:,:91],self.out_fields[3].field[-1][:,:91]),axis=1) , 100, 'y', 'z',
                                  'Tt [K]', self.outdir + '/adapted_Tt' + step,
                                  figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                        nplt.close('adapted_Tt' + step)

                else:
                    nplt.contourf(np.min(self.out_fields[0].field[0]), np.max(self.out_fields[0].field[0]),
                                  'adapted_u' + step, self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                                  self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                                  self.out_fields[0].field[0], 100, 'y', 'z',
                                  'u [m/s]', self.outdir + '/adapted_u' + step,
                                  figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                    nplt.close('adapted_u' + step)

                    nplt.contourf(np.min(self.out_fields[1].field[0]), np.max(self.out_fields[1].field[0]),
                                  'adapted_v' + step, self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                                  self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                                  self.out_fields[1].field[0], 100, 'y', 'z',
                                  'v [m/s]', self.outdir + '/adapted_v' + step,
                                  figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                    nplt.close('adapted_v' + step)

                    nplt.contourf(np.min(self.out_fields[2].field[0]), np.max(self.out_fields[2].field[0]),
                                  'adapted_w' + step, self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                                  self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                                  self.out_fields[2].field[0], 100, 'y', 'z',
                                  'w [m/s]', self.outdir + '/adapted_w' + step,
                                  figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                    nplt.close('adapted_w' + step)

                    if 'Tt' in self.fluct_comps:
                        nplt.contourf(np.min(self.out_fields[3].field[0]), np.max(self.out_fields[3].field[0]),
                                  'adapted_Tt' + step, self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                                  self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                                  self.out_fields[3].field[0], 100, 'y', 'z',
                                  'Tt [K]', self.outdir + '/adapted_Tt' + step,
                                  figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                        nplt.close('adapted_Tt' + step)

                # nplt.contourf(np.min(self.out_fields[0].field[-1][:,:91]), np.max(self.out_fields[0].field[-1][:,:91]),
                #               'adapted' + step, y,
                #               z,
                #               np.concatenate((self.out_fields[0].field[-1][:,:91], self.out_fields[0].field[-1][:,:91],self.out_fields[0].field[-1][:,:91]),axis=1) , 100, 'y', 'z',
                #               'u', self.outdir + '/adapted' + step,
                #               figsize=(16, 16 * ((self.profile.res_z-3) / (self.profile.res_y-3)/3)))
                # nplt.close('adapted' + step)

        # store first fields in last fields for temporal periodic conditions
        if periodic and self.dt_split > 1:
            for el in self.out_fields:
                el.add_last_fields()
            if self.verbose:
                self.save_field(timesteps, hdf5, last=True)

        # save settings file
        if self.verbose or self.plot_flag[0]:
            self.save_stats(timesteps, periodic)
            # save convergence plot
            if conv_stat:
                fig.savefig(self.outdir + '/PODFS/convStats.png')

    def save_stats(self, timesteps, periodic):
        with open(self.outdir + '/PODFS/00_params_{:06d}.txt'.format(timesteps * self.dt_split), 'w') as out_file:
            out_file.write('Datei erstellt: ' + time.ctime() + '\n')
            out_file.write('Zeitschritt Filter: {:14e}\n'.format(self.dt_filter))
            out_file.write('Interpolation Filterwerte: {:d}\n'.format(self.dt_split))
            out_file.write('Zeitschritt LES: {:14e}\n'.format(self.dt_filter / self.dt_split))
            out_file.write('Inputprofil: ' + self.profile.name + '\n')
            out_file.write('Aufloesung: {:.8f}\n'.format(self.res))
            out_file.write('Filterweite: {:d}\n'.format(self.fwidth))
            out_file.write('Lengthscale: {:8e} m\n'.format(self.profile.lengthscale_x))
            out_file.write('Lengthscale / Resolution: {:d}\n'.format(int(self.profile.lnx[0, 0])))
            out_file.write('Anzahl Zeitschritte Filter: {:d}\n'.format(timesteps))
            out_file.write('Anzahl Zeitschritte Dateien: {:d}\n'.format(timesteps * self.dt_split))
            out_file.write('Periodisch: {:b}\n'.format(periodic))
            out_file.write('Anzahl Zeilen: {:d}\n'.format(self.out_fields[0].get_last_field().size))
            out_file.write('\n')
            u_mean_f = self.get_mean_diff()
            k_mean_f = self.get_k_mean_diff()
            u_max_f = self.get_max_diff()
            k_max_f = self.get_k_max_diff()
            val_range = [(self.profile.get_data(el).min(), self.profile.get_data(el).max()) for el in ['u', 'v',
                                                                                                       'w', 'k']]
            out_file.write('Difference to Input Profile:\n'
                           'Norm: <val> / (Max_Profil - Min_Profil)\n'
                           'ID \tMean:     \tMax:      \tMin_Profil\tMax_Profil\tNorm_Mean_diff\tNorm_Max_diff\n'
                           'u: \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e} \n'
                           'v: \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e} \n'
                           'w: \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e} \n'
                           'k: \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e}  \t{:8e} \n'.
                           format(u_mean_f[0], u_max_f[0], val_range[0][0], val_range[0][1],
                                  u_mean_f[0] / (val_range[0][1] - val_range[0][0]),
                                  u_max_f[0] / (val_range[0][1] - val_range[0][0]),
                                  u_mean_f[1], u_max_f[1], val_range[1][0], val_range[1][1],
                                  u_mean_f[1] / (val_range[1][1] - val_range[1][0]),
                                  u_max_f[1] / (val_range[1][1] - val_range[1][0]),
                                  u_mean_f[2], u_max_f[2], val_range[2][0], val_range[2][1],
                                  u_mean_f[2] / (val_range[2][1] - val_range[2][0]),
                                  u_max_f[2] / (val_range[2][1] - val_range[2][0]),
                                  k_mean_f, k_max_f, val_range[3][0], val_range[3][1],
                                  k_mean_f / (val_range[3][1] - val_range[3][0]),
                                  k_max_f / (val_range[3][1] - val_range[3][0])
                                  ))
        fig_k = self.plot_k_field()
        fig_k.savefig(self.outdir + '/PODFS/00_k_field_{:06d}.png'.format(timesteps * self.dt_split))
        plt.close(fig_k)

    def save_field(self, timestep, hdf5, last=False):
    ######################################## save snapshots ############################################################
        if hdf5:
            # for HYDRA export use r and phi as coordinates
            coord = [self.profile.dataarrays[self.profile.get_datanames().index(ident)].data for ident in ['phi', 'r']]
            coord = [el.reshape(el.size) for el in coord]
            coord[0] = coord[0] * 180 / np.pi
        else:
            # for cfx export use x,y,z as coordinates
            coord = [self.profile.dataarrays[self.profile.get_datanames().index(ident)].data for ident in ['x', 'y', 'z']]
            coord = [el.reshape(el.size) for el in coord]
        if timestep == 0:
            start_ind = -1
        else:
            start_ind = -self.dt_split
        if not last:
            end_ind = 0
        else:
            end_ind = -1

        for n in range(start_ind, end_ind):
            # save HDF5 file
            if hdf5:
                whirl = self.out_fields_2[0].field[n]
                whirl = whirl.reshape(whirl.size) * 180 / np.pi
                pitch = self.out_fields_2[1].field[n]
                pitch = pitch.reshape(pitch.size) * 180 / np.pi
                ptotal = self.out_fields_2[2].field[n]
                ptotal = ptotal.reshape(ptotal.size)
                if 'Tt' in self.fluct_comps:
                    ttotal = self.out_fields[3].field[n]
                else:
                    ttotal = self.profile.dataarrays[self.profile.get_datanames().index('Tt')].data
                ttotal = ttotal.reshape(ttotal.size)
                turbke = self.profile.dataarrays[self.profile.get_datanames().index('k')].data
                turbke = turbke.reshape(ttotal.size)
                turbepsilon = self.profile.dataarrays[self.profile.get_datanames().index('e')].data
                turbepsilon = turbepsilon.reshape(ttotal.size)
                turbomega = turbepsilon / turbke
                data = np.transpose(np.array((coord[0],
                                            coord[1],
                                            ptotal,
                                            ttotal,
                                            whirl,
                                            pitch,
                                            turbke,
                                            turbomega)))
                f = h5py.File(self.outdir+ '/PODFS/Inlet.hdf5', 'r+')
                # write in 'main' group
                main = f.require_group("main")
                data = main.create_dataset("snapshot_"+str(timestep * self.dt_split + n + 2), data = data)

                # # export of snapshots for HYDRA
                # with open(self.outdir + '/PODFS/inlet.dat' + str(timestep * self.dt_split + n + 2), 'w') as out_file:
                #     out_file.write('theta, radius, ptotal, ttotal, awhirl, apitch, turbke, turbomega \n')
                #     for m in range(0, coord[0].size):
                #         out_file.write('{:.12f} {:.12f} {:.12f} {:.12f} {:.12f} {:.12f} {:.12f} {:.12f}\n'.format(
                #            coord[0][m], coord[1][m], ptotal[m], ttotal[m], whirl[m], pitch[m], turbke[m], turbomega[m]))

            else:
            # save .dat file
                velo = [el.get_indexed_field(n) for el in self.out_fields]
                velo = [el.reshape(el.size) for el in velo]

                with open(self.outdir + '/PODFS/inlet.dat' + str(timestep * self.dt_split + n + 2), 'w') as out_file:
                    for m in range(0, coord[0].size):
                        outline = ' '.join(['{:.12f}'.format(coord[comp][m]) for comp in range(3)] +
                                           ['{:.12f}'.format(velo[self.fluct_comps.index(name)][m])
                                            for name in self.fluct_comps] + ['\n'])
                        #out_file.write('{:.12f} {:.12f} {:.12f} {:.12f} {:.12f} {:.12f}\n'.format(
                        #    coord[0][m], coord[1][m], coord[2][m], velo[0][m], velo[1][m], velo[2][m]
                        #))
                        out_file.write(outline)

    def create_hdf5_file(self, timesteps):
        f = h5py.File(self.outdir+ '/PODFS/Inlet.hdf5', 'w')
        main = f.create_group("main")
        main.attrs.create("Number of Snapshots", timesteps * self.dt_split)
        main.attrs.create("Number of Points",len(self.profile.fr.data_columns[0].data))
        main.attrs.create("Number of Variables",8)
        main.attrs.create("Variables", np.string_("theta, radius, ptotal, ttotal, awhirl, apitch, turbke, turbomega"))
        main.attrs.create("Timestep width", self.dt_filter/self.dt_split)
        f.close()

    def plot_avg_field(self, fig='new', timestep=0):
        if fig is 'new':
            fig = plt.figure()
        h = plt.quiver(self.profile.get_data('y'), self.profile.get_data('z'),
                       self.out_fields[1].get_average_until_timestep(timestep),
                       self.out_fields[2].get_average_until_timestep(timestep),
                       headwidth=2, headlength=1, scale=2000, color='k', figure=fig)
        return fig, h

    def plot_instant_field(self, index, fig='new'):
        if fig is 'new':
            fig = plt.figure()
        h = plt.quiver(self.profile.get_data('y'), self.profile.get_data('z'),
                       self.out_fields[1].get_indexed_field(index),
                       self.out_fields[2].get_indexed_field(index),
                       headwidth=2, headlength=1, scale=2000, color='k', figure=fig)
        plt.title(str(index))
        return fig, h

    def plot_k_field(self, timestep=0):
        fig = plt.figure(figsize=(15, 4))
        plt.subplot(1, 2, 1)
        k_prof = self.profile.get_data('k')
        crange = (np.percentile(k_prof, 10) * 0.9, np.percentile(k_prof, 90) * 1.1)
        plt.scatter(self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                    self.profile.dataarrays[self.profile.get_datanames().index('z')].data, 10,
                    k_prof)
        plt.clim(crange[0], crange[1])
        plt.axis('equal')
        plt.subplot(1, 2, 2)
        plt.scatter(self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                    self.profile.dataarrays[self.profile.get_datanames().index('z')].data, 10,
                    self.get_k_from_outfields(timestep))
        plt.clim(crange[0], crange[1])
        plt.axis('equal')
        plt.colorbar()
        fig.suptitle('Statistisches k Filter: Zeitschritt {:d}'.format(timestep))
        return fig

    def get_diff_mean_field_profile(self, timestep=0):
        diff_fields = [Field(el.name) for el in self.out_fields]
        for i, el in enumerate(diff_fields):
            mean_field = self.out_fields[i].get_average_until_timestep(timestep)
            el.add_field(mean_field - self.profile.get_data(el.name))
        return diff_fields

    def get_mean_diff(self, timestep=0):
        diff_fields = self.get_diff_mean_field_profile(timestep)
        return np.array([np.mean(np.abs(el.field)) for el in diff_fields])

    def get_max_diff(self, timestep=0):
        diff_fields = self.get_diff_mean_field_profile(timestep)
        return np.array([np.max(np.abs(el.field)) for el in diff_fields])

    def get_k_from_outfields(self, timestep=0):
        return 0.5 * sum([el.get_var_until_timestep(timestep)
                          for el in self.out_fields if el.name in ['u', 'v', 'w']])

    def get_k_diff_profile(self, timestep=0):
        return self.get_k_from_outfields(timestep) - self.profile.get_data('k')

    def get_k_mean_diff(self, timestep=0):
        return np.mean(np.abs(self.get_k_diff_profile(timestep)))

    def get_k_max_diff(self, timestep=0):
        return np.max(np.abs(self.get_k_diff_profile(timestep)))

    def diff_stats(self, timestep=0):
        self.mean_diff = np.append(self.mean_diff,
                                   [np.append(self.get_mean_diff(timestep), self.get_k_mean_diff(timestep))],
                                   axis=0)
        self.max_diff = np.append(self.max_diff,
                                  [np.append(self.get_max_diff(timestep), self.get_k_max_diff(timestep))],
                                  axis=0)

    def plot_diff_stats(self, fig='new', noshow=False):
        if fig is 'new':
            fig = plt.figure()
        else:
            fig.clf()
        cols = ['k', 'b', 'r', 'c']
        names = ['u', 'v', 'w', 'k']
        for i in range(4):
            normval = np.max(self.profile.get_data(names[i])[np.logical_not(self.profile.outside_flag)]) - \
                      np.min(self.profile.get_data(names[i])[np.logical_not(self.profile.outside_flag)])
            if normval == 0:
                normval = self.u_bulk
            plt.plot(self.mean_diff[:, i] / normval, color=cols[i], linestyle='-', figure=fig, label=names[i] + ' mean')
            plt.plot(self.max_diff[:, i] / normval, color=cols[i], linestyle='-.', figure=fig, label=names[i] + ' max')
        plt.legend(loc='upper right')
        plt.grid(b=True)
        if not noshow:
            plt.draw()
            plt.pause(1e-17)
            time.sleep(0.00001)
        return fig

    def diff_stats_complete(self):
        self.mean_diff = np.zeros((0, 4))
        self.max_diff = np.zeros((0, 4))
        for i in range(1, self.out_fields[0].field.shape[0] + 1):
            self.diff_stats(i)
            print('Calculating convergence statistics {:d} / {:d}'.format(i, self.out_fields[0].field.shape[0] + 1))


def plot_field(y, z, field, hwid,hlen, scale, fwidth, ind, typ, cmin, cmax):
    fig = plt.figure()
    if fwidth==0:
        plt.contourf(y, z, field[0].field[-1,:,:],20,vmin=cmin,vmax=cmax)
        plt.quiver(y,z,field[1].field[-1,:,:],field[2].field[-1,:,:],
                   headwidth=hwid,headlength=hlen,scale=scale)
    else:
        plt.contourf(y, z, field[0].field[-fwidth-1, fwidth:-fwidth, fwidth:-fwidth],20,vmin=cmin,vmax=cmax)
        plt.quiver(y,z,field[1].field[-fwidth-1,fwidth:-fwidth,fwidth:-fwidth],
                   field[2].field[-fwidth - 1, fwidth:-fwidth, fwidth:-fwidth],
                   headwidth=hwid,headlength=hlen,scale=scale)
    #plt.colorbar()
    #plt.clim(cmin,cmax)
    fig.savefig('/home_elbe/goertz/Bilder/fields/{:s}_{:03d}.png'.format(typ, ind))
    plt.close(fig)


class Profil:

    def __init__(self, name, lengthscale_x=0, lengthscale_y=0, lengthscale_z=0):
        self.dataarrays = []
        self.res_y = 0
        self.res_z = 0
        self.lengthscale_x = lengthscale_x
        self.lengthscale_y = lengthscale_y
        self.lengthscale_z = lengthscale_z
        self.lnx = np.zeros(0)
        self.lny = np.zeros(0)
        self.lnz = np.zeros(0)
        self.outside_flag = 0
        self.rst_transform = []
        self.res = 0
        self.name = name

    def prepare_for_filter(self):
        pass

    def plot_quiver(self):
        y_ind = self.get_datanames().index('y')
        z_ind = self.get_datanames().index('z')
        v_ind = self.get_datanames().index('v')
        w_ind = self.get_datanames().index('w')
        plt.figure()
        plt.quiver(self.dataarrays[y_ind].data,
                   self.dataarrays[z_ind].data,
                   self.dataarrays[v_ind].data,
                   self.dataarrays[w_ind].data,
                   headwidth=2, headlength=1, scale=2000)

    def get_datanames(self):
        return [el.standardname for el in self.dataarrays]

    def get_data(self, name, rst_flag = False):
        if self.polar and name == 'v' and rst_flag:
            name = 'vtan'
        if self.polar and name == 'w' and rst_flag:
            name = 'vrad'
        return self.dataarrays[self.get_datanames().index(name)].data

    def set_outside_flag(self):
        # Mark points outside of BC area where U_i = u_ij = 0
        # --> To save computation time
        u_ind = self.get_datanames().index('u')
        v_ind = self.get_datanames().index('v')
        w_ind = self.get_datanames().index('w')

        self.outside_flag = (abs(self.dataarrays[u_ind].data) +
                             abs(self.dataarrays[v_ind].data) +
                             abs(self.dataarrays[w_ind].data)) == 0

    def init_length_fields(self):
        self.lnx = math.ceil(self.lengthscale_x / self.res) * np.ones((self.res_z, self.res_y))
        if self.lengthscale_y==0:
            self.lny = self.lnx
        else:
            self.lny = math.ceil(self.lengthscale_y / self.res) * np.ones((self.res_z, self.res_y))
        if self.lengthscale_z==0:
            self.lnz = self.lnx
        else:
            self.lnz = math.ceil(self.lengthscale_z / self.res) * np.ones((self.res_z, self.res_y))

    def plot_input_profile(self, outdir):
        if self.polar:
            # for polar plot 3 sectors
            phi = np.concatenate((self.dataarrays[self.get_datanames().index('phi')].data[:,:91]-22.5*np.pi/180, self.dataarrays[self.get_datanames().index('phi')].data[:,:91], self.dataarrays[self.get_datanames().index('phi')].data[:,:91]+22.5*np.pi/180), axis=1)
            r = np.concatenate((self.dataarrays[self.get_datanames().index('r')].data[:,:91], self.dataarrays[self.get_datanames().index('r')].data[:,:91], self.dataarrays[self.get_datanames().index('r')].data[:,:91]), axis=1)
            y = r * -np.sin(phi)
            z = r * np.cos(phi)

            for el in self.dataarrays:
                if el.standardname not in ['x', 'y', 'z']:
                    plot_id = el.standardname
                    nplt.contourf2(plot_id, el.unit,
                                   y,
                                   z,
                                   np.concatenate((el.data[:,:91], el.data[:,:91], el.data[:,:91]), axis = 1), 100, 'y', 'z',
                                   plot_id, outdir + '/' + plot_id, figsize=(16, 16 * ((self.res_z-3) / (self.res_y-3)/3)))
                    nplt.close(plot_id)
        else:
        # for cartesian plot 1 sector
            for el in self.dataarrays:
                if el.standardname not in ['x', 'y', 'z']:
                    plot_id = el.standardname
                    nplt.contourf2(plot_id, el.unit,
                                   self.get_data('y'),
                                   self.get_data('z'),
                                   el.data, 100, 'y', 'z',
                                   plot_id, outdir + '/' + plot_id, figsize=(8, 8 * self.res_z / self.res_y))
                    nplt.close(plot_id)

class Profil1DProf(Profil):
# build fields for U, V, W, (Tt, pt) from 1D profile or from scratch

    def __init__(self, res, profile='hyperbolic-tangent', turb_prof='top-hat', u_bulk=0, tu_intensity=0,
                 filename='none', res_y=10, res_z=10, oy=0, oz=0, lengthscale=1, y_range=(0, 0), z_range=(0, 0),
                 inner_diameter_norm=0, prof1d=False):

        ################################################################################################################
        ####################################### Setting Parameters #####################################################
        ################################################################################################################

        # define profile name from filename and velocity profile. Set to "none" if no file is defined (prof1d = False)
        if prof1d:
            name = filename + ' ' + profile
        else:
            name = profile
            filename = 'none'
        Profil.__init__(self, name)

        # bulk velocity profile shape
        self.profile = profile

        # geometric resolution in m
        self.res = res

        # File name of input file
        self.filename = filename

        # polar coordinates not possible yet
        self.polar = False

        # define geometric range of inlet field
        if not y_range == (0, 0):
            self.y_range = y_range
        else:
            self.y_range = (oy - res_y / 2 * res, oy + res_y / 2 * res)

        if not z_range == (0, 0):
            self.z_range = z_range
        else:
            self.z_range = (oz - res_z / 2 * res, oz + res_z / 2 * res)

        # define lenght scale
        self.lengthscale_x = lengthscale * self.res

        # shape of turbulent profile
        self.turb_prof = turb_prof

        # bulk velocity in m/s
        self.u_bulk = u_bulk

        # turbulent intensity
        self.tu_intensity = tu_intensity

        # coordinate origin
        self.oy = np.mean(self.y_range)
        self.oz = np.mean(self.z_range)

        # no file reader
        self.fr = None

        #
        self.inner_diameter_norm = inner_diameter_norm

        ################################################################################################################
        ####################################### Create profiles ########################################################
        ################################################################################################################

        self.prof_functions = self.make_prof_functions()

    def prepare_for_filter(self):
        self.build_coord_grid()
        self.build_profile()

        self.set_outside_flag()

        self.init_length_fields()

        self.rst_transform = RstTransform(self.dataarrays, self.res_y, self.res_z)

    def build_coord_grid(self):
        yi = np.arange(self.y_range[0], self.y_range[1], self.res)
        zi = np.arange(self.z_range[0], self.z_range[1], self.res)

        y, z = np.meshgrid(yi, zi)

        self.dataarrays.append(preptrav.DataColumn('x', y * 0))
        self.dataarrays.append(preptrav.DataColumn('y', y))
        self.dataarrays.append(preptrav.DataColumn('z', z))

        self.res_y = len(yi)
        self.res_z = len(zi)

    def make_prof_functions(self):
        prof_functions = []

        # build profiles for u,v,w and RST
        if self.filename == 'none':

            # from bulk velocity and turbulence intensity if no 1D profile is defined
            f = lambda x: (1 + np.tanh(10 * (-np.abs(x) + 0.5))) / 2
            prof_functions.append(ProfFunc('u', lambda u, y, z: u * f(y) * f(z)))
            prof_functions.append(ProfFunc('v', lambda u, y, z: 0 * y))
            prof_functions.append(ProfFunc('w', lambda u, y, z: 0 * y))
            prof_functions.append(ProfFunc('uu', lambda u, y, z: u ** 2 * self.tu_intensity ** 2 *
                                                                 f(y) ** 2 * f(z) ** 2))
            prof_functions.append(ProfFunc('vv', lambda u, y, z: u ** 2 * self.tu_intensity ** 2 *
                                                                 f(y) ** 2 * f(z) ** 2))
            prof_functions.append(ProfFunc('ww', lambda u, y, z: u ** 2 * self.tu_intensity ** 2 *
                                                                 f(y) ** 2 * f(z) ** 2))
            prof_functions.append(ProfFunc('uv', lambda u, y, z: 0 * y))
            prof_functions.append(ProfFunc('uw', lambda u, y, z: 0 * y))
            prof_functions.append(ProfFunc('vw', lambda u, y, z: 0 * y))
        else:

            # read 1D profile if exists
            # create file reader
            fr = preptrav.FileReader(self.filename)

            # open file and read data from file
            fr.open_file()

            # write data in data_columns array
            fr.extract_column_data()
            self.fr = fr

            # prepare y data
            y_data = fr.get_data('y')
            y_data = (y_data - min(y_data)) / (max(y_data) - min(y_data))
            y_data = np.abs(y_data / 2 - 0.5)

            # define field names
            fieldnames = ['u', 'v', 'w', 'uu', 'vv', 'ww', 'uv', 'uw', 'vw']

            # creation of profiles for all elements in fieldnames
            for el in fr.data_columns:
                if el.standardname in fieldnames:
                    val_data = el.data

                    # replace 0 values
                    if val_data[-1] == 0:
                        val_data[-1] = -1e-5
                    ref_val = val_data[-1]

                    # normalise values
                    val_data = val_data / ref_val

                    # create profile funktions
                    prof_functions.append(ProfFunc(el.standardname,
                                                   interp=interpolate.interp1d(y_data, val_data,
                                                                               fill_value='extrapolate'),
                                                   ref_val=ref_val))

            # create profile funktions vor variables not listed in fieldnames
            colnames = [el.name for el in prof_functions]
            for name in fieldnames:
                if name not in colnames:
                    prof_functions.append(ProfFunc(name, lambda u, y, z: 0 * y))

        return prof_functions

    def build_profile(self):
        y_data = self.get_data('y')
        y_data = y_data / (y_data.max() - y_data.min())
        y_data = y_data - y_data.min() - 0.5
        z_data = self.get_data('z')
        z_data = z_data / (z_data.max() - z_data.min())
        z_data = z_data - z_data.min() - 0.5

        if self.profile == 'hyperbolic-tangent':
            z_data = 0
            outside = y_data * 0
        elif self.profile == 'double-hyperbolic-tangent':
            outside = y_data * 0
        elif self.profile == 'circular-hyperbolic-tangent':
            y_data = np.sqrt(y_data ** 2 + z_data ** 2)
            outside = y_data > 0.5
            z_data = 0
        elif self.profile == 'ring-hyperbolic-tangent':
            mid_radius = (1 + self.inner_diameter_norm) / 4
            d = (1 - self.inner_diameter_norm) / 2
            radius = np.sqrt(y_data ** 2 + z_data ** 2)
            y_data = (radius - mid_radius) / d
            outside = np.abs(y_data) > 0.5
            z_data = 0
        else:
            outside = y_data * 0
        outside = outside.astype('bool')
        for func in self.prof_functions:
            self.dataarrays.append(preptrav.DataColumn(func.name,
                                                       func.get_func(self.profile)(self.u_bulk, y_data, z_data)))
            self.dataarrays[-1].data[outside] = 0
        self.dataarrays.append(preptrav.DataColumn('k', 0.5 * (self.get_data('uu') +
                                                               self.get_data('vv') +
                                                               self.get_data('ww'))))


class ProfFunc:
# Object ProfFunktion creates pfrofiles
    def __init__(self, name, func=None, interp=None, ref_val=None):
        self.name = name
        self.interp = interp
        self.ref_val = ref_val
        if func is not None:
            self.func_1d = func
            self.func_2d = func
        else:
            self.func_2d = lambda u, y, z: self.ref_val * self.interp(np.abs(y)) * self.interp(np.abs(z))
            self.func_1d = lambda u, y, z: self.ref_val * self.interp(np.abs(y))

    def get_func(self, profile=None):
        if profile == 'double-hyperbolic-tangent':
            return self.func_2d
        else:
            return self.func_1d


class Profil2D(Profil):
# preperation of 2D field data: extract mean fields of U, V, W, (Tt, pt), k and epsilon/omega from data file

    def __init__(self, filename, res, mdot, bulk_velocity, density, polar, fluct_comps=('u', 'v', 'w'),
                 temp_fluct=0, lengthscale_x=0, lengthscale_y=0, lengthscale_z=0):

        ################################################################################################################
        ####################################### Setting Parameters #####################################################
        ################################################################################################################

        Profil.__init__(self, filename, lengthscale_x, lengthscale_y, lengthscale_z)
        # file name of input data
        self.filename = filename
        # grid resulution in m
        self.res = res
        # massflow in kg/s for scaling
        self.mdot = mdot
        # bulk velocity for scaling
        self.bulk_velocity = bulk_velocity
        # density for massflow scaling
        self.den = density
        self.valid = False
        self.mode = 'none'
        self.rad_flag = 0
        self.gradients = []
        self.plane = None
        # components fo fluctuating fields
        self.fluct_comps = fluct_comps
        # standard deviation for temperature fluctuations
        self.temp_fluct = temp_fluct
        # for polar coordinates = True, cartesian = False
        self.polar = polar

        ################################################################################################################
        ####################################### Import from File #######################################################
        ################################################################################################################
        # Create object FileReader
        self.fr = preptrav.FileReader(self.filename)

        # open file and read data from file
        self.fr.open_file()

        # write data in data_columns array
        self.fr.extract_column_data()

        # for polar coordinates, calculate polar coordinates and velocity components
        if self.polar:
            self.fr.make_pol_coord()
            self.fr.make_pol_velo()


    def prepare_for_filter(self):
        # prepare data for filter
        self.check_validity()

        ################################# Turbulence Settings ##########################################################
        # calculate epsilon from omega if necessary
        if self.mode == '2-equation' and 'e' not in self.fr.get_standardnames():
            self.convert_sdr_to_eps()

        # use RST from data if available
        elif self.mode == 'reynolds-stresses':
            self.check_reynolds_fields()

        # initialise correlations with temperature if temperature fluctuations are active
        if 'T' in self.fluct_comps:
            name_tcorr = ['uT', 'vT', 'wT', 'TT']
            self.check_temperature_fluct(name_tcorr)
        elif 'Tt' in self.fluct_comps:
            name_tcorr = ['uTt', 'vTt', 'wTt', 'TtTt']
            self.check_temperature_fluct(name_tcorr)

        ############################# Interpolation of data to new grid ################################################
        self.interpolate_data()

        ################################################### Scaling#####################################################
        # massflow or bulkvelocity scaling
        if self.mdot != 0.0 or self.bulk_velocity != 0.0:
            self.scale_vel()

        self.set_outside_flag()

        ############################################### Calculation of RST #############################################
        if self.mode == '2-equation':
            self.calculate_rst_from_k_eps()

        # calculation of length scale for filter
        self.calculate_lengthscale()

        ############################################# Calculation of RST transformation (Cholesky Matrix) ##############
        if self.polar:
            self.rst_transform = RstTransform(self.dataarrays, self.res_phi, self.res_r, fids=self.fluct_comps)
        else:
            self.rst_transform = RstTransform(self.dataarrays, self.res_y, self.res_z, fids=self.fluct_comps)

        if 'k' not in self.get_datanames():
            self.dataarrays.append(preptrav.DataColumn('k', 0.5 * (self.get_data('uu') + self.get_data('vv') +
                                                              self.get_data('ww'))))

    def check_validity(self):
        # check if coordinates and velocity exist in the profile
        fields = self.fr.get_standardnames()
        coord_check = all([val in fields for val in ['x', 'y', 'z']])
        velo_check = all([val in fields for val in ['u', 'v', 'w']])

        # check which turbulent data exist in profile
        keps_check = all(['k' in fields] +
                         [any([val in fields for val in ['e', 'sdr']])])
        rs_check1 = all([val in fields for val in ['uu', 'vv', 'ww']])
        turb_check = any([keps_check, rs_check1])

        self.valid = all([coord_check, velo_check, turb_check])

        # set mode depending on existing turbulent values
        if rs_check1:
            self.mode = 'reynolds-stresses'
        elif turb_check:
            self.mode = '2-equation'

    def interpolate_data(self):
        # caching of coordinates
        x_data = self.fr.data_columns[self.fr.get_standardnames().index('x')].data
        y_data = self.fr.data_columns[self.fr.get_standardnames().index('y')].data
        z_data = self.fr.data_columns[self.fr.get_standardnames().index('z')].data

        if self.polar:
            ############################################################################################################
            ############################# polar coordinates ############################################################
            ############################################################################################################

            # expand data columns to three sectors
            self.fr.expand_rot_periodic(-5.625)

            # get polar coordinates
            r_data = self.fr.data_columns[self.fr.get_standardnames().index('radius')].data
            phi_data = self.fr.data_columns[self.fr.get_standardnames().index('phi_rad')].data

            # calculate phi from original data (one sector)
            phi_data_original = np.arctan(-y_data / z_data)

            # get x coordinate = constant
            xval = np.mean(x_data)

            # calculate channel high from r field
            channel_hight = np.max(r_data)-np.min(r_data)

            ############################################ r coordinate ##################################################
            # calculate resolution in r direction
            l = channel_hight / self.res
            lr = round(l)
            self.disc_r = channel_hight / lr

            # maximal and minimal radius
            rmax = np.max(r_data)
            rmin = np.min(r_data)

            # define new r vector from min to max with resolution disc_r
            ri = np.arange(np.min(r_data), np.max(r_data), self.disc_r)
            # ri =np.append(ri, np.max(r_data))

            Rmax = np.max(ri)
            Rmin = np.min(ri)

            ############################################## phi coordinate ##############################################
            # get maximum angle phi_max
            left_phi = max(phi_data_original)#+0.03

            # calculate min angle from sector width or read from data
            # right_phi = min(phi_data_original)-0.03
            right_phi = left_phi - 22.5*np.pi/180

            # check sector
            sector = left_phi - right_phi
            sector_deg = sector * 180 / np.pi

            # calculate resolution in phi direction at mid channel
            rm = np.mean(ri)
            res_phi1 = np.arcsin(self.res / rm)

            # calculate resolution in phi direction
            q = sector / res_phi1
            qr = round(q)
            self.disc_phi = sector / qr

            # define new phi vector from min to max with resolution disc_phi
            phii = np.arange(right_phi, left_phi, self.disc_phi)
            phii = np.append(phii, left_phi)

            ############################################ construct new mesh in r and phi ###############################
            phi, r = np.meshgrid(phii, ri)

            # save resolution
            self.res_r = len(ri)
            self.res_phi = len(phii)
            self.res_y = len(phii)
            self.res_z = len(ri)

            ########################################## interpolation of data onto new grid #############################
            # interpolation on to grid (with NaN on outer boundaries) using linear method
            self.dataarrays = [preptrav.DataColumn(el.standardname, interpolate.griddata(
                (r_data, phi_data), el.data, (r, phi), method='linear'))
                               for el in self.fr.data_columns
                               if el.standardname not in ['x', 'y', 'z', 'radius', 'phi_deg', 'phi_rad']]

            # replace the most inner and outer line with the mean value of the original data on hub and shroud
            for el in self.fr.data_columns:
                if el.standardname not in ['x', 'y', 'z', 'radius', 'phi_deg', 'phi_rad']:
                    self.dataarrays[self.get_datanames().index(el.standardname)].mean_inner = \
                        np.mean(self.fr.data_columns[self.fr.get_standardnames().index(el.standardname)].data
                                [np.round(r_data, 6) == np.round(rmin, 6)])

                    self.dataarrays[self.get_datanames().index(el.standardname)].data[np.round(r, 6) ==
                                                                                      np.round(rmin, 6)] = \
                        self.dataarrays[self.get_datanames().index(el.standardname)].mean_inner

                    self.dataarrays[self.get_datanames().index(el.standardname)].mean_outer = \
                        np.mean(self.fr.data_columns[self.fr.get_standardnames().index(el.standardname)].data
                                [np.round(r_data, 6) == np.round(rmax, 6)])

                    self.dataarrays[self.get_datanames().index(el.standardname)].data[np.round(r, 6) ==
                                                                                      np.round(rmax, 6)] =\
                        self.dataarrays[self.get_datanames().index(el.standardname)].mean_outer

            # calculate cartesion coordinates of new grid
            y = -r * np.sin(phi)
            z = r*np.cos(phi)
            x = y * 0 + xval

            # append coordinates to data array
            self.dataarrays.append(preptrav.DataColumn('x', x))
            self.dataarrays.append(preptrav.DataColumn('y', y))
            self.dataarrays.append(preptrav.DataColumn('z', z))
            self.dataarrays.append(preptrav.DataColumn('r', r))
            self.dataarrays.append(preptrav.DataColumn('phi', phi))

        else:
            ############################################################################################################
            ################################### cartesian coordinates ##################################################
            ############################################################################################################
            # create plane object --> includes transform on y-z Plane
            self.plane = Plane(x_data, y_data, z_data)
            # cache coordinates
            x_data = self.plane.x_data_rot
            y_data = self.plane.y_data_rot
            z_data = self.plane.z_data_rot

            # calculation of x-coordinate (constant)
            xval = np.mean(x_data)

            # determination of inner radius of inlet sector, setting velocity to zero inside this radius
            inner_radius = min(np.sqrt(y_data ** 2 + z_data ** 2))

            # determination of value range of y and z, creation of cartesian coordinate vectors
            yi = np.arange(np.min(y_data), np.max(y_data), self.res)
            zi = np.arange(np.min(z_data), np.max(z_data), self.res)

            ########################################## creation of grid ################################################
            y, z = np.meshgrid(yi, zi)

            # matrix form:
            #     | y_min ... y_max |        | z_min ... z_min |
            # y = |                 | ,  z = |                 |
            #     | y_min ... y_max |        | z_max ... z_max |
            # y(i, j) : 1. Index row --> varying z, 2. Index collumn --> varying y

            # size of matrix
            self.res_y = len(yi)
            self.res_z = len(zi)

            ############################## interpolation of imported fields on new grid ################################
            self.dataarrays = [preptrav.DataColumn(el.standardname, interpolate.griddata(
            (y_data, z_data), el.data, (y, z), fill_value=0, method='linear'))
                           for el in self.fr.data_columns
                           if el.standardname not in ['x', 'y', 'z']]

            # determination of values inside the inner radius
            self.rad_flag = np.sqrt(y ** 2 + z ** 2) < inner_radius

            # setting of values inside the inner radius to 0
            for el in self.dataarrays:
                el.data[self.rad_flag] = 0

            # rotate interpolated plane back on original location
            x, y, z = self.plane.rotate_back(y * 0 + xval, y, z)

            # append coordinates to data array
            self.dataarrays.append(preptrav.DataColumn('x', x))
            self.dataarrays.append(preptrav.DataColumn('y', y))
            self.dataarrays.append(preptrav.DataColumn('z', z))

    def scale_vel(self):
        # scaling of velocity and turbulent values to target mass flow or bulk velocity
        xn = 1  # get out of coordinate transformation
        yn = 0
        zn = 0

        # calc old mdot
        c_area = self.res ** 2
        A = c_area * (self.res_y - 1) * (self.res_z - 1)
        meanVel = [np.mean(self.dataarrays[self.get_datanames().index(u_id)].data)
                   for u_id in ['u', 'v', 'w']]
        udotn = meanVel[0] * xn + meanVel[1] * yn + meanVel[2] * zn
        mdot_old = udotn * A * self.den

        # calculate old turbulent intensity
        TI = np.sqrt(2. / 3. * self.dataarrays[self.get_datanames().index('k')].data) / \
             np.sqrt(self.dataarrays[self.get_datanames().index('u')].data ** 2 + \
                     self.dataarrays[self.get_datanames().index('v')].data ** 2 + \
                     self.dataarrays[self.get_datanames().index('w')].data ** 2)

        # find old length scales
        L = self.dataarrays[self.get_datanames().index('k')].data ** (3 / 2) / \
            self.dataarrays[self.get_datanames().index('e')].data

        # determine scaling factor for new velocities
        if self.mdot != 1 and self.bulk_velocity != 1:
            print('Scaling of massflow and bulk velocity at the same time not possible!')
            return
        elif self.mdot != 1:
            # from new massflow
            scale = self.mdot / mdot_old
            print('Scaling to Massflow: ', self.mdot, 'kg/s')
        elif self.bulk_velocity != 1:
            # from new bulk velocity
            scale = self.bulk_velocity / udotn
            print('Scale to bulk velocity: ', self.bulk_velocity, 'm/s')

        # calculate new velocities
        U = self.dataarrays[self.get_datanames().index('u')].data * scale
        V = self.dataarrays[self.get_datanames().index('v')].data * scale
        W = self.dataarrays[self.get_datanames().index('w')].data * scale

        # find new k
        k = TI ** 2 * (U ** 2 + W ** 2 + V ** 2)

        # find new epsilon
        e = k ** (3 / 2) / L

        # write into data array
        self.dataarrays[self.get_datanames().index('u')].data = U
        self.dataarrays[self.get_datanames().index('v')].data = V
        self.dataarrays[self.get_datanames().index('w')].data = W
        self.dataarrays[self.get_datanames().index('k')].data = np.nan_to_num(k)
        self.dataarrays[self.get_datanames().index('e')].data = np.nan_to_num(e)

    def calculate_rst_from_k_eps(self):
        ################################################################################################################
        #################################### polar coordinates #########################################################
        ################################################################################################################

        if self.polar:

            ########################### calculation of gradients in polar coordinates ##################################

            dVaxdr = np.gradient(self.dataarrays[self.get_datanames().index('u')].data, self.disc_r, axis=0)
            dVaxdphi = np.gradient(self.dataarrays[self.get_datanames().index('u')].data, self.disc_phi, axis=1) /\
                       self.dataarrays[self.get_datanames().index('r')].data
            dVraddr = np.gradient(self.dataarrays[self.get_datanames().index('vrad')].data, self.disc_r, axis=0)
            dVraddphi = np.gradient(self.dataarrays[self.get_datanames().index('vrad')].data, self.disc_phi, axis=1) /\
                        self.dataarrays[self.get_datanames().index('r')].data
            dVtandr = np.gradient(self.dataarrays[self.get_datanames().index('vtan')].data, self.disc_r, axis=0)
            dVtandphi = np.gradient(self.dataarrays[self.get_datanames().index('vtan')].data, self.disc_phi, axis=1) /\
                        self.dataarrays[self.get_datanames().index('r')].data

            # approximation dVaxdx from incompressibility
            dVaxdx = - dVraddr - dVtandphi

            # set other gradients in flow direction to zero
            dVraddx = np.zeros((self.res_r, self.res_phi), dtype=np.float64)
            dVtandx = np.zeros((self.res_r, self.res_phi), dtype=np.float64)

            # write gradients to array
            self.gradients = [[dVaxdx, dVaxdphi, dVaxdr],
                          [dVtandx, dVtandphi, dVtandr],
                          [dVraddx, dVraddphi, dVraddr]]

            # append gradients to data array
            name_rst = [['dVaxdx', 'dVaxdphi', 'dVaxdr'], ['dVtandx', 'dVtandphi', 'dVtandr'],
                        ['dVraddx', 'dVraddphi', 'dVraddr']]
            for i in range(3):
                for j in range(3):
                    self.dataarrays.append(preptrav.DataColumn(name_rst[i][j], self.gradients[i][j]))

            #################################### eddy viscosity model ##################################################

            k = self.dataarrays[self.get_datanames().index('k')].data
            eps = self.dataarrays[self.get_datanames().index('e')].data
            flag = np.where(eps > 0)
            nu_t = np.zeros((self.res_r, self.res_phi), dtype=np.float64)

            # calculate turbulent viscosity
            nu_t[flag] = 0.09 * k[flag] ** 2 / eps[flag]

            # calculate Reynolds Stresses with Boussinesq hypothesis and append to data array
            self.dataarrays.append(preptrav.DataColumn('uu', -2 * nu_t * dVaxdx + 2 / 3 * k))
            self.dataarrays.append(preptrav.DataColumn('vv', -2 * nu_t * dVtandphi + 2 / 3 * k))
            self.dataarrays.append(preptrav.DataColumn('ww', -2 * nu_t * dVraddr + 2 / 3 * k))
            self.dataarrays.append(preptrav.DataColumn('uv', -nu_t * (dVaxdphi + dVtandx)))
            self.dataarrays.append(preptrav.DataColumn('uw', -nu_t * (dVaxdr + dVraddx)))
            self.dataarrays.append(preptrav.DataColumn('vw', -nu_t * (dVtandr + dVraddphi)))

        ################################################################################################################
        ########################################### cartesian coordinates ##############################################
        ################################################################################################################

        else:
            ########################################## gradients in cartesian coordinates ##############################
            # format: gradients = [[dudz, dudy], [dvdz,dvdy], [dwdz,dwdy]]
            # --> Change z-direction to axis 0, change y-direction to axis 1
            gradients = [npc.gradient(self.dataarrays[self.get_datanames().index(u_id)].data, self.res, self.outside_flag)
                      for u_id in ['u', 'v', 'w']]

            # for compatibility with old script
            dudy = gradients[0][1]
            dudz = gradients[0][0]
            dvdy = gradients[1][1]
            dvdz = gradients[1][0]
            dwdy = gradients[2][1]
            dwdz = gradients[2][0]

            # approximation dudx from incompressibility
            dudx = - dvdy - dwdz

            # setting all other gradients in flow direction to 0
            dvdx = np.zeros((self.res_z, self.res_y), dtype=np.float64)
            dwdx = np.zeros((self.res_z, self.res_y), dtype=np.float64)

            # caching gradients in object
            self.gradients = [[dudx, dudy, dudz],
                            [dvdx, dvdy, dvdz],
                            [dwdx, dwdy, dwdz]]

            # write gradient fields into dataarrays
            name_rst = [['dudx', 'dudy', 'dudz'], ['dvdx', 'dvdy', 'dvdz'], ['dwdx', 'dwdy', 'dwdz']]
            for i in range(3):
                for j in range(3):
                    self.dataarrays.append(preptrav.DataColumn(name_rst[i][j], self.gradients[i][j]))

            #################################### eddy viscosity model ##################################################

            k = self.dataarrays[self.get_datanames().index('k')].data
            eps = self.dataarrays[self.get_datanames().index('e')].data
            flag = np.where(eps > 0)
            nu_t = np.zeros((self.res_z, self.res_y), dtype=np.float64)

            # calculate turbulent viscosity
            nu_t[flag] = 0.09 * k[flag] ** 2 / eps[flag]

            # calculate Reynolds Stresses with Boussinesq hypothesis and append to data array
            self.dataarrays.append(preptrav.DataColumn('uu', -2 * nu_t * dudx + 2 / 3 * k))
            self.dataarrays.append(preptrav.DataColumn('vv', -2 * nu_t * dvdy + 2 / 3 * k))
            self.dataarrays.append(preptrav.DataColumn('ww', -2 * nu_t * dwdz + 2 / 3 * k))
            self.dataarrays.append(preptrav.DataColumn('uv', -nu_t * (dudy + dvdx)))
            self.dataarrays.append(preptrav.DataColumn('uw', -nu_t * (dudz + dwdx)))
            self.dataarrays.append(preptrav.DataColumn('vw', -nu_t * (dvdz + dwdy)))

            # values inside the inner radius to 0
            for el in self.dataarrays:
                if el.standardname in ['uu', 'vv', 'ww', 'uv', 'uw', 'vw']:
                    el.data[self.rad_flag] = 0

            # setting negative stresses to 0
            for el in self.dataarrays:
                if el.standardname in ['uu', 'vv', 'ww']:
                    el.data[np.where(el.data < 0)] = 0

    def convert_sdr_to_eps(self):
        # conversion of omega into epsilon
        self.fr.insert_column('e', 0.09 * self.fr.data_columns[self.fr.get_standardnames().index('k')].data *
                              self.fr.data_columns[self.fr.get_standardnames().index('sdr')].data)

    def check_reynolds_fields(self):
        # check if shear stresses exist in data
        for name in ['uv', 'uw', 'vw']:
            if name not in self.fr.get_standardnames():
                self.fr.insert_column(name, self.fr.data_columns[0].data * 0)
                print('Reynoldsspannung ' + name + ' nicht gefunden -> Mit 0 ersetzt')

    def check_temperature_fluct(self, name_tcor):
        # set temperature fluctuation according to standard deviation
        for name in name_tcor:
            if name not in self.fr.get_standardnames():
                if 'TtTt' in name or 'TT' in name:
                    self.fr.insert_column(name, self.fr.data_columns[0].data * 0 +self.temp_fluct **2)
                    print('Temperature correlation {:s} not found. Replaced by standard deviation '.format(name)+str(self.temp_fluct)+'K')
                else:
                    self.fr.insert_column(name, self.fr.data_columns[0].data * 0)
                    print('Temperature correlation {:s} not found. Replaced by 0'.format(name))

    def calculate_lengthscale(self):
        # determining length scale for filter
        # at the moment: constant length scale for whole plane and all components, isotropic
        # output: matrix, length scale could be defined individually for every grid point
        if self.lengthscale_x == 0:
            b = np.max(self.fr.data_columns[self.fr.get_standardnames().index('y')].data) - \
                np.min(self.fr.data_columns[self.fr.get_standardnames().index('y')].data)
            c = np.max(self.fr.data_columns[self.fr.get_standardnames().index('z')].data) - \
                np.min(self.fr.data_columns[self.fr.get_standardnames().index('z')].data)
            self.lengthscale_x = 0.07 * 2 * b * c / (b + c)
        self.init_length_fields()


class Plane:

    def __init__(self, x_data, y_data, z_data):
        self.x_data = x_data
        self.y_data = y_data
        self.z_data = z_data
        self.axis = None
        self.angle = None
        self.norm = self.get_normal()
        self.center = self.get_centre_point()
        self.x_data_rot, self.y_data_rot, self.z_data_rot = self.rotate()

    def get_normal(self):
        # calculate vectors in plane
        # maximum and minimum values for each coordiante --> indizes to build 3 vectors for each value pair
        x_inds = (self.x_data.argmin(),
                  self.x_data.argmax())
        y_inds = (self.y_data.argmin(),
                  self.y_data.argmax())
        z_inds = (self.z_data.argmin(),
                  self.z_data.argmax())

        # build vectors between minimum and maximum points (in plane)
        vectors = [np.array([self.x_data[ind[1]] - self.x_data[ind[0]],
                             self.y_data[ind[1]] - self.y_data[ind[0]],
                             self.z_data[ind[1]] - self.z_data[ind[0]]])
                   for ind in [x_inds, y_inds, z_inds]]

        norms = []

        # form cross product for each pair of vectors (1 or 3 products possible, skip if (0,0,0) vector)
        for i in range(3):
            for j in range(i):
                if i != j and not all(vectors[i] == np.array([0, 0, 0])) and not all(vectors[j] == np.array([0, 0, 0])):
                    n = np.cross(vectors[i], vectors[j])
                    norms.append(np.round(n / np.linalg.norm(n), 5))

        # return first vector in list that has length 1
        return norms[[np.round(np.linalg.norm(v)) for v in norms].index(1)]

    def get_centre_point(self):
        xc = (np.amax(self.x_data) + np.amin(self.x_data)) / 2
        yc = (np.amax(self.y_data) + np.amin(self.y_data)) / 2
        zc = (np.amax(self.z_data) + np.amin(self.z_data)) / 2
        return [xc, yc, zc]

    def rotate(self):
        if all(np.cross(np.array([1, 0, 0]), self.norm) == 0):
            self.angle = 0
            return self.x_data, self.y_data, self.z_data
        else:
            # rotation with pointer (axis, angle)
            axis = npc.unit_vector(np.cross(np.array([1, 0, 0]), self.norm))
            angle = npc.angle_between(np.array([1, 0, 0]), self.norm)
            x_rot, y_rot, z_rot = rotate(self.x_data, self.y_data, self.z_data, axis, angle)
            self.axis = axis
            self.angle = angle
            return x_rot, y_rot, z_rot

    def rotate_back(self, x, y, z):
        if self.angle == 0:
            return x, y, z
        else:
            return rotate(x, y, z, self.axis, -self.angle)


def rotate(x_data, y_data, z_data, axis, angle):
    ax_mat = np.array([[0, -axis[2], axis[1]],
                       [axis[2], 0, -axis[0]],
                       [-axis[1], axis[0], 0]])
    rot_mat = np.eye(3) + np.sin(-angle) * ax_mat + (1 - np.cos(-angle)) * np.matmul(ax_mat, ax_mat)
    x_rot = x_data * 0
    y_rot = y_data * 0
    z_rot = z_data * 0
    for index in np.ndindex(x_data.shape):
        point = np.matmul(rot_mat, np.array([[x_data[index]],
                                             [y_data[index]],
                                             [z_data[index]]]))
        x_rot[index] = point[0]
        y_rot[index] = point[1]
        z_rot[index] = point[2]
    return x_rot, y_rot, z_rot


def calccoeff(n, ln):
    # determination of filter coefficients for one direction (from length scale and filter width)
    a = np.zeros(n * 2 + 1)

    for i in range(n * 2 + 1):
        k = i - n
        a[i] = np.exp(-np.pi * k ** 2 / (2 * ln ** 2))

    return a / np.sqrt(np.sum(a ** 2))


class FilterMatrix:
# object filter matrix
    def __init__(self, nfx, nfy, nfz, lnx, lny, lnz):
        self.nfx = int(nfx)
        self.nfy = int(nfy)
        self.nfz = int(nfz)
        self.lnx = int(lnx)
        self.lny = int(lny)
        self.lnz = int(lnz)
        # format: x_dim, z_dim, y_dim --> s. format profile
        # indexing: self.matrix[x_ind, z_ind, y_ind]
        self.matrix = np.zeros((self.nfx * 2 + 1,
                                self.nfz * 2 + 1,
                                self.nfy * 2 + 1))
        self.make_matrix()

    def make_matrix(self):
        bx = calccoeff(self.nfx, self.lnx)
        by = calccoeff(self.nfy, self.lny)
        bz = calccoeff(self.nfz, self.lnz)

        # 3D filter matrix
        for i in range(self.nfx * 2 + 1):
            for j in range(self.nfz * 2 + 1):
                for k in range(self.nfy * 2 + 1):
                    self.matrix[i, j, k] = bx[i] * bz[j] * by[k]


class Field:

    def __init__(self, name):
        self.name = name
        self.field = np.zeros((0, 0, 0))

    def add_field(self, field):
        self.field = field

    def get_last_field(self):
        return self.field[-1, :, :]

    def get_indexed_field(self, index):
        return self.field[index, :, :]

    def get_time_average(self):
        return np.mean(self.field, axis=0)

    def get_time_fluctuation(self):
        return self.field - self.get_time_average()

    def get_average_until_timestep(self, index=0):
        if index == 0:
            index = self.field.shape[0]
        return np.mean(self.field[:index, :, :], axis=0)

    def get_time_fluctuation_until_timestep(self, index=0):
        if index == 0:
            index = self.field.shape[0]
        return self.field[:index, :, :] - self.get_average_until_timestep(index)

    def get_var_until_timestep(self, index=0):
        return np.mean(self.get_time_fluctuation_until_timestep(index) ** 2, axis=0)

    def get_average(self, outside='none'):
        if outside is 'none':
            outside = self.field * 0
        return np.mean(self.field[np.logical_not(outside)])

    def plot_field(self, folder):
        pass


class RandomField(Field):

    def __init__(self, name, nfx, nfy, nfz, pdfr, polar, keep=False):
        Field.__init__(self, name)
        self.nfx = nfx
        self.nfy = nfy
        self.nfz = nfz
        self.pdfr = pdfr
        self.polar = polar
        self.extendx = int(self.nfx.max())
        self.extendy = int(self.nfy.max())
        self.extendz = int(self.nfz.max())
        self.res_z = self.nfx.shape[0]
        self.res_y = self.nfy.shape[1]

        if self.polar:
        # for polar coordinates: create random fields for 3 sectors + filter width overhang
            # field 1 sector
            field = np.random.uniform(low=-self.pdfr, high=self.pdfr,
                                       size=(self.extendx * 2 + 1,
                                             self.extendz * 2 + self.res_z,
                                             self.res_y))

            # overhang of 1 filter width
            overhang = np.random.uniform(low=-self.pdfr, high=self.pdfr,
                                       size=(self.extendx * 2 + 1,
                                             self.extendz * 2 + self.res_z,
                                             self.extendy *2))

            # allocation of 3 sec + overhang field
            self.field = np.zeros((self.extendx * 2 + 1, self.extendz * 2 + self.res_z, self.extendy * 2 +  (self.res_y)*3))

            # composition of field
            for i in range(np.shape(self.field)[0]):
                self.field[i] = np.concatenate([field[i], field[i], field[i], overhang[i]], axis=1)
            self.startfield = self.field

        # cartesien coordinates
        else:
            self.field = np.random.uniform(low=-self.pdfr, high=self.pdfr,
                                       size=(self.extendx * 2 + 1,
                                             self.extendy * 2 + self.res_z,
                                             self.extendz * 2 + self.res_y))
            self.startfield = self.field
        self.keep = keep

    def roll_one_timestep(self):
        # move random fields one position back
        self.field = np.roll(self.field, -1, axis=0)

        if self.polar:
            field_ex2 = np.random.uniform(low=-self.pdfr, high=self.pdfr,
                                                     size=(self.extendy * 2 + self.res_z,
                                                           self.res_y))
            overhang = np.random.uniform(low=-self.pdfr, high=self.pdfr,
                                       size=(self.extendz * 2 + self.res_z,
                                             self.extendy *2))
            self.field[-1, :, :] = np.concatenate([field_ex2, field_ex2,field_ex2, overhang], axis=1)
        else:
            self.field[-1, :, :] = np.random.uniform(low=-self.pdfr, high=self.pdfr,
                                                 size=(self.extendy * 2 + self.res_z,
                                                       self.extendz * 2 + self.res_y))
    def fill_from_startfield(self, index):
        self.field = np.roll(self.field, -1, axis=0)
        self.field[-1, :, :] = self.startfield[index, :, :]

    def next_field(self, period=0, timestep=0, periodic=True):
    # create new field depending on periodic
        if period == 0 or not periodic:
            self.roll_one_timestep()
        else:
            s_field_ind = -(period - timestep - self.extendx) + self.extendx + 1
            if s_field_ind < 0:
                self.roll_one_timestep()
            else:
                # for last steps: fill from start field
                self.fill_from_startfield(s_field_ind)


class CorrField(Field):

    def __init__(self, name, rand_field, corr_matrix, outside_flag, polar, keep=False):
        Field.__init__(self, name)
        ############################################ initialising field ################################################
        self.polar = polar
        self.field = np.zeros((0, rand_field.res_z, rand_field.res_y))
        self.outside = outside_flag
        self.corr_matrix = corr_matrix
        self.keep = keep
        self.rand_field = rand_field

    # manual correlation (not used, very slow)
    def add_field(self):
        new_field = np.zeros((1, self.rand_field.res_z, self.rand_field.res_y))
        for k in range(self.rand_field.res_y):
            for j in range(self.rand_field.res_z):
                if not self.outside[j, k]:
                    new_field[0, j, k] = np.sum(self.rand_field.field[:, j:j + 2 * int(self.rand_field.nfz[j, k]) + 1,
                                                k:k + 2 * int(self.rand_field.nfy[j, k]) + 1] *
                                                self.corr_matrix[k][j].matrix)
        if self.keep:
            self.field = np.concatenate((self.field, new_field), axis=0)
        else:
            self.field = new_field

    # manual correlation (not used, very slow)
    def add_field_conv_indiv(self):
        new_field = np.zeros((1, self.rand_field.res_z, self.rand_field.res_y))
        for k in range(self.rand_field.res_y):
            for j in range(self.rand_field.res_z):
                if not self.outside[j, k]:
                    new_field[0, j, k] = np.squeeze(convolve(self.rand_field.field[:, j:j + 2 * int(self.rand_field.nfz[j, k]) + 1,
                                                k:k + 2 * int(self.rand_field.nfy[j, k]) + 1],
                                                self.corr_matrix[k][j].matrix, mode='valid'))
        self.field = np.concatenate((self.field, new_field), axis=0)

    def add_field_conv(self):
    ############################################ filter process ########################################################
        if self.rand_field.polar:
        # for polar coordinates: filtering of 3 sectors and than extracting middle sector
            # convolve of 3 sector random field
            new_field_ex = convolve(self.rand_field.field, self.corr_matrix[0][0].matrix, mode='valid')

            # allocation of 1 sector field
            new_field = np.zeros((1,self.rand_field.res_z,self.rand_field.res_y))

            # write middle sector in 1 sector field
            new_field[0,:,:] = new_field_ex[0,:,self.rand_field.res_y:self.rand_field.res_y+self.rand_field.res_y]
        else:
        # cartesian coordinates
            new_field = convolve(self.rand_field.field, self.corr_matrix[0][0].matrix, mode='valid')

        # save all fields for keep = True
        if self.keep:
            self.field = np.concatenate((self.field, new_field), axis=0)
        else:
            self.field = new_field


class OutField(Field):
# object output field

    def __init__(self, name, rand_field, dt_split=1, keep=False):
        Field.__init__(self, name)
        # input parameter
        self.keep = keep
        self.res_z = rand_field.res_z
        self.res_y = rand_field.res_y
        self.field = np.zeros((0, rand_field.res_z, rand_field.res_y))
        self.dt_split = dt_split
        self.initialised = 0
        self.first_field = self.field
        self.running_sum = np.zeros((self.res_z, self.res_y))
        self.running_square_sum = np.zeros((self.res_z, self.res_y))
        self.steps = 0

    def add_field(self, field):
    # add fields for fluctuating variables to out field array
        # clip temperatre fluctuations
        if self.name == 'T' or self.name == 'Tt':
            field[field < 900] = 900
            field[field > 2450] = 2450
        if not self.keep and self.field.shape[0] > 0:
            self.field = self.field[-1, :, :].reshape((1, self.res_z, self.res_y))
        if self.dt_split > 1 and self.field.shape[0] > 0:
        #for timestep split, ad interpolated intermediate fields
            self.add_intermediate_fields(self.field[-1, :, :], field, self.dt_split)
        self.field = np.concatenate((self.field, field.reshape((1, self.res_z,
                                                                self.res_y))), axis=0)

        # statistics
        if self.initialised == 0:
            self.first_field = self.field
            self.running_sum = np.squeeze(self.field)
            self.running_square_sum = np.squeeze(self.field**2)
            self.steps = 1
            self.initialised = 1
        else:
            self.running_sum = self.running_sum + np.sum(self.field[-self.dt_split:, :, :], axis=0)
            self.running_square_sum = self.running_square_sum + np.sum(self.field[-self.dt_split:, :, :]**2, axis=0)
            self.steps = self.steps + self.dt_split

    def add_intermediate_fields(self, startfield, endfield, dt_split):
    # interpolation between timesteps
        endfield = endfield.reshape((1, self.res_z, self.res_y))
        startfield = startfield.reshape((1, self.res_z, self.res_y))
        for i in range(1, dt_split):
            self.field = np.concatenate((self.field, i / dt_split * endfield +
                                         (dt_split - i) / dt_split * startfield),
                                        axis=0)

    def add_last_fields(self):
        self.add_intermediate_fields(self.field[-1, :, :], self.first_field[0, :, :], self.dt_split)
        self.running_sum = self.running_sum + np.sum(self.field[-(self.dt_split - 1):, :, :], axis=0)
        self.running_square_sum = self.running_square_sum + np.sum(self.field[-(self.dt_split - 1):, :, :]**2, axis=0)
        self.steps = self.steps + self.dt_split - 1

    def get_average_until_timestep(self, index=0):
        if self.keep:
            return Field.get_average_until_timestep(self, index)
        else:
            return self.running_sum / self.steps

    def get_var_until_timestep(self, index=0):
        if self.keep:
            return Field.get_var_until_timestep(self, index)
        else:
            out = np.zeros(self.running_square_sum.shape)
            avg = self.get_average_until_timestep()
            flag = np.logical_not(avg == 0)
            out[flag] = (self.running_square_sum[flag] / self.steps) - avg[flag]**2
            return out


class RstTransform:
# create RST trasform matrix
    def __init__(self, dataarrays, res_phi_y, res_r_z, fids=('u', 'v', 'w')):
        self.fids = fids
        self.rst = np.zeros((res_r_z, res_phi_y, len(fids), len(fids)))
        self.rst_transform = np.zeros((res_r_z, res_phi_y, len(fids), len(fids)))
        ids = [el.standardname for el in dataarrays]
        rst_comps = [[dataarrays[ids.index(id_switch(i, j, ids))].data
                      for i in fids]
                     for j in fids]

        for i in range(res_r_z):
            for j in range(res_phi_y):
                self.rst[i, j, :, :] = np.array([[rst_comps[k][l][i, j]
                                                  for k in range(len(fids))]
                                                 for l in range(len(fids))])

        for i in range(res_r_z):
            for j in range(res_phi_y):
                self.rst_transform[i, j, :, :] = npc.cholesky(self.rst[i, j, :, :])

    def adapt_2_profile(self, corrfields, outfields, profile, outfields_2):
    # adapt filtered field to RST and mean value
        # identify names of fluctuating fields
        corr_names = [el.name for el in corrfields]
        out_names = [el.name for el in outfields]

        corr_inds = [corr_names.index(ident) for ident in self.fids]
        out_inds = [out_names.index(ident) for ident in self.fids]

        # multiply RST transform matrix with correlated fields
        for i, eli in enumerate(self.fids):
            summe = 0
            for j, elj in enumerate(self.fids):
                summe = summe + self.rst_transform[:, :, i, j] * corrfields[corr_inds[j]].get_last_field()
            rst_flag = True

            # superpose fluctuating field with mean fields
            summe = summe + profile.get_data(eli, rst_flag)
            outfields[out_inds[i]].add_field(summe)

        # for polar coordinates, calculate angles and total pressure
        if profile.polar:
            self.calc_hydra_input(outfields, profile)

            outfields_2[0].add_field(self.whirl)
            outfields_2[1].add_field(self.pitch)
            outfields_2[2].add_field(self.ptotal)


    def calc_hydra_input(self, outfields, profile):
    # calculate flow angles and total pressure fluctuations for HYDRA export
        u = outfields[0].field[-1]
        v = outfields[1].field[-1]
        w = outfields[2].field[-1]
        vel_mag = np.sqrt(u**2 + v**2 + w**2)
        self.whirl = -np.arctan(v / u)
        self.pitch = np.arctan(w / u)
        p = 820758
        rho = 1.71066
        self.ptotal = p + rho / 2 * vel_mag ** 2


def id_switch(id1, id2, ids):
    if id1 + id2 in ids:
        return id1 + id2
    else:
        return id2 + id1

class Prepare_POD:

    def __init__(self, hdf5, timesteps, res_y, res_z, fluct_comps, fluct_comps_2, dt, coord, coord_comps):
        self.hdf5 = hdf5
        self.timesteps = timesteps
        self.res_y = res_y
        self.res_z = res_z
        self.fluct_comps = fluct_comps
        self.fluct_comps_2 = fluct_comps_2
        self.num_points = self.res_y * self.res_z
        self.num_snapshots = self.timesteps
        self.num_components = self.fluct_comps.__len__()
        self.num_components_2 = self.fluct_comps_2.__len__()
        self.dt = dt
        self.coordinates = coord
        self.coord_comps = coord_comps
        self.num_coord_comps = self.coord_comps.__len__()

    def prepare_for_POD(self, out_fields, out_fields_2, i):
        if self.hdf5:

            velo = [el.get_indexed_field(-1) for el in out_fields_2]
            velo = [el.reshape(el.size) for el in velo]
            velo[0] = velo[0] * 180 / np.pi
            velo[1] = velo[1] * 180 / np.pi

            if 'Tt' in self.fluct_comps:
                if i == 0:
                    self.fluct_comps_pod = self.fluct_comps_2[:]
                    self.fluct_comps_pod.append(self.fluct_comps[-1])
                    self.num_components_pod = self.fluct_comps_pod.__len__()
                    self.A = np.transpose(np.zeros((self.timesteps, self.res_y * self.res_z * self.num_components_pod)))
                velo1 = [el.get_indexed_field(-1) for el in out_fields]
                velo1 = [el.reshape(el.size) for el in velo1]
                velo.append(velo1[-1])
                self.A[:, i] = np.concatenate([velo[self.fluct_comps_pod.index(name)][:] for name in self.fluct_comps_pod])

            else:
                if i == 0:
                    self.A = np.transpose(np.zeros((self.timesteps, self.res_y * self.res_z * self.num_components_2)))
                self.fluct_comps_pod = self.fluct_comps_2[:]
                self.num_components = self.fluct_comps_pod.__len__()
                self.A[:, i] = np.concatenate([velo[self.fluct_comps_2.index(name)][:] for name in self.fluct_comps_2])

        else:
            if i == 0:
                self.A = np.transpose(np.zeros((self.timesteps, self.res_y * self.res_z * self.num_components)))
            velo = [el.get_indexed_field(-1) for el in out_fields]
            velo = [el.reshape(el.size) for el in velo]
            self.A[:,i] = np.concatenate([velo[self.fluct_comps.index(name)][:] for name in self.fluct_comps])
            self.fluct_comps_pod = self.fluct_comps[:]

        if i == self.timesteps -1:
            self.mean_field = np.mean(self.A, axis = 1)
            for j in range(0, self.timesteps):
                self.A[:, j] = self.A[:, j] - self.mean_field[:]

