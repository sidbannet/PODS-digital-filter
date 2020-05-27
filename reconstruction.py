import numpy as np
from numpy import linalg
import math
import nplotlib as nplt
import h5py
from scipy import interpolate
import prepare_traverse as preptrav

class PODFS_reader:
    def __init__(self, filename):
        self.filename = filename

        self.read_hdf5()

    def read_hdf5(self):

        print('########### read HDF5 file ##############')
        print(self.filename)
        f = h5py.File(self.filename,'r')

        main = f.require_group('main')

        # read attributes
        self.num_modes = main.attrs.__getitem__('number_POD_modes')
        print('            number of POD modes: ', self.num_modes)
        self.period = main.attrs.__getitem__('period')
        print('            period: ', self.period)
        self.coordinate_comp = main.attrs.__getitem__('coordinates')
        print('            coordinates', self.coordinate_comp)
        self.num_coord_comp = main.attrs.__getitem__('number_coord_components')
        self.num_points = main.attrs.__getitem__('number_points')
        print('            number of points', self.num_points)
        self.num_variables = main.attrs.__getitem__('number_variables')
        self.variables = main.attrs.__getitem__('variables')
        print('            variables: ', self.variables)

        # read fourier coefficients
        self.fourier_coeff = main.__getitem__('FC')[()]

        # read fourier coefficients ranking
        self.FC_ind = main.__getitem__('Rank_FC')[()]

        # read coordinates
        self.coordinates = main.__getitem__('coordinates')[()]

        # read mean field
        self.mean_field = main.__getitem__('mean')[()]

        # read number of fourier coefficients per mode
        self.num_FC = main.__getitem__('number_FC_per_mode')[()]
        print('            number of Fourier coefficients per mode:', self.num_FC)

        # read spartial modes
        self.mode = np.array(np.zeros((self.num_points * self.num_variables, self.num_modes)))

        spatial_modes = main.require_group('modes')
        for i in range (0,self.num_modes):
            counter = '%4.4i'% (i+1)
            self.mode[:,i] = spatial_modes.__getitem__('mode_'+counter)[()]


class Reconstruct:
    def __init__(self, podfs, outputfolder, plot_flag, time):
        self.podfs = podfs
        self.outputfolder = outputfolder
        self.plot_flag = plot_flag

        # set time vector for evaluation of velocity fields
        self.time = time
        # self.time = np.linspace(0, 0.5*self.podfs.period, 200)

        # get velocity fields from modes and fourier coefficients
        self.get_velocity()

        # plot velocity fields
        self.plot_fields()


    def get_velocity(self):

        # allocate arrays
        a_temp = np.array(np.zeros((self.podfs.num_modes, self.time.__len__())), dtype=complex)
        a = np.array(np.zeros((self.podfs.num_modes, self.time.__len__())))
        A_temp = np.array(np.zeros((self.podfs.num_points * self.podfs.num_variables, self.podfs.num_modes)))
        A_temp1 = np.array(np.zeros((self.podfs.num_points * self.podfs.num_variables)))
        A = np.array(np.zeros((self.podfs.num_points * self.podfs.num_variables, self.time.__len__())))

        self.BC_fields = [BC_Field(name, self.podfs.num_points, self.time.__len__()) for name in self.podfs.variables]
        self.mean_fields = [mean_Field(name, self.podfs.num_points) for name in self.podfs.variables]
        self.mode_fields = [mode_Field(name, self.podfs.num_points, self.podfs.num_modes) for name in self.podfs.variables]

        # loop over time values
        for i in range(0, self.time.__len__()):

            # loop over pod modes
            for x in range(0, self.podfs.num_modes):

                # loop over fourier coefficients
                for n in range(0, self.podfs.num_FC[x]):
                    k = self.podfs.FC_ind[x, n] - self.podfs.FC_ind.shape[1] / 2
                    a_temp[x,i] += self.podfs.fourier_coeff[self.podfs.FC_ind[x,n],x] * np.exp(2 * np.pi * 1j * self.time[i] * k / self.podfs.period)
                    # temporal modes
                    a[x,i] = a_temp[x,i].real

                # spartial mode * temporal mode
                A_temp[:,x] = self.podfs.mode[:,x] *  a[x,i]

                # write spartial modes to fieldarray
                if i == 0:
                    for h in range(self.podfs.num_variables):
                        self.mode_fields[h].field[:,x] = self.podfs.mode[self.podfs.num_points * h:self.podfs.num_points * (h+1), x]

            # sum of all modes
            A_temp1[:] = np.sum(A_temp, axis=1)

            # superposition of fluctuation and mean field
            A[:,i] = A_temp1[:] + self.podfs.mean_field

            # build velocity components
            for j in range(self.podfs.num_variables):
                self.BC_fields[j].field[:,i] = A[self.podfs.num_points * j:self.podfs.num_points * (j+1), i]
                if i == self.time.__len__()-1:
                    self.BC_fields[j].min = np.nanmin(self.BC_fields[j].field)
                    self.BC_fields[j].max = np.nanmax(self.BC_fields[j].field)

        print('finish reconstruction of velocity fields')

        # plot reconstructed temporal modes
        if plot_flag[4]:
            for kk in range(0, self.podfs.num_modes):
                nplt.timeseries('Temp_mode_recon' + str(kk), a[kk, :], self.time, '\,',
                                self.outputfolder + 'POD_tmode_reconstructed' + str(kk))

        # calculate time average fields
        for l in range(self.podfs.num_variables):
            self.mean_fields[l].field[:] = np.mean(self.BC_fields[l].field, axis=1)



    def plot_fields(self):
        # get coordinates

        if 'r' in self.podfs.coordinate_comp:
            radius = self.podfs.coordinates[:,1]
            phi_rad = self.podfs.coordinates[:,0]*np.pi/180
        else:
            self.x = self.podfs.coordinates[:,0]
            self.y = self.podfs.coordinates[:,1]
            self.z = self.podfs.coordinates[:,2]

            # calculate polar coordinates
            radius = np.sqrt(self.y ** 2 + self.z ** 2)
            phi_rad = np.arctan(-self.y / self.z)

        # define new grid
        ri = np.linspace(np.min(radius), np.max(radius), 90)
        phii = np.linspace(np.min(phi_rad), np.max(phi_rad), 110)

        # calculate new grid
        phi, r = np.meshgrid(phii, ri)

        # calculate cartesian coordinates of new grid
        y = -r * np.sin(phi)
        z = r * np.cos(phi)

        # for 3 sectors expand coordinate arrays
        if plot_flag[2]:
            r1 = np.concatenate((r[:,:108], r[:,:108], r[:,:108]), axis=1)
            phi1 = np.concatenate((phi[:,:108]-22.5*np.pi/180, phi[:,:108], phi[:,:108]+22.5*np.pi/180), axis=1)
            y1 = -r1 * np.sin(phi1)
            z1 = r1 * np.cos(phi1)

        ###################################### plot mean fields ########################################################
        if self.plot_flag[0]:

            # interpolate mean fields to new grid
            self.dataarrays_mean = [preptrav.DataColumn(self.mean_fields[el].name, interpolate.griddata(
                (radius, phi_rad),self.mean_fields[el].field, (r, phi), method='linear'))
                               for el in range(0,self.podfs.num_variables)]

            # for 3 sectors
            if self.plot_flag[2]:

                for f in range(self.podfs.num_variables):
                    # plot mean fields
                    nplt.contourf2(self.dataarrays_mean[f].headerfield+'_mean', self.dataarrays_mean[f].unit, y1, z1, np.concatenate((self.dataarrays_mean[f].data[:,:108], self.dataarrays_mean[f].data[:,:108], self.dataarrays_mean[f].data[:,:108]), axis=1), 100, 'y', 'z', self.dataarrays_mean[f].headerfield, self.outputfolder + self.dataarrays_mean[f].headerfield+ '_mean',
                                   figsize=(16, 16 * ((90-3) / (110-3)/3)))
                    nplt.close(self.dataarrays_mean[f].headerfield+'_mean')

            # for 1 sector
            else:
                for f in range(self.podfs.num_variables):
                    # plot mean fields
                    nplt.contourf2(self.dataarrays_mean[f].headerfield+'_mean', self.dataarrays_mean[f].unit, y, z, self.dataarrays_mean[f].data, 100, 'y', 'z', self.dataarrays_mean[f].headerfield, self.outputfolder + self.dataarrays_mean[f].headerfield+ '_mean',
                                   figsize=(8, 8 * 90 / 110))
                    nplt.close(self.dataarrays_mean[f].headerfield+'_mean')

        ######################################## plot instantan velocity fields ########################################
        if self.plot_flag[1]:
            # loop over output timesteps
            for i in range(0, self.time.__len__()):
                step = str("{:04d}".format(i + 1))
                # interpolate mean fields to new grid
                self.dataarrays = [preptrav.DataColumn(self.BC_fields[el].name, interpolate.griddata(
                    (radius, phi_rad), self.BC_fields[el].field[:,i], (r, phi), method='linear'))
                                        for el in range(0, self.podfs.num_variables)]

                # for 3 sectors
                if self.plot_flag[2]:

                    for f in range(self.podfs.num_variables):
                        # plot mean fields
                        nplt.contourf(self.BC_fields[f].min, self.BC_fields[f].max, self.dataarrays[f].headerfield + '_reconstructed_' + step, y1,
                                       z1, np.concatenate((self.dataarrays[f].data[:, :108],
                                                           self.dataarrays[f].data[:, :108],
                                                           self.dataarrays[f].data[:, :108]), axis=1), 100, 'y',
                                       'z', self.dataarrays[f].headerfield + ' in ' + self.dataarrays[f].unit,
                                       self.outputfolder + self.dataarrays[f].headerfield + '_reconstructed_' + step,
                                       figsize=(16, 16 * ((90 - 3) / (110 - 3) / 3)))
                        nplt.close(self.dataarrays_mean[f].headerfield + '_reconstructed_' + step)

                # for 1 sector
                else:
                    for f in range(self.podfs.num_variables):
                        # plot instant fields
                        nplt.contourf(self.BC_fields[f].min, self.BC_fields[f].max, self.dataarrays[f].headerfield + '_reconstructed_' + step, y,
                                       z, self.dataarrays[f].data, 100, 'y', 'z',
                                       self.dataarrays[f].headerfield + ' in ' + self.dataarrays[f].unit,
                                       self.outputfolder + self.dataarrays[f].headerfield + '_reconstructed_' + step,
                                       figsize=(8, 8 * 90 / 110))
                        nplt.close(self.dataarrays_mean[f].headerfield + '_reconstructed_' + step)


        ########################################## plot of spatial modes ###############################################
        if plot_flag[3]:
            for k in range(0, self.podfs.num_modes):
                # interpolate mean fields to new grid
                self.dataarrays_mode = [preptrav.DataColumn(self.mode_fields[el].name, interpolate.griddata(
                    (radius, phi_rad), self.mode_fields[el].field[:, k], (r, phi), method='linear'))
                                   for el in range(0, self.podfs.num_variables)]

                step1 = str("{:04d}".format(k))

                if plot_flag[2]:

                    for f in range(0, self.podfs.num_variables):
                        nplt.contourf2('mode_' + self.dataarrays_mode[f].headerfield + step1, '-', y1, z1, np.concatenate((self.dataarrays_mode[f].data[:, :108],
                                                               self.dataarrays_mode[f].data[:, :108],
                                                               self.dataarrays_mode[f].data[:, :108]), axis=1), 100, 'y', 'z', 'mode_' + self.dataarrays_mode[f].headerfield + step1,
                                          self.outputfolder + '/mode_' + self.dataarrays_mode[f].headerfield + step1,
                                          figsize=(16, 16 * ((90 - 3) / (110 - 3) / 3)))
                        nplt.close('mode_' + self.dataarrays_mode[f].headerfield + step1)

                else:
                    for f in range(0, self.podfs.num_variables):
                        # plot mode fields
                        nplt.contourf2('mode_' + self.dataarrays_mode[f].headerfield + step1, '-', y, z,
                                       self.dataarrays_mode[f].data, 100, 'y', 'z', 'mode_' + self.dataarrays_mode[f].headerfield + step1,
                                       self.outputfolder + '/mode_' + self.dataarrays_mode[f].headerfield + step1,
                                       figsize=(8, 8 * 90 / 110))
                        nplt.close('mode_' + self.dataarrays_mode[f].headerfield + step1)


class Field:
    def __init__(self, name):
        self.name = name

class BC_Field(Field):
    def __init__(self, name, num_points, num_timesteps):
        Field.__init__(self, name)

        self.field = np.array(np.zeros((num_points, num_timesteps)))

class mean_Field(Field):
    def __init__(self, name, num_points):
        Field.__init__(self, name)

        self.field = np.array(np.zeros((num_points,)))

class mode_Field(Field):
    def __init__(self, name, num_points, num_modes):
        Field.__init__(self, name)

        self.field = np.array(np.zeros((num_points, num_modes)))

# podfs reader
podfs = PODFS_reader('../Output/test_polar/HYDRA/PODFS/PODFS_HYDRA.hdf5')

# plot input
plot_mean = True
plot_inst = False
plot_3_sec = True
plot_modes = False
plot_temp_modes = True
plot_flag = [plot_mean, plot_inst, plot_3_sec, plot_modes, plot_temp_modes]
timestep = 2.0425e-4
num_timesteps = 200
time = np.linspace(0, timestep * num_timesteps, num_timesteps)

# reconstruction of flow fields
reconstruct = Reconstruct(podfs, '../Output/test_polar/HYDRA/PODFS/', plot_flag, time)