import numpy as np
import numpy_custom as npc
import prepare_traverse as preptrav
from scipy import interpolate
import matplotlib.pyplot as plt
import math
import os
import time
import nplotlib as nplt


class Filter:

    def __init__(self, filename='none', res=0.1, outdir='cwd', dt_split=1, profile='hyperbolic-tangent',
                 fwidth=2, verbose=True, pdfr=np.sqrt(3), keep_all=False,
                 res_y=10, res_z=10, oy=0, oz=0, turb_prof='top-hat',
                 tu_intensity=0, lengthscale=3, y_range=(0, 0), z_range=(0, 0),
                 inner_diameter_norm=0, prof2d=False, mdot=0, bulk_velocity=0,
                 density=1, plot_flag=False):

        # setting parameters
        # time step split
        self.dt_split = dt_split
        self.plot_flag = plot_flag

        # filter width
        self.fwidth = fwidth

        # save profiles
        self.verbose = verbose
        if outdir is not 'cwd':
            self.outdir = outdir
        else:
            self.outdir = os.getcwd()

        # preparation of input field
        if filename.endswith('.prf') or prof2d:
            self.profile = Profil2D(filename, res, mdot, bulk_velocity, density)
        else:
            self.profile = Profil1DProf(res, profil=profile, turb_prof=turb_prof, u_bulk=bulk_velocity,
                                        tu_intensity=tu_intensity, filename=filename,
                                        res_y=res_y, res_z=res_z, oy=oy, oz=oz,
                                        lengthscale=lengthscale, y_range=y_range, z_range=z_range,
                                        inner_diameter_norm=inner_diameter_norm)
        self.profile.prepare_for_filter()
        if self.plot_flag[0] and self.plot_flag[1]:
            self.profile.plot_input_profile(outdir=self.outdir)

        # setting size of filter matrix
        self.nfx = self.profile.lnx * self.fwidth
        self.nfy = self.profile.lny * self.fwidth
        self.nfz = self.profile.lnz * self.fwidth

        # initialising field of filter matrix
        # format: self.filter_coeffs[y_index][z_index]
        self.filter_coeffs = [[[] for i in range(self.nfy.shape[0])] for j in range(self.nfy.shape[1])]
        self.fill_coeff_field()

        # initialising of fields for: 1. Random field, 2. Correlated field, 3. Output field
        self.rand_fields = [RandomField(name, self.nfx, self.nfy, self.nfz, pdfr) for name in ['u', 'v', 'w']]
        self.corr_fields = [CorrField(el.name, el, self.filter_coeffs, self.profile.outside_flag, keep_all)
                            for el in self.rand_fields]
        self.out_fields = [OutField(el.name, el, self.dt_split, keep_all) for el in self.rand_fields]

        # resolution
        self.res = res

        # getting filter timestep from resolution and velocity
        self.dt_filter = self.res / \
                         np.mean(self.profile.dataarrays[self.profile.get_datanames().index('u')].data[
                                     np.logical_not(self.profile.outside_flag)])
        self.mean_diff = np.zeros((0, 4))
        self.max_diff = np.zeros((0, 4))

        if self.verbose:
            if not os.path.exists(self.outdir + '/PODFS'):
                os.makedirs(self.outdir + '/PODFS')
        self.u_bulk = bulk_velocity

    def fill_coeff_field(self):
        # At the moment same coefficient matrix for all grid points ==> maybe make different for lenght scale distribution
        matrix = FilterMatrix(self.nfx[0, 0], self.nfy[0, 0], self.nfz[0, 0],
                              self.profile.lnx[0, 0], self.profile.lny[0, 0], self.profile.lnz[0, 0])
        for i in range(len(self.filter_coeffs)):
            for j in range(len(self.filter_coeffs[i])):
                self.filter_coeffs[i][j] = matrix

    def filtering(self, timesteps, periodic=True, conv_stat=False):
        if conv_stat:
            fig = plt.figure()

        if self.dt_split > 1:
            timesteps = math.ceil(timesteps / self.dt_split)

        starttime = time.time()
        for i in range(0, timesteps):
            print('Processing Inlet {:d}/{:d}'.format(i + 1, timesteps))

            # creation of spatial correlations
            for el in self.corr_fields:
                el.add_field()

            # adapt profile to Reynolds stresses
            self.profile.rst_transform.adapt_2_profile(self.corr_fields,
                                                       self.out_fields,
                                                       self.profile)
            # save profiles
            if self.verbose:
                self.save_field(i)

            # prepare random fields for next time step
            for el in self.rand_fields:
                el.next_field(timesteps, i, periodic)

            if conv_stat:
                self.diff_stats()
                self.plot_diff_stats(fig)
                if self.verbose and (i + 1) % 10 == 0:
                    fig2 = self.plot_k_field(i + 1).savefig(self.outdir + '/PODFS/k_field_{:04d}.png'.format(i + 1),
                                                            dpi=200)
                    plt.close(fig2)

            # approximation of remaining time
            endtime = time.time()
            avg_time = (endtime - starttime) / (i + 1)
            remain_time = (timesteps - (i + 1)) * avg_time
            print('Verbleibende Zeit: {:02d}h {:02d}m {:02d}s'.format(int(remain_time // 3600),
                                                                      int((remain_time % 3600) // 60),
                                                                      int(remain_time % 60)))

            # plot random fields
            if self.plot_flag[0] and self.plot_flag[2]:
                step = str(i + 1)
                nplt.contourf2('rand' + step, self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                               self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                               self.rand_fields[0].field[0][0:self.profile.res_z, 0:self.profile.res_y], 100, 'y', 'z',
                               'r', self.outdir + '/rand' + step, figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                nplt.close('rand' + step)

            # plot correlated fields
            if self.plot_flag[0] and self.plot_flag[3]:
                step = str(i + 1)
                nplt.contourf(np.min(self.corr_fields[0].field[0]), np.max(self.corr_fields[0].field[0]),
                              'filtered' + step, self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                              self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                              self.corr_fields[0].field[i], 100, 'y', 'z',
                              'u', self.outdir + '/filtered' + step,
                              figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                nplt.close('filtered' + step)

            # plot adapted fields
            if self.plot_flag[0] and self.plot_flag[4]:
                step = str(i + 1)
                nplt.contourf(np.min(self.out_fields[0].field[0]), np.max(self.out_fields[0].field[0]),
                              'adapted' + step, self.profile.dataarrays[self.profile.get_datanames().index('y')].data,
                              self.profile.dataarrays[self.profile.get_datanames().index('z')].data,
                              self.out_fields[0].field[i], 100, 'y', 'z',
                              'u', self.outdir + '/adapted' + step, figsize=(8, 8 * self.profile.res_z / self.profile.res_y))
                nplt.close('adapted' + step)

        if periodic and self.dt_split > 1:
            for el in self.out_fields:
                el.add_last_fields()
            if self.verbose:
                self.save_field(timesteps, last=True)

        # save settings file
        if self.verbose:
            with open(self.outdir + '/PODFS/00_params.txt', 'w') as out_file:
                out_file.write('Datei erstellt: ' + time.ctime() + '\n')
                out_file.write('Zeitschritt Filter: {:14e}\n'.format(self.dt_filter))
                out_file.write('Interpolation Filterwerte: {:d}\n'.format(self.dt_split))
                out_file.write('Zeitschritt LES: {:14e}\n'.format(self.dt_filter / self.dt_split))
                out_file.write('Inputprofil: ' + self.profile.name + '\n')
                out_file.write('Aufloesung: {:.8f}\n'.format(self.res))
                out_file.write('Filterweite: {:d}\n'.format(self.fwidth))
                out_file.write('Anzahl Zeitschritte Filter: {:d}\n'.format(timesteps))
                out_file.write('Anzahl Zeitschritte Dateien: {:d}\n'.format(timesteps * self.dt_split))
                out_file.write('Periodisch: {:b}\n'.format(periodic))
                out_file.write('Anzahl Zeilen: {:d}\n'.format(self.out_fields[0].get_last_field().size))

    def save_field(self, timestep, last=False):
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
            velo = [el.get_indexed_field(n) for el in self.out_fields]
            velo = [el.reshape(el.size) for el in velo]

            with open(self.outdir + '/PODFS/inlet.dat' + str(timestep * self.dt_split + n + 2), 'w') as out_file:
                for m in range(0, coord[0].size):
                    out_file.write('{:.12f} {:.12f} {:.12f} {:.12f} {:.12f} {:.12f}\n'.format(
                        coord[0][m], coord[1][m], coord[2][m], velo[0][m], velo[1][m], velo[2][m]
                    ))

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
        crange = (np.percentile(k_prof, 10)*0.9, np.percentile(k_prof, 90)*1.1)
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
        return 0.5 * sum([np.mean(el.get_time_fluctuation_until_timestep(timestep) ** 2, axis=0)
                          for el in self.out_fields if el.name in ['u', 'v', 'w']])

    def get_k_diff_profile(self, timestep=0):
        return self.get_k_from_outfields(timestep) - self.profile.get_data('k')

    def get_k_mean_diff(self, timestep=0):
        return np.mean(np.abs(self.get_k_diff_profile(timestep)))

    def get_k_max_diff(self, timestep=0):
        return np.max(np.abs(self.get_k_diff_profile(timestep)))

    def diff_stats(self, timestep=0):
        self.mean_diff = np.append(self.mean_diff, [np.append(self.get_mean_diff(), self.get_k_mean_diff(timestep))],
                                   axis=0)
        self.max_diff = np.append(self.max_diff, [np.append(self.get_max_diff(), self.get_k_max_diff(timestep))],
                                  axis=0)

    def plot_diff_stats(self, fig):
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
        plt.draw()
        plt.pause(1e-17)
        time.sleep(0.00001)


class Profil:

    def __init__(self, name):
        self.dataarrays = []
        self.res_y = 0
        self.res_z = 0
        self.lengthscale_x = 0
        self.lengthscale_y = 0
        self.lengthscale_z = 0
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

    def get_data(self, name):
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
        self.lny = self.lnx
        self.lnz = self.lnx

    def plot_input_profile(self, outdir):
        for el in self.dataarrays:
            if el.standardname not in ['x', 'y', 'z']:
                plot_id = el.standardname
                nplt.contourf2(plot_id,
                               self.get_data('y'),
                               self.get_data('z'),
                               el.data, 100, 'y', 'z',
                               plot_id, outdir + '/' + plot_id, figsize=(8, 8 * self.res_z / self.res_y))
                nplt.close(plot_id)


class Profil1DProf(Profil):

    def __init__(self, res, profil='hyperbolic-tangent', turb_prof='top-hat', u_bulk=0, tu_intensity=0,
                 filename='none', res_y=10, res_z=10, oy=0, oz=0, lengthscale=1,
                 y_range=(0, 0), z_range=(0, 0), inner_diameter_norm=0):
        if filename is not 'none':
            name = filename + ' ' + profil
        else:
            name = profil
        Profil.__init__(self, name)
        self.profil = profil
        self.res = res
        self.filename = filename
        if not y_range == (0, 0):
            self.y_range = y_range
        else:
            self.y_range = (oy - res_y / 2 * res, oy + res_y / 2 * res)

        if not z_range == (0, 0):
            self.z_range = z_range
        else:
            self.z_range = (oz - res_z / 2 * res, oz + res_z / 2 * res)

        self.lengthscale_x = lengthscale * self.res
        self.turb_prof = turb_prof
        self.u_bulk = u_bulk
        self.tu_intensity = tu_intensity
        self.oy = np.mean(self.y_range)
        self.oz = np.mean(self.z_range)
        self.fr = None
        self.inner_diameter_norm = inner_diameter_norm
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
        if self.filename == 'none':
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
            fr = preptrav.FileReader(self.filename)
            fr.open_file()
            fr.extract_column_data()

            self.fr = fr

            y_data = fr.get_data('y')
            y_data = (y_data - min(y_data)) / (max(y_data) - min(y_data))
            y_data = np.abs(y_data / 2 - 0.5)
            fieldnames = ['u', 'v', 'w', 'uu', 'vv', 'ww', 'uv', 'uw', 'vw']
            for el in fr.data_columns:
                if el.standardname in fieldnames:
                    val_data = el.data
                    if val_data[-1] == 0:
                        val_data[-1] = -1e-5
                    ref_val = val_data[-1]

                    val_data = val_data / ref_val
                    prof_functions.append(ProfFunc(el.standardname,
                                                   interp=interpolate.interp1d(y_data, val_data,
                                                                               fill_value='extrapolate'),
                                                   ref_val=ref_val))
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

        if self.profil == 'hyperbolic-tangent':
            z_data = 0
            outside = y_data * 0
        elif self.profil == 'double-hyperbolic-tangent':
            outside = y_data * 0
        elif self.profil == 'circular-hyperbolic-tangent':
            y_data = np.sqrt(y_data ** 2 + z_data ** 2)
            outside = y_data > 0.5
            z_data = 0
        elif self.profil == 'ring-hyperbolic-tangent':
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
                                                       func.get_func(self.profil)(self.u_bulk, y_data, z_data)))
            self.dataarrays[-1].data[outside] = 0
        self.dataarrays.append(preptrav.DataColumn('k', 0.5 * (self.get_data('uu') +
                                                               self.get_data('vv') +
                                                               self.get_data('ww'))))


class ProfFunc:

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

    def __init__(self, filename, res, mdot, bulk_velocity, density):
        Profil.__init__(self, filename)
        self.filename = filename
        self.res = res
        self.mdot = mdot
        self.bulk_velocity = bulk_velocity
        self.den = density
        self.fr = preptrav.FileReader(self.filename)
        self.fr.open_file()
        self.fr.extract_column_data()
        self.valid = False
        self.mode = 'none'
        self.rad_flag = 0
        self.gradients = []
        self.plane = None

    def prepare_for_filter(self):
        # prepare data for filter
        self.check_validity()

        if self.mode == '2-equation' and 'e' not in self.fr.get_standardnames():
            self.convert_sdr_to_eps()
        elif self.mode == 'reynolds-stresses':
            self.check_reynolds_fields()

        self.interpolate_data()

        if self.mdot != 0.0 or self.bulk_velocity != 0.0:
            self.scale_vel()

        self.set_outside_flag()

        if self.mode == '2-equation':
            self.calculate_rst_from_k_eps()

        self.calculate_lengthscale()

        self.rst_transform = RstTransform(self.dataarrays, self.res_y, self.res_z)

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

        if rs_check1:
            self.mode = 'reynolds-stresses'
        elif turb_check:
            self.mode = '2-equation'

    def interpolate_data(self):
        # caching of coordinates
        x_data = self.fr.data_columns[self.fr.get_standardnames().index('x')].data
        y_data = self.fr.data_columns[self.fr.get_standardnames().index('y')].data
        z_data = self.fr.data_columns[self.fr.get_standardnames().index('z')].data

        # create plane object --> includes transform on y-z Plane
        self.plane = Plane(x_data, y_data, z_data)
        x_data = self.plane.x_data_rot
        y_data = self.plane.y_data_rot
        z_data = self.plane.z_data_rot

        # calculation of x-coordinate
        xval = np.mean(x_data)

        # determination of inner radius of inlet sector, setting velocity to zero inside this radius
        inner_radius = min(np.sqrt(y_data ** 2 + z_data ** 2))

        # determination of value range of y and z, creation of cartesian coordinate vectors
        yi = np.arange(np.min(y_data), np.max(y_data), self.res)
        zi = np.arange(np.min(z_data), np.max(z_data), self.res)

        # creation of grid
        y, z = np.meshgrid(yi, zi)
        # matrix form:
        #     | y_min ... y_max |        | z_min ... z_min |
        # y = |                 | ,  z = |                 |
        #     | y_min ... y_max |        | z_max ... z_max |
        # y(i, j) : 1. Index row --> varying z, 2. Index collumn --> varying y

        # size of matrix
        self.res_y = len(yi)
        self.res_z = len(zi)

        # interpolation of imported fields on new grid
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

        # find TI
        TI = np.sqrt(2. / 3. * self.dataarrays[self.get_datanames().index('k')].data) / \
             np.sqrt(self.dataarrays[self.get_datanames().index('u')].data ** 2 + \
                     self.dataarrays[self.get_datanames().index('v')].data ** 2 + \
                     self.dataarrays[self.get_datanames().index('w')].data ** 2)

        # find length scales
        L = self.dataarrays[self.get_datanames().index('k')].data ** (3 / 2) / \
            self.dataarrays[self.get_datanames().index('e')].data

        # find new velocities
        if self.mdot != 1 and self.bulk_velocity != 1:
            print('Scaling of massflow and bulk velocity at the same time not possible!')
            return
        elif self.mdot != 1:
            scale = self.mdot / mdot_old
            print('Scaling to Massflow: ', self.mdot, 'kg/s')
        elif self.bulk_velocity != 1:
            scale = self.bulk_velocity / udotn
            print('Scale to bulk velocity: ', self.bulk_velocity, 'm/s')

        U = self.dataarrays[self.get_datanames().index('u')].data * scale
        V = self.dataarrays[self.get_datanames().index('v')].data * scale
        W = self.dataarrays[self.get_datanames().index('w')].data * scale

        # find new k
        k = TI ** 2 * (U ** 2 + W ** 2 + V ** 2)

        # find new epsilon
        e = k ** (3 / 2) / L

        self.dataarrays[self.get_datanames().index('u')].data = U
        self.dataarrays[self.get_datanames().index('v')].data = V
        self.dataarrays[self.get_datanames().index('w')].data = W
        self.dataarrays[self.get_datanames().index('k')].data = np.nan_to_num(k)
        self.dataarrays[self.get_datanames().index('e')].data = np.nan_to_num(e)

    def calculate_rst_from_k_eps(self):
        # perhaps smothing of interpolated velocity fields to avoid high gradients at the boundaries
        # calculation of velocity gradients
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

        # smoothing of gradients ==> problems with ring sector
        #         for i in range(1, self.res_z - 1):
        #             for j in range(1, self.res_y - 1):
        #                 dudy[i, j] = np.mean(dudy[i - 1:i + 1, j - 1:j + 1])
        #                 dudz[i, j] = np.mean(dudz[i - 1:i + 1, j - 1:j + 1])
        #                 dvdy[i, j] = np.mean(dvdy[i - 1:i + 1, j - 1:j + 1])
        #                 dvdz[i, j] = np.mean(dvdz[i - 1:i + 1, j - 1:j + 1])
        #                 dwdy[i, j] = np.mean(dwdy[i - 1:i + 1, j - 1:j + 1])
        #                 dwdz[i, j] = np.mean(dwdz[i - 1:i + 1, j - 1:j + 1])

        # approximation dudx from incompressibility
        dudx = - dvdy - dwdz

        # setting all other gradients in flow direction to 0
        dvdx = np.zeros((self.res_z, self.res_y), dtype=np.float64)
        dwdx = np.zeros((self.res_z, self.res_y), dtype=np.float64)

        # caching gradients in object
        self.gradients = [[dudx, dudy, dudz],
                          [dvdx, dvdy, dvdz],
                          [dwdx, dwdy, dwdz]]

        # eddy viscosity model
        k = self.dataarrays[self.get_datanames().index('k')].data
        eps = self.dataarrays[self.get_datanames().index('e')].data
        flag = np.where(eps > 0)
        nu_t = np.zeros((self.res_z, self.res_y), dtype=np.float64)
        nu_t[flag] = 0.09 * k[flag] ** 2 / eps[flag]
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

        # write gradient fields into dataarrays
        name_rst = [['dudx', 'dudy', 'dudz'], ['dvdx', 'dvdy', 'dvdz'], ['dwdx', 'dwdy', 'dwdz']]
        for i in range(3):
            for j in range(3):
                self.dataarrays.append(preptrav.DataColumn(name_rst[i][j], self.gradients[i][j]))

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

    def calculate_lengthscale(self):
        # determining length scale for filter
        # at the moment: constant length scale for whole plane and all components, isotropic
        # output: matrix, length scale could be defined individually for every grid point
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

    def get_average(self, outside='none'):
        if outside is 'none':
            outside = self.field * 0
        return np.mean(self.field[np.logical_not(outside)])


class RandomField(Field):

    def __init__(self, name, nfx, nfy, nfz, pdfr, keep=False):
        Field.__init__(self, name)
        self.nfx = nfx
        self.nfy = nfy
        self.nfz = nfz
        self.pdfr = pdfr
        self.extendx = int(self.nfx.max())
        self.extendy = int(self.nfy.max())
        self.extendz = int(self.nfz.max())
        self.res_z = self.nfx.shape[0]
        self.res_y = self.nfy.shape[1]
        self.field = np.random.uniform(low=-self.pdfr, high=self.pdfr,
                                       size=(self.extendx * 2 + 1,
                                             self.extendy * 2 + 1 + self.res_z,
                                             self.extendz * 2 + 1 + self.res_y))
        self.startfield = self.field
        self.keep = keep

    def roll_one_timestep(self):
        self.field = np.roll(self.field, -1, axis=0)
        self.field[-1, :, :] = np.random.uniform(low=-self.pdfr, high=self.pdfr,
                                                 size=(self.extendy * 2 + 1 + self.res_z,
                                                       self.extendz * 2 + 1 + self.res_y))

    def fill_from_startfield(self, index):
        self.field = np.roll(self.field, -1, axis=0)
        self.field[-1, :, :] = self.startfield[index, :, :]

    def next_field(self, period=0, timestep=0, periodic=True):
        if period == 0 or not periodic:
            self.roll_one_timestep()
        else:
            s_field_ind = -(period - timestep - self.extendx) + self.extendx + 1
            if s_field_ind < 0:
                self.roll_one_timestep()
            else:
                self.fill_from_startfield(s_field_ind)


class CorrField(Field):

    def __init__(self, name, rand_field, corr_matrix, outside_flag, keep=False):
        Field.__init__(self, name)
        # initialising field
        self.field = np.zeros((0, rand_field.res_z, rand_field.res_y))
        self.outside = outside_flag
        self.corr_matrix = corr_matrix
        self.keep = keep
        self.rand_field = rand_field

    def add_field(self):
        new_field = np.zeros((1, self.rand_field.res_z, self.rand_field.res_y))
        for k in range(self.rand_field.res_y):
            for j in range(self.rand_field.res_z):
                if not self.outside[j, k]:
                    new_field[0, j, k] = np.sum(self.rand_field.field[:, j:j + 2 * int(self.rand_field.nfz[j, k]) + 1,
                                                k:k + 2 * int(self.rand_field.nfy[j, k]) + 1] *
                                                self.corr_matrix[k][j].matrix)
        self.field = np.concatenate((self.field, new_field), axis=0)


class OutField(Field):

    def __init__(self, name, rand_field, dt_split=1, keep=False):
        Field.__init__(self, name)
        self.keep = keep
        self.res_z = rand_field.res_z
        self.res_y = rand_field.res_y
        self.field = np.zeros((0, rand_field.res_z, rand_field.res_y))
        self.dt_split = dt_split

    def add_field(self, field):
        if self.dt_split > 1 and self.field.shape[0] > 0:
            self.add_intermediate_fields(self.field[-1, :, :], field, self.dt_split)
        self.field = np.concatenate((self.field, field.reshape((1, self.res_z,
                                                                self.res_y))), axis=0)

    def add_intermediate_fields(self, startfield, endfield, dt_split):
        endfield = endfield.reshape((1, self.res_z, self.res_y))
        startfield = startfield.reshape((1, self.res_z, self.res_y))
        for i in range(1, dt_split):
            self.field = np.concatenate((self.field, i / dt_split * endfield +
                                         (dt_split - i) / dt_split * startfield),
                                        axis=0)

    def add_last_fields(self):
        self.add_intermediate_fields(self.field[-1, :, :], self.field[0, :, :], self.dt_split)


class RstTransform:

    def __init__(self, dataarrays, res_y, res_z, fids=('u', 'v', 'w')):
        self.fids = fids
        self.rst = np.zeros((res_z, res_y, len(fids), len(fids)))
        self.rst_transform = np.zeros((res_z, res_y, len(fids), len(fids)))
        ids = [el.standardname for el in dataarrays]
        rst_comps = [[dataarrays[ids.index(id_switch(i, j, ids))].data
                      for i in fids]
                     for j in fids]

        for i in range(res_z):
            for j in range(res_y):
                self.rst[i, j, :, :] = np.array([[rst_comps[k][l][i, j]
                                                  for k in range(len(fids))]
                                                 for l in range(len(fids))])

        for i in range(res_z):
            for j in range(res_y):
                self.rst_transform[i, j, :, :] = npc.cholesky(self.rst[i, j, :, :])

    def adapt_2_profile(self, corrfields, outfields, profile):
        corr_names = [el.name for el in corrfields]
        out_names = [el.name for el in outfields]

        corr_inds = [corr_names.index(ident) for ident in self.fids]
        out_inds = [out_names.index(ident) for ident in self.fids]

        for i, eli in enumerate(self.fids):
            summe = 0
            for j, elj in enumerate(self.fids):
                summe = summe + self.rst_transform[:, :, i, j] * corrfields[corr_inds[j]].get_last_field()
            summe = summe + profile.get_data(eli)
            outfields[out_inds[i]].add_field(summe)


def id_switch(id1, id2, ids):
    if id1 + id2 in ids:
        return id1 + id2
    else:
        return id2 + id1
