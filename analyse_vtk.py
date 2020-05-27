import vtk
import numpy as np
import os
import fnmatch
import matplotlib.pyplot as plt
from scipy import interpolate
import prepare_traverse as preptrav
import h5py


class VtkField:
    def __init__(self, name, conv_check=-1, keep=False):
        self.conv_check = conv_check
        self.name = name
        self.count = 0
        self.inst_field = np.zeros(0)
        self.sum_field = np.zeros(0)
        self.conv_array = []
        self.time_array = []
        self.history = np.zeros(0)
        self.keep = keep
        self.ly = np.zeros(0)
        self.lz = np.zeros(0)
        self.lnx = np.zeros(0)
        self.lr = np.zeros(0)
        self.lphi = np.zeros(0)

    def add_field(self, field):
        if self.inst_field.size == 0:
            self.inst_field = field
            self.sum_field = field
            if self.keep:
                self.history = np.array([field])
        else:
            self.inst_field = field
            self.sum_field = self.sum_field + field
            if self.keep:
                self.history = np.concatenate((self.history, np.array([field])))
        self.count = self.count + 1
        if self.conv_check != -1:
            self.conv_array.append(self.sum_field[self.conv_check] / self.count)
            self.time_array.append(self.inst_field[self.conv_check])

    def get_avg_field(self):
        return self.sum_field / self.count

    def get_inst_field(self):
        return self.inst_field

    def get_hist_field(self, index, split=1, splitindex=0):
        if split == 1 or index == self.count-1:
            return self.history[index, :]
        elif split > 1 and splitindex < split:
            return self.history[index, :] * (split-splitindex)/split + self.history[index + 1, :] * splitindex / split

    def get_cov_field(self):
        return np.cov(self.history, rowvar=False)

    def calc_lengthscale(self, vtk_handle, res=1e-4, decay=1 / np.e):
        print('Processing Lengthscale Field {:s}'.format(self.name))
        cov_field = self.get_cov_field()
        self.ly = np.zeros(self.inst_field.shape)
        self.lz = np.zeros(self.inst_field.shape)
        self.lnx = np.zeros(self.inst_field.shape)
        coords = vtk_handle.get_cellpoints()
        # r = np.sqrt(coords[:, 1] ** 2 + coords[:, 2] ** 2)
        # index = np.logical_and(r > 0.2, r < 0.265)
        index = vtk_handle.get_type_inner_index()
        cy = coords[:, 1]
        cz = coords[:, 2]
        vy = np.arange(cy.min(), cy.max(), res)
        vz = np.arange(cz.min(), cz.max(), res)
        y, z = np.meshgrid(vy, vz)
        for i in range(cov_field.shape[0]):
            print('   Processing Index {:4d} / {:4d}'.format(i + 1, cov_field.shape[0]))
            if index[i]:
                cov_img = interpolate.griddata((cy, cz), cov_field[:, i] / cov_field[i, i],
                                               (y, z), fill_value=0)
                bw_img = cov_img > decay
                y_ind, z_ind = get_indizes(cy[i], cz[i], vy, vz)
                lny = np.sum(bw_img[z_ind, :]) / 2
                lnz = np.sum(bw_img[:, y_ind]) / 2

                self.ly[i] = lny * res
                self.lz[i] = lnz * res
                avg_val = self.get_avg_field()[i]
                tcorr = np.correlate(self.history[:, i] - avg_val, self.history[:, i] - avg_val, mode='full')
                tcorr_norm = tcorr / max(tcorr)
                lnx = (tcorr_norm[tcorr_norm > decay].size - 1) / 2

                self.lnx[i] = lnx

    def calc_lengthscale_2(self, vtk_handle, res=1e-4, decay=1/np.e, calcpoints=None):
        print('Processing Lengthscale Field {:s}'.format(self.name))
        cov_field = self.get_cov_field()
        self.ly = np.zeros(self.inst_field.shape)
        self.lz = np.zeros(self.inst_field.shape)
        self.lnx = np.zeros(self.inst_field.shape)
        self.lr = np.zeros(self.inst_field.shape)
        self.lphi = np.zeros(self.inst_field.shape)
        r, phi = vtk_handle.get_polar_cellpoints()
        r = np.squeeze(r)
        phi = np.squeeze(phi)
        index = vtk_handle.get_type_inner_index()
        r_min = r.min()
        r_max = r.max()
        rvect = np.arange(r_min, r_max+res, res)
        phi_res = res / ((r_min + r_max) / 2)
        phivect = np.arange(phi.min(), phi.max() + phi_res, phi_res)
        rg, phig = np.meshgrid(rvect, phivect)
        if calcpoints is None:
            calcpoints = range(cov_field.shape[0])
        for i in calcpoints:
            print('   Processing Field {:s} | Index {:4d} / {:4d}'.format(self.name, i+1, cov_field.shape[0]))
            if index[i]:
                cov_img = interpolate.griddata((r,phi), cov_field[:,i] / cov_field[i,i], (rg,phig),
                                               fill_value=0)
                bw_img = cov_img > decay
                r_ind, phi_ind = get_indizes(r[i],phi[i],rvect,phivect)
                lnr = np.sum(bw_img[phi_ind, :]) / 2
                lnphi = np.sum(bw_img[:,r_ind]) / 2
                self.lr[i] = lnr * res
                self.lphi[i] = lnphi * phi_res * rvect[r_ind]

                avg_val = self.get_avg_field()[i]
                tcorr = np.correlate(self.history[:, i] - avg_val, self.history[:, i] - avg_val, mode='full')
                tcorr_norm = tcorr / max(tcorr)
                lnx = (tcorr_norm[tcorr_norm > decay].size - 1) / 2

                self.lnx[i] = lnx




def get_indizes(x, y, x_array, y_array):
    x_ind = np.nanargmin((x_array - x) ** 2)
    y_ind = np.nanargmin((y_array - y) ** 2)
    return x_ind, y_ind


class CorrField:
    def __init__(self, field1, field2, conv_check=-1):
        self.name = field1.name + field2.name
        self.field1 = field1
        self.field2 = field2
        self.sum_field = np.zeros(0)
        self.count = 0
        self.conv_check = conv_check
        self.conv_array = []

    def add_field(self):
        if self.sum_field.size == 0:
            self.sum_field = self.field1.inst_field * self.field2.inst_field
        else:
            self.sum_field = self.sum_field + self.field1.inst_field * self.field2.inst_field
        self.count = self.count + 1
        if self.conv_check != -1:
            self.conv_array.append(self.sum_field[self.conv_check] / self.count)

    def get_avg_field(self):
        return self.sum_field / self.count

    def get_cov_field(self):
        return self.get_avg_field() - self.field1.get_avg_field() * self.field2.get_avg_field()


class FieldHandler:
    def __init__(self, ids, conv_check=-1, keep=False):
        self.field_list = [VtkField(name, conv_check, keep=keep) for name in ids]
        self.corr_list = [CorrField(el1, el2, conv_check) for el1 in self.field_list for el2 in self.field_list]

    def get_field_names(self):
        return [el.name for el in self.field_list]

    def get_field_names_corr(self):
        return [el.name for el in self.corr_list]

    def get_field_el(self, name):
        return self.field_list[self.get_field_names().index(name)]

    def get_field_el_corr(self, name):
        return self.corr_list[self.get_field_names_corr().index(name)]

    def add_fields(self, vtk_handle):
        for el in self.field_list:
            el.add_field(vtk_handle.get_datas(el.name))
        for el in self.corr_list:
            el.add_field()

    def get_count(self):
        return self.field_list[0].count


class VtkHandler:
    def __init__(self):
        self.reader = vtk.vtkGenericDataObjectReader()
        self.reader.ReadAllFieldsOn()
        self.reader.ReadAllVectorsOn()
        self.reader.ReadAllScalarsOn()

    def add_file(self, filename):
        self.reader.SetFileName(filename)
        self.reader.Update()

    def get_celldata(self):
        return self.reader.GetOutput().GetCellData()

    def get_pointdata(self):
        return self.reader.GetOutput().GetPoints()

    def get_cellpoints(self, scale=1):
        return self.get_cellpoints_orig(scale)[self.get_order(), :]

    def get_cellpoints_orig(self, scale=1):
        cellfilt = vtk.vtkCellCenters()
        cellfilt.SetInputData(self.reader.GetOutput())
        cellfilt.Update()
        out = np.array(cellfilt.GetOutput().GetPoints().GetData())
        out[:, 1:3] = out[:, 1:3] * scale
        return out

    def get_order(self):
        points = self.get_cellpoints_orig()
        return np.lexsort(points.T)

    def get_polar_cellpoints(self, scale=1, index=None):
        coords = self.get_cellpoints(scale)
        r = np.sqrt(coords[:, 1] ** 2 + coords[:, 2] ** 2)
        phi = np.arctan(-coords[:, 1] / coords[:, 2])
        return r[index], phi[index]

    def get_inner_index(self, min_rad=0.2, max_rad=0.265):
        r, phi = self.get_polar_cellpoints()
        return np.squeeze(np.logical_and(r > min_rad, r < max_rad))

    def get_type_index(self, typ=1):
        return self.get_datas('type_id') == typ

    def get_type_inner_index(self, min_rad=0.2, max_rad=0.265, typ=1):
        return np.logical_and(self.get_inner_index(min_rad, max_rad), self.get_type_index(typ))

    def get_datas(self, name):
        if name in ['u', 'v', 'w']:
            modname = 'velocity'
        else:
            modname = name
        scalnames = [self.get_celldata().GetArrayName(i) for i in range(self.reader.GetNumberOfScalarsInFile())]
        data = np.array(self.get_celldata().GetArray(scalnames.index(modname)))
        if name in ['u', 'v', 'w']:
            return data[self.get_order(), ['u', 'v', 'w'].index(name)]
        else:
            return data[self.get_order()]


class PostProcessor:
    def __init__(self, directory, ids, conv_check=-1, keep=False):
        self.vtk_handle = VtkHandler()
        self.fields = FieldHandler(ids, conv_check, keep=keep)
        self.directory = directory
        self.fr = None

    def process_start_stop(self, start, stop):
        for i in range(start, stop, 10):
            filename = '/home_elbe/goertz/mnt_hdd1/data/cutplane_1_Not4FluxCalc_0000{:05d}.vtk'.format(i)
            self.vtk_handle.add_file(filename)
            self.fields.add_fields(self.vtk_handle)

    def process_directory(self, nr_files=None, debug=False):
        if debug:
            fig = plt.figure(figsize=(12,8))
            fig_u = plt.figure(figsize=(12,8))
        files = sorted(os.listdir(self.directory))
        files = fnmatch.filter(files, '*.vtk')
        for i, f in enumerate(files):
            print('Processing File {:s}, {:4d}/{:4d}'.format(f, i + 1, len(files)))
            filename = self.directory + '/' + f
            self.vtk_handle.add_file(filename)
            self.fields.add_fields(self.vtk_handle)
            if nr_files is not None:
                if i > nr_files:
                    break
            if debug and (i+1) % 50 == 0:
                xyz = self.vtk_handle.get_cellpoints()
                # fig = plt.figure()
                # ax = fig.add_subplot(111)
                fig.clf()
                ax = fig.add_subplot(111)
                h = ax.scatter(xyz[:,1],xyz[:,2],8,self.fields.get_field_el('T').get_avg_field())
                fig.colorbar(h, ax=ax)
                ax.set_title(str(i+1))
                fig.savefig('/home_elbe/goertz/mnt_hdd1/data/convergence/convT{:06d}.png'.format(i+1))
                fig_u.clf()
                ax_u = fig_u.add_subplot(111)
                h_u = ax_u.scatter(xyz[:,1],xyz[:,2],8,self.fields.get_field_el('u').get_avg_field())
                fig_u.colorbar(h_u, ax=ax_u)
                ax_u.set_title(str(i+1))
                fig_u.savefig('/home_elbe/goertz/mnt_hdd1/data/convergence/conv_u{:06d}.png'.format(i+1))
                #plt.draw()
                #plt.waitforbuttonpress()

    def process_lengthscales(self, ids='all', res=1e-4, decay=1/np.e):
        if ids is 'all':
            for f in self.fields.field_list:
                f.calc_lengthscale(self.vtk_handle, res, decay)
        else:
            for name in ids:
                self.fields.get_field_el(name).calc_lengthscale(self.vtk_handle, res, decay)

    def summary_plot_velo(self):
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        h1 = ax1.scatter(self.vtk_handle.get_cellpoints()[:, 1],
                         self.vtk_handle.get_cellpoints()[:, 2],
                         10,
                         self.fields.get_field_el('u').get_avg_field())
        quiv1 = ax1.quiver(self.vtk_handle.get_cellpoints()[:, 1],
                           self.vtk_handle.get_cellpoints()[:, 2],
                           self.fields.get_field_el('v').get_avg_field(),
                           self.fields.get_field_el('w').get_avg_field(),
                           scale=2000, headlength=3, headwidth=2)
        ax1.set_aspect('equal', 'box')
        fig.colorbar(h1, ax=ax1)
        ax2 = fig.add_subplot(122)
        h2 = ax2.scatter(self.vtk_handle.get_cellpoints()[:, 1],
                         self.vtk_handle.get_cellpoints()[:, 2],
                         10,
                         0.5 * (self.fields.corr_list[0][0].get_cov_field() +
                                self.fields.corr_list[1][1].get_cov_field() +
                                self.fields.corr_list[2][2].get_cov_field()))
        ax2.set_aspect('equal', 'box')
        fig.colorbar(h2, ax=ax2)
        return h1, h2, quiv1

    def add_file_reader_instant(self, timestep, splitindex, split, scale=1, filename='generic'):
        coords = self.vtk_handle.get_cellpoints(scale)
        index = self.vtk_handle.get_type_inner_index()
        headerfields = ['x', 'y', 'z']
        data = coords
        for el in self.fields.field_list:
            headerfields.append(el.name)
            hist_field = el.get_hist_field(timestep, split, splitindex)
            data = np.concatenate((data, hist_field.reshape(hist_field.size, 1)), axis=1)
        data = data[index, :]
        self.add_file_reader(headerfields, data, filename)

    def add_file_reader_avg(self, scale=1, filename='generic'):
        coords = self.vtk_handle.get_cellpoints(scale)
        index = self.vtk_handle.get_type_inner_index()
        headerfields = ['x', 'y', 'z']
        data = coords
        for el in self.fields.field_list:
            headerfields.append(el.name)
            avg_field = el.get_avg_field()
            data = np.concatenate((data, avg_field.reshape(avg_field.size, 1)), axis=1)
        for el in self.fields.corr_list:
            headerfields.append(el.name)
            cov_field = el.get_cov_field()
            data = np.concatenate((data, cov_field.reshape(cov_field.size, 1)), axis=1)
        data = data[index, :]
        self.add_file_reader(headerfields, data, filename)

    def add_file_reader(self, headerfields, data, filename):
        self.fr = preptrav.FileReader(filename)
        self.fr.add_data(headerfields, data)
        self.fr.extract_column_data()
        self.fr.expand_rot_periodic()
        self.fr.calc_angles()

    def change_x_value(self, x_val=None):
        if x_val is not None:
            self.fr.change_data_to_constant('x', x_val)

    def save_traverses(self, filename, mode='avgcfx', scale=1.0, min_phi=None, max_phi=None, split=1, xval=None):
        if mode.startswith('avg'):
            self.add_file_reader_avg(scale, filename)
            if min_phi is not None and max_phi is not None:
                self.fr.set_phi_range(min_phi, max_phi)
            self.change_x_value(xval)
            if mode.endswith('cfx'):
                self.fr.write_cfx_file(0)
            elif mode.endswith('prf'):
                self.fr.write_prf_file(0)
        elif mode == 'inst':
            os.makedirs(filename, exist_ok=True)
            for i in range(self.fields.get_count()):
                for j in range(split):
                    self.add_file_reader_instant(i, j, split, scale, filename + '/inlet.dat{:d}'.format(i*split+j+1))
                    if min_phi is not None and max_phi is not None:
                        self.fr.set_phi_range(min_phi, max_phi)
                    self.change_x_value(xval)
                    self.fr.write_inlet_file()
        elif mode == 'hdf5':
            f = h5py.File(filename+ 'Inlet.hdf5', 'w')
            main = f.create_group("main")
            main.attrs.create("Number of Snapshots",self.fields.get_count()*split+split+1)
            main.attrs.create("Number of Points",len(self.fr.data_columns[0].data))
            main.attrs.create("Number of Variables",7)
            main.attrs.create("Variables", np.string_("phi, r, pt, Tt, whirl, Pitch"))
            main.attrs.create("Timestep width", 1.001e-5/split)
            f.close()
            for i in range(self.fields.get_count()):
                for j in range(split):
                    self.count = i * split + j + 1
                    self.add_file_reader_instant(i, j, split, scale, filename + 'Inlet.hdf5')
                    if min_phi is not None and max_phi is not None:
                        self.fr.set_phi_range(min_phi, max_phi)
                    self.change_x_value(xval)
                    self.fr.write_hdf5(self.count)
            return


    def get_polar_field(self, typ='avg'):
        if typ is 'avg':
            v_feld = np.array([self.fields.get_field_el('v').get_avg_field()])
            w_feld = np.array([self.fields.get_field_el('w').get_avg_field()])
        elif typ is 'hist':
            v_feld = self.fields.get_field_el('v').history
            w_feld = self.fields.get_field_el('w').history

        r, phi = self.vtk_handle.get_polar_cellpoints()

        vrad = w_feld * np.cos(phi) - v_feld * np.sin(phi)
        vtan = -v_feld * np.cos(phi) - w_feld * np.sin(phi)

        return np.squeeze(vrad), np.squeeze(vtan)


def autocorr(x):
    x_avg = np.mean(x)
    cov_x = np.correlate(x - x_avg, x - x_avg, mode='full')
    cov_norm = cov_x / np.max(cov_x)
    return cov_norm[len(x) - 1:]


def meandiffhistory(histfield):
    # timestep history on index 0
    steps, _ = np.meshgrid(np.arange(histfield.shape[0])+1, np.arange(histfield.shape[1]), indexing='ij')
    convhist = np.cumsum(histfield, axis=0)/steps
    refvals, _ = np.meshgrid(convhist[-1,:], np.arange(convhist.shape[0]))
    normdiff = (convhist-refvals)/refvals
    return np.mean(np.abs(normdiff),axis=1)


def get_convhist(histfield):
    steps, _ = np.meshgrid(np.arange(histfield.shape[0])+1, np.arange(histfield.shape[1]), indexing='ij')
    return np.cumsum(histfield, axis=0)/steps


def get_corrhistory(histfield1, histfield2):
    covhist = get_convhist(histfield1*histfield2)-get_convhist(histfield1)*get_convhist(histfield2)
    refvals, _ = np.meshgrid(covhist[-1,:], np.arange(covhist.shape[0]))
    normdiff = (covhist-refvals)/refvals
    return np.mean(np.abs(normdiff),axis=1)


#pp = PostProcessor('/mnt/raid/Working/Precise/data_small/', ['u', 'v', 'w', 'T', 'P'], conv_check=10855, keep=True)
# pp.process_start_stop(8560,8660)
pp.process_directory()
pp.save_traverses('/mnt/raid/Working/Precise/Output/mean.csv',
                  mode='avgcfx',
                  scale=1.0086,
                  min_phi=-15,
                  max_phi=15,
                  xval=-0.053612)
pp.save_traverses('/mnt/raid/Working/Precise/Output/mean.csv',
                  mode='avgprf',
                  scale=1.0086,
                  min_phi=-15,
                  max_phi=15,
                  xval=-0.053612)
# pp.save_traverses('/mnt/raid/Working/Precise/Output/instat/',
#                   mode='inst',
#                   scale=1.0086,
#                   min_phi=-15,
#                   max_phi=15,
#                   split=33,
#                   xval=-0.053612)
pp.save_traverses('/mnt/raid/Working/Precise/Output/',
                  mode='hdf5',
                  scale=1.0086,
                  min_phi=-15,
                  max_phi=15,
                  split=32,
                  xval=-0.053612)


# pp.process_lengthscales()
# pp.fields.field_list[0].calc_lengthscale(pp.vtk_handle)
# pp.fields.field_list[1].calc_lengthscale(pp.vtk_handle)
# pp.fields.field_list[2].calc_lengthscale(pp.vtk_handle)
# pp.fields.field_list[3].calc_lengthscale(pp.vtk_handle)
# h1, h2, q1 = pp.summary_plot_velo()
