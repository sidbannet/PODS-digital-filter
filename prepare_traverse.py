import numpy as np
import fnmatch
import h5py


class DataColumn:

    def __init__(self, headerfield, data, input_col=False, ignoreunknown=True):
        self.headerfield = headerfield
        self.standardname = 'other'
        self.unit = '-'
        self.input_col = input_col
        self.data = data
        self.filterfield = False
        self.cart_field = False
        self.csv_header = None
        self.csvfield = True
        self.ignore = ignoreunknown
        self.fill_standard_name()

    def fill_standard_name(self):
        if self.headerfield in ['X [ m ]', 'x [m]', 'x']:
            self.standardname = 'x'
            self.unit = 'm'
            self.filterfield = True
            self.csv_header = 'x [ m ]'
        elif self.headerfield in ['Y [ m ]', 'y [m]', 'y']:
            self.standardname = 'y'
            self.unit = 'm'
            self.filterfield = True
            self.cart_field = True
            self.csv_header = 'y [ m ]'
        elif self.headerfield in ['Z [ m ]', 'z [m]', 'z']:
            self.standardname = 'z'
            self.unit = 'm'
            self.filterfield = True
            self.cart_field = True
            self.csv_header = 'z [ m ]'
        elif self.headerfield in ['Velocity u [ m s^-1 ]', 'Vax [m s^-1]', 'u', 'U',
                                  'Velocity u.Trnavg [ m s^-1 ]']:
            self.standardname = 'u'
            self.unit = 'm s^{-1}'
            self.filterfield = True
            self.csv_header = 'Velocity u [ m s^-1 ]'
        elif self.headerfield in ['Velocity v [ m s^-1 ]', 'V [m s^-1]', 'v',
                                  'Velocity v.Trnavg [ m s^-1 ]']:
            self.standardname = 'v'
            self.unit = 'm s^{-1}'
            self.filterfield = True
            self.cart_field = True
            self.csv_header = 'Velocity v [ m s^-1 ]'
        elif self.headerfield in ['Velocity w [ m s^-1 ]', 'W [m s^-1]', 'w',
                                  'Velocity w.Trnavg [ m s^-1 ]']:
            self.standardname = 'w'
            self.unit = 'm s^{-1}'
            self.filterfield = True
            self.cart_field = True
            self.csv_header = 'Velocity w [ m s^-1 ]'
        elif self.headerfield in ['Vrad [m s^-1]', 'V rad', 'vrad']:
            self.standardname = 'vrad'
            self.unit = 'm s^{-1}'
        elif self.headerfield in ['Vtan [m s^-1]', 'V tan', 'vtan']:
            self.standardname = 'vtan'
            self.unit = 'm s^{-1}'
        elif self.headerfield in ['Turbulence Kinetic Energy [ m^2 s^-2 ]', 'k [m^2 s^-2]', 'k']:
            self.standardname = 'k'
            self.unit = 'm^2 s^{-2}'
            self.filterfield = True
        elif self.headerfield in ['Turbulence Eddy Dissipation [ m^2 s^-3 ]', 'e']:
            self.standardname = 'e'
            self.unit = 'm^2 s^{-3}'
            self.filterfield = True
        elif self.headerfield in ['w [s^-1]', 'sdr']:
            self.standardname = 'sdr'
            self.unit = 's^{-1}'
            self.filterfield = True
        elif self.headerfield in ['Radius']:
            self.standardname = 'radius'
            self.unit = 'm'
        elif self.headerfield in ['Phi_rad']:
            self.standardname = 'phi_rad'
            self.unit = 'rad'
        elif self.headerfield in ['Phi_deg']:
            self.standardname = 'phi_deg'
            self.unit = 'deg'
        elif self.headerfield in ['uu', 'Statistical Reynolds Stress uu [ m^2 s^-2 ]']:
            self.standardname = 'uu'
            self.filterfield = True
            self.unit = 'm^2 s^{-2}'
            self.csv_header = 'Statistical Reynolds Stress uu [ m^2 s^-2 ]'
        elif self.headerfield in ['vv', 'Statistical Reynolds Stress vv [ m^2 s^-2 ]']:
            self.standardname = 'vv'
            self.filterfield = True
            self.unit = 'm^2 s^{-2}'
            self.csv_header = 'Statistical Reynolds Stress vv [ m^2 s^-2 ]'
        elif self.headerfield in ['ww', 'Statistical Reynolds Stress ww [ m^2 s^-2 ]']:
            self.standardname = 'ww'
            self.filterfield = True
            self.unit = 'm^2 s^{-2}'
            self.csv_header = 'Statistical Reynolds Stress ww [ m^2 s^-2 ]'
        elif self.headerfield in ['uv', 'Statistical Reynolds Stress uv [ m^2 s^-2 ]']:
            self.standardname = 'uv'
            self.filterfield = True
            self.unit = 'm^2 s^{-2}'
            self.csv_header = 'Statistical Reynolds Stress uv [ m^2 s^-2 ]'
        elif self.headerfield in ['uw', 'Statistical Reynolds Stress uw [ m^2 s^-2 ]']:
            self.standardname = 'uw'
            self.filterfield = True
            self.unit = 'm^2 s^{-2}'
            self.csv_header = 'Statistical Reynolds Stress uw [ m^2 s^-2 ]'
        elif self.headerfield in ['vw', 'Statistical Reynolds Stress vw [ m^2 s^-2 ]']:
            self.standardname = 'vw'
            self.unit = 'm^2 s^{-2}'
            self.filterfield = True
            self.csv_header = 'Statistical Reynolds Stress vw [ m^2 s^-2 ]'
        elif self.headerfield in ['TT']:
            self.standardname = 'TT'
            self.filterfield = True
        elif self.headerfield in ['uT']:
            self.standardname = 'uT'
            self.filterfield = False
        elif self.headerfield in ['vT']:
            self.standardname = 'vT'
            self.filterfield = False
        elif self.headerfield in ['wT']:
            self.standardname = 'wT'
            self.filterfield = False
        elif self.headerfield in ['TtTt']:
            self.standardname = 'TtTt'
            self.filterfield = True
        elif self.headerfield in ['uTt']:
            self.standardname = 'uTt'
            self.filterfield = False
        elif self.headerfield in ['vTt']:
            self.standardname = 'vTt'
            self.filterfield = False
        elif self.headerfield in ['wTt']:
            self.standardname = 'wTt'
            self.filterfield = False
        elif self.headerfield in ['Total Pressure [ Pa ]', 'Pt [ Pa ]', 'pt']:
            self.standardname = 'pt'
            self.unit = 'Pa'
        elif self.headerfield in ['P']:
            self.standardname = 'p'
            self.unit = 'Pa'
            self.csv_header = 'Pressure [ Pa ]'
        elif self.headerfield in ['T']:
            self.standardname = 'T'
            self.unit = 'K'
            self.filterfield = True
            self.csv_header = 'Temperature [ K ]'
        elif self.headerfield in ['Total Temperature [ K ]', 'Tt [ K ]', 'Tt']:
            self.standardname = 'Tt'
            self.unit = 'K'
            self.filterfield = True
        else:
            self.standardname = self.headerfield
            if self.ignore:
                self.csvfield = False

        if self.csv_header is None:
            self.csv_header = self.headerfield


class FileReader:
# this class creates a file reader object.

    def __init__(self, filename):
        self.filename = filename
        self.headerfields = []
        self.data = []
        self.data_columns = []
        self.phi_range = [0, 0]
        self.header = []

    def set_phi_range(self, min_phi, max_phi):
        self.phi_range = [min_phi, max_phi]

    def open_file(self):
        # open input file depending on data extension .csv, .prf, or 1d profile
        with open(self.filename, 'r') as dat_file:
            lines = dat_file.readlines()
            if fnmatch.fnmatch(self.filename, '*.csv'):
                idline = fnmatch.filter(lines, '[Data*')[0]
                startline = lines.index(idline) + 2
                headerline = lines[lines.index(idline) + 1]
                headerfields = headerline.strip().split(',')
                delim = ','
            elif fnmatch.fnmatch(self.filename, '*.prf'):
                idline = fnmatch.filter(lines, 'data,*')[0]
                startline = lines.index(idline) + 1
                headerline = lines[lines.index(idline)]
                headerfields = headerline.strip().split(',')[1:]
                delim = ','
            else:
                idline = fnmatch.filter(lines, '*[yUu]*[yUu]*')[0]
                startline = lines.index(idline) + 1
                headerfields = idline.split()
                delim = None

            # read header fields
            self.headerfields = [i.strip() for i in headerfields]

            # read data
            self.data = np.loadtxt(self.filename, skiprows=startline, delimiter=delim)

            # read header
            self.header = lines[:startline - 1]

    def add_data(self, headerfields, data):
        self.headerfields = headerfields
        self.data = data
        self.header = '[Name]\n' \
                      'BK Austritt mean\n' \
                      '\n' \
                      '[Spatial Fields]\n' \
                      'x,y,z\n' \
                      '\n' \
                      '[Data]\n'

    def change_data_to_constant(self, name, const):
        self.data_columns[self.get_standardnames().index(name)].data = self.data_columns[
            self.get_standardnames().index(name)].data * 0 + const

    def write_prf_file(self, extend=1.5):
        if not self.phi_range == [0, 0]:
            phi_ind = self.get_standardnames().index('phi_deg')

            out_ind = np.logical_or(self.data_columns[phi_ind].data < self.phi_range[0] - extend,
                                    self.data_columns[phi_ind].data > self.phi_range[1] + extend)
        else:
            out_ind = np.ones(self.data_columns[0].data.shape)

        with open(self.filename[:-4] + '_new.prf', 'w') as out_file:
            headstr = 'data,' + ''.join([el.standardname + ','
                                         for el in self.data_columns if el.filterfield])[:-1] + '\n'
            out_file.write(headstr)

            for i in range(self.data_columns[0].data.shape[0]):
                if not out_ind[i]:
                    dataline = ''.join(['{:8e}, '.format(el.data[i])
                                        for el in self.data_columns if el.filterfield])[:-2] + '\n'
                    out_file.write(dataline)

    def write_cfx_file(self, extend=1.5):
        if not self.phi_range == [0, 0]:
            phi_ind = self.get_standardnames().index('phi_deg')

            out_ind = np.logical_or(self.data_columns[phi_ind].data < self.phi_range[0] - extend,
                                    self.data_columns[phi_ind].data > self.phi_range[1] + extend)
        else:
            out_ind = np.ones(self.data_columns[0].data.shape)

        with open(self.filename[:-4] + '_new.csv', 'w') as out_file:
            out_file.writelines(self.header)
            headstr = ''.join([el.csv_header + ', ' for el in self.data_columns if el.csvfield])[:-2] + '\n'
            out_file.write(headstr)

            for i in range(self.data_columns[0].data.shape[0]):
                if not out_ind[i]:
                    dataline = ''.join(['{:8e}, '.format(el.data[i])
                                        for el in self.data_columns if el.csvfield])[:-2] + '\n'
                    out_file.write(dataline)

    def write_inlet_file(self):
        if not self.phi_range == [0, 0]:
            phi_ind = self.get_standardnames().index('phi_deg')

            out_ind = np.logical_or(self.data_columns[phi_ind].data < self.phi_range[0],
                                    self.data_columns[phi_ind].data > self.phi_range[1])
        else:
            out_ind = np.ones(self.data_columns[0].data.shape)

        with open(self.filename, 'w') as out_file:
            for i in range(self.data_columns[0].data.shape[0]):
                if not out_ind[i]:
                    out_file.write('{:.12f} {:.12f} {:.12f} {:.12f} {:.12f} {:.12f} {:.12f}\n'.format(
                        self.get_data('x')[i],
                        self.get_data('y')[i],
                        self.get_data('z')[i],
                        self.get_data('u')[i],
                        self.get_data('v')[i],
                        self.get_data('w')[i],
                        self.get_data('T')[i]
                    ))

    def write_hdf5(self, count):
        data = np.transpose(np.array((self.get_data('phi_deg'),
                self.get_data('radius'),
                self.get_data('p'),
                self.get_data('T'),
                self.get_data('whirl_deg'),
                self.get_data('pitch_deg'))))
        f = h5py.File(self.filename, 'r+')
        # create 'main' group
        main = f.require_group("main")
        data = main.create_dataset("snapshot_"+str(count), data = data)


        f.close()

    # write data columns into array data_columns
    def extract_column_data(self):
        for i, fieldn in enumerate(self.headerfields):
            self.data_columns.append(DataColumn(fieldn, self.data[:, i], input_col=True))

    def insert_column(self, fieldn, data):
        self.data_columns.append(DataColumn(fieldn, data))

    def get_standardnames(self):
        return [i.standardname for i in self.data_columns]

    def get_data(self, name):
        return self.data_columns[self.get_standardnames().index(name)].data

    def get_headerfields(self):
        return [i.headerfield for i in self.data_columns]

    def make_pol_coord(self):
        y_ind = self.get_standardnames().index('y')
        z_ind = self.get_standardnames().index('z')

        radius = np.sqrt(self.data_columns[y_ind].data ** 2 + self.data_columns[z_ind].data ** 2)
        phi_rad = np.arctan(-self.data_columns[y_ind].data / self.data_columns[z_ind].data)

        self.data_columns.append(DataColumn('Radius', radius))
        self.data_columns.append(DataColumn('Phi_rad', phi_rad))
        self.data_columns.append(DataColumn('Phi_deg', phi_rad * 180 / np.pi))

    def make_pol_velo(self):
        v_ind = self.get_standardnames().index('v')
        w_ind = self.get_standardnames().index('w')
        phi_ind = self.get_standardnames().index('phi_rad')

        vrad = self.data_columns[w_ind].data * np.cos(self.data_columns[phi_ind].data) - \
               self.data_columns[v_ind].data * np.sin(self.data_columns[phi_ind].data)
        vtan = - self.data_columns[v_ind].data * np.cos(self.data_columns[phi_ind].data) - \
               self.data_columns[w_ind].data * np.sin(self.data_columns[phi_ind].data)

        self.data_columns.append(DataColumn('V rad', vrad))
        self.data_columns.append(DataColumn('V tan', vtan))

    def make_cart_coord(self):
        phi_ind = self.get_standardnames().index('phi_rad')
        rad_ind = self.get_standardnames().index('radius')

        y_new = self.data_columns[rad_ind].data * -np.sin(self.data_columns[phi_ind].data)
        if 'y' in self.get_standardnames():
            self.data_columns[self.get_standardnames().index('y')].data = y_new
        else:
            self.data_columns.append(DataColumn('y [m]', y_new))

        z_new = self.data_columns[rad_ind].data * np.cos(self.data_columns[phi_ind].data)
        if 'z' in self.get_standardnames():
            self.data_columns[self.get_standardnames().index('z')].data = z_new
        else:
            self.data_columns.append(DataColumn('z [m]', z_new))

    def make_cart_velo(self):
        phi_ind = self.get_standardnames().index('phi_rad')
        v_tan_ind = self.get_standardnames().index('vtan')
        v_rad_ind = self.get_standardnames().index('vrad')

        v_new = - self.data_columns[v_tan_ind].data * \
                np.cos(self.data_columns[phi_ind].data) - \
                self.data_columns[v_rad_ind].data * \
                np.sin(self.data_columns[phi_ind].data)
        if 'v' in self.get_standardnames():
            self.data_columns[self.get_standardnames().index('v')].data = v_new
        else:
            self.data_columns.append(DataColumn('V [m s^-1]', v_new))

        w_new = self.data_columns[v_rad_ind].data * \
                np.cos(self.data_columns[phi_ind].data) - \
                self.data_columns[v_tan_ind].data * \
                np.sin(self.data_columns[phi_ind].data)
        if 'w' in self.get_standardnames():
            self.data_columns[self.get_standardnames().index('w')].data = w_new
        else:
            self.data_columns.append(DataColumn('W [m s^-1]', w_new))

    def expand_rot_periodic(self, delta_phi=0):
        if 'phi_deg' not in self.get_standardnames():
            self.make_pol_coord()
        if 'vtan' not in self.get_standardnames() and 'v' in self.get_standardnames()\
                and 'w' in self.get_standardnames():
            self.make_pol_velo()

        # self.data_columns = [el for el in self.data_columns if not el.cart_field]
        rad_ind = self.get_standardnames().index('radius')
        phi_ind = self.get_standardnames().index('phi_deg')

        phi_row = self.data_columns[phi_ind].data[np.round(self.data_columns[rad_ind].data, 4) ==
                                                  np.round(self.data_columns[rad_ind].data.max(), 4)]

        phi_period = np.round(phi_row.max() - phi_row.min(), 1)

        self.phi_range = (self.data_columns[phi_ind].data.min(),
                          self.data_columns[phi_ind].data.max())

        for el in self.data_columns:
            if el.standardname == 'phi_deg':
                el.data = np.concatenate((el.data - phi_period, el.data,
                                          el.data + phi_period), axis=0)
            elif el.standardname == 'phi_rad':
                el.data = np.concatenate((el.data - phi_period * np.pi / 180,
                                          el.data,
                                          el.data + phi_period * np.pi / 180))
            elif el.cart_field:
                pass
            else:
                el.data = np.concatenate((el.data, el.data, el.data), axis=0)

        self.make_cart_coord()
        if 'vtan' in self.get_standardnames() and 'vrad' in self.get_standardnames():
            self.make_cart_velo()

    def calc_angles(self):
        whirl_angle_rad = -np.arctan(self.data_columns[self.get_standardnames().index('vtan')].data/self.data_columns[self.get_standardnames().index('u')].data)
        whirl_angle_deg = whirl_angle_rad * 180 / np.pi
        pitch_angle_rad = np.arctan(self.data_columns[self.get_standardnames().index('vrad')].data/self.data_columns[self.get_standardnames().index('u')].data)
        pitch_angle_deg = pitch_angle_rad * 180 / np.pi
        self.data_columns.append(DataColumn('whirl_rad', whirl_angle_rad))
        self.data_columns.append(DataColumn('pitch_rad', pitch_angle_rad))
        self.data_columns.append(DataColumn('whirl_deg', whirl_angle_deg))
        self.data_columns.append(DataColumn('pitch_deg', pitch_angle_deg))

# fr = FileReader("/home_elbe/goertz/mnt_hdd4/traversen/ISSUE09_from_sim_P.csv")
# fr = FileReader("/home_elbe/goertz/profile_3.prf")
# fr.open_file()
# fr.extract_column_data()
# fr.expand_rot_periodic(-5.625)
# fr.write_prf_file()
# fr.write_cfx_file()
