import numpy as np
import fnmatch


class DataColumn:

    def __init__(self, headerfield, data, input_col=False):
        self.headerfield = headerfield
        self.standardname = 'other'
        self.input_col = input_col
        self.data = data
        self.filterfield = False
        self.cart_field = False
        self.fill_standard_name()

    def fill_standard_name(self):
        if self.headerfield in ['X [ m ]', 'x [m]', 'x']:
            self.standardname = 'x'
            self.filterfield = True
        elif self.headerfield in ['Y [ m ]', 'y [m]', 'y']:
            self.standardname = 'y'
            self.filterfield = True
            self.cart_field = True
        elif self.headerfield in ['Z [ m ]', 'z [m]', 'z']:
            self.standardname = 'z'
            self.filterfield = True
            self.cart_field = True
        elif self.headerfield in ['Velocity u [ m s^-1 ]', 'Vax [m s^-1]', 'u', 'U']:
            self.standardname = 'u'
            self.filterfield = True
        elif self.headerfield in ['Velocity v [ m s^-1 ]', 'V [m s^-1]', 'v']:
            self.standardname = 'v'
            self.filterfield = True
            self.cart_field = True
        elif self.headerfield in ['Velocity w [ m s^-1 ]', 'W [m s^-1]', 'w']:
            self.standardname = 'w'
            self.filterfield = True
            self.cart_field = True
        elif self.headerfield in ['Vrad [m s^-1]', 'V rad']:
            self.standardname = 'vrad'
        elif self.headerfield in ['Vtan [m s^-1]', 'V tan']:
            self.standardname = 'vtan'
        elif self.headerfield in ['Turbulence Kinetic Energy [ m^2 s^-2 ]', 'k [m^2 s^-2]', 'k']:
            self.standardname = 'k'
            self.filterfield = True
        elif self.headerfield in ['Turbulence Eddy Dissipation [ m^2 s^-3 ]', 'e']:
            self.standardname = 'e'
            self.filterfield = True
        elif self.headerfield in ['w [s^-1]', 'sdr']:
            self.standardname = 'sdr'
            self.filterfield = True
        elif self.headerfield in ['Radius']:
            self.standardname = 'radius'
        elif self.headerfield in ['Phi_rad']:
            self.standardname = 'phi_rad'
        elif self.headerfield in ['Phi_deg']:
            self.standardname = 'phi_deg'
        elif self.headerfield in ['uu', 'Statistical Reynolds Stress uu [ m^2 s^-2 ]']:
            self.standardname = 'uu'
        elif self.headerfield in ['vv', 'Statistical Reynolds Stress vv [ m^2 s^-2 ]']:
            self.standardname = 'vv'
        elif self.headerfield in ['ww', 'Statistical Reynolds Stress ww [ m^2 s^-2 ]']:
            self.standardname = 'ww'
        elif self.headerfield in ['uv', 'Statistical Reynolds Stress uv [ m^2 s^-2 ]']:
            self.standardname = 'uv'
        elif self.headerfield in ['uw', 'Statistical Reynolds Stress uw [ m^2 s^-2 ]']:
            self.standardname = 'uw'
        elif self.headerfield in ['vw', 'Statistical Reynolds Stress vw [ m^2 s^-2 ]']:
            self.standardname = 'vw'
        else:
            self.standardname = self.headerfield


class FileReader:

    def __init__(self, filename):
        self.filename = filename
        self.headerfields = []
        self.data = []
        self.data_columns = []
        self.phi_range = [0, 0]
        self.header = []

    def open_file(self):
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

            self.headerfields = [i.strip() for i in headerfields]
            self.data = np.loadtxt(self.filename, skiprows=startline, delimiter=delim)
            self.header = lines[:startline - 1]

    def write_prf_file(self):
        if not self.phi_range == [0, 0]:
            phi_ind = self.get_standardnames().index('phi_deg')

            out_ind = np.logical_or(self.data_columns[phi_ind].data < self.phi_range[0] - 1.5,
                                    self.data_columns[phi_ind].data > self.phi_range[1] + 1.5)
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

    def write_cfx_file(self):
        with open(self.filename[:-4] + '_new.csv', 'w') as out_file:
            out_file.writelines(self.header)
            headstr = ''.join([el.headerfield + ', ' for el in self.data_columns])[:-2] + '\n'
            out_file.write(headstr)

            for i in range(self.data_columns[0].data.shape[0]):
                dataline = ''.join(['{:8e}, '.format(el.data[i])
                                    for el in self.data_columns])[:-2] + '\n'
                out_file.write(dataline)

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
        if 'vtan' not in self.get_standardnames():
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
                el.data = np.concatenate((el.data - phi_period + delta_phi, el.data + delta_phi,
                                          el.data + phi_period + delta_phi), axis=0)
            elif el.standardname == 'phi_rad':
                el.data = np.concatenate((el.data - (phi_period - delta_phi) * np.pi / 180,
                                          el.data + delta_phi * np.pi / 180,
                                          el.data + (phi_period + delta_phi) * np.pi / 180))
            elif el.cart_field:
                pass
            else:
                el.data = np.concatenate((el.data, el.data, el.data), axis=0)

        self.make_cart_coord()
        self.make_cart_velo()


# fr = FileReader("/home_elbe/goertz/mnt_hdd4/traversen/ISSUE09_from_sim_P.csv")
# fr = FileReader("/home_elbe/goertz/profile_3.prf")
# fr.open_file()
# fr.extract_column_data()
# fr.expand_rot_periodic(-5.625)
# fr.write_prf_file()
# fr.write_cfx_file()
