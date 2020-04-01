import numpy as np
from numpy import linalg
import math
import nplotlib as nplt
import h5py


class POD:
    def __init__(self, digfilter, num_modes, outputfolder, correct_for_cell_volumes):
        self.num_modes_trunc = num_modes
        self.outputfolder = outputfolder+'/PODFS/'
        self.correct_for_cell_volumes = correct_for_cell_volumes
        self.cell_volume = 1
        self.A = digfilter.A
        self.dt = digfilter.dt
        self.mean_field = digfilter.mean_field
        self.num_components = digfilter.num_components
        self.fluct_comps = digfilter.fluct_comps_pod
        self.num_points = digfilter.num_points
        self.num_snapshots = digfilter.num_snapshots
        self.tol_CN = 1e-15
        self.coordinates = digfilter.coordinates
        self.coord_comps = digfilter.coord_comps
        self.num_coord_comps = digfilter.num_coord_comps

        print('\n==============================================================')
        print('Calculating POD modes ...')

        print('\n   Calculating covariance matrix ...')

        # calculate covariance matrix
        self.calculate_correlation_matrix()

        print('\n   Solving eigenvalue problem ...')
        # solve eigenvalue problem
        self.eigenvalue_problem()

        # sort eigenvalues
        self.sort_eigenvalues()

        # determine valid modes with positive energy
        self.valid_temporal_modes()

        # scaling of temporal modes
        self.scale_temp_modes()

        # calculate spatial modes
        self.spatial_modes()

        print('\n==============================================================')
        print('Writing POD output ...')
        # write POD output
        self.pod_output()




    def calculate_correlation_matrix(self):
        # allocate covariance matrix
        self.C = np.array(np.zeros((self.num_snapshots, self.num_snapshots), dtype=np.float64))
        if self.correct_for_cell_volumes:
            for i in range(0, self.num_snapshots):
                print("      Correlations with snapshot ", i + 1, " of ", self.num_snapshots)
            for j in range(0, self.num_snapshots):
                if (j >= i):
                    for k in range(0, self.num_components):
                        self.C[i, j] = self.C[i, j] + np.dot(self.A[k * self.num_points:(k + 1) * self.num_points, i].T * self.cell_volume[:],
                                                   self.A[k * self.num_points:(k + 1) * self.num_points, j]) / self.num_snapshots
                if (j < i):
                    self.C[i, j] = self.C[j, i]
        else:
            self.C[:, :] = np.dot(self.A[:, 0:self.num_snapshots].T, self.A[:, 0:self.num_snapshots]) / self.num_snapshots

    def eigenvalue_problem(self):
        # allocate arrays
        self.energy = np.array(np.zeros((self.num_snapshots), dtype=np.complex64))
        self.temporal_modes = np.array(np.zeros((self.num_snapshots, self.num_snapshots), dtype=np.complex64))

        # determine eigenvalues
        self.energy, self.temporal_modes = linalg.eig(self.C)

    def sort_eigenvalues(self):
        # allocate arrays
        energy_sorted = np.array(np.zeros((self.num_snapshots), dtype=np.float64))
        mode_index = np.array(np.zeros((self.num_snapshots), dtype=np.int))

        # replace invalid modes with zero and small value
        for k in range(0, self.num_snapshots):
            mode_index[k] = k
            if ((math.isnan(self.energy[k].real)) or (math.isnan(self.energy[k].imag))):
                energy_sorted[k] = -1e10
                self.temporal_modes[:, k] = 0
            else:
                energy_sorted[k] = self.energy[k].real

        # sort energy according to size from high to low
        energy_sorted[0:self.num_snapshots], mode_index[0:self.num_snapshots] = zip(*sorted(zip(energy_sorted[:], mode_index[:]), reverse=True))

        # plot histogramms of energy
        nplt.bar('hist_energy', np.linspace(0, self.num_snapshots - 1, self.num_snapshots), self.energy,
                 'Mode Number', '\lambda', self.outputfolder + 'hist_energy', 'original')
        nplt.close('hist_energy')

        nplt.bar('hist_energy_sorted', np.linspace(0, self.num_snapshots - 1, self.num_snapshots), energy_sorted,
                 'Mode Number', '\lambda', self.outputfolder + 'hist_energy_sorted', 'sorted')
        nplt.close('hist_energy_sorted')

        nplt.bar('hist_energy_comp', np.linspace(0, self.num_snapshots - 1, self.num_snapshots), self.energy,
                 'Mode Number', '\lambda', self.outputfolder + 'hist_energy_comp', 'original')
        nplt.bar('hist_energy_comp', np.linspace(0, self.num_snapshots - 1, self.num_snapshots), energy_sorted,
                 'Mode Number', '\lambda', self.outputfolder + 'hist_energy_comp', 'sorted')
        nplt.close('hist_energy_comp')

        nplt.eigs('POD_mode_energies',self.energy,'\lambda',self.outputfolder+'POD_mode_energies', 'original')
        nplt.eigs('POD_mode_energies',energy_sorted,'\lambda',self.outputfolder+'POD_mode_energies', 'sorted')
        nplt.close('POD_mode_energies')

        self.energy[0:self.num_snapshots] = energy_sorted[0:self.num_snapshots]

        # sort temporal modes according to mode index (energy)
        temporal_modes0 = np.array(np.zeros((self.num_snapshots, self.num_snapshots), dtype=np.float64))
        temporal_modes0[0:self.num_snapshots, 0:self.num_snapshots] = self.temporal_modes[0:self.num_snapshots, 0:self.num_snapshots].real.copy()
        for k in range(0, self.num_snapshots):
            self.temporal_modes[0:self.num_snapshots, k] = temporal_modes0[0:self.num_snapshots, mode_index[k]]

    def valid_temporal_modes(self):
        num_valid_modes = 0
        while ((self.energy[num_valid_modes].real / self.energy[0].real > pow(self.tol_CN, 2)) and (num_valid_modes < self.num_snapshots - 2) \
               and (self.energy[num_valid_modes].real > 0)):
            num_valid_modes += 1
            if ((self.energy[num_valid_modes].real / self.energy[0].real > pow(self.tol_CN, 2)) and (
                    self.energy[num_valid_modes].real > 0)):
                num_valid_modes += 1
        print('      Number of valid POD modes with positive energies = ', num_valid_modes)
        self.num_valid_modes = num_valid_modes
        if ((self.num_modes_trunc < 0) or (self.num_modes_trunc > num_valid_modes)):
            self.num_modes_trunc = num_valid_modes

        # calculate resolved energy in %
        total_energy = np.sum(self.energy.real)
        trunc_energy = np.sum(self.energy[0:self.num_modes_trunc])
        resolved_energy = trunc_energy / total_energy * 100

        print('      Resolved energy: = ', resolved_energy, '%')

        # plot histogram with truncated modes
        nplt.bar('hist_energy_comp_trunc', np.linspace(0, self.num_modes_trunc - 1, self.num_modes_trunc),
                 self.energy[0:self.num_modes_trunc],
                 'Mode Number', '\lambda', self.outputfolder + 'hist_energy_comp_trunc',
                 'truncated (' + str("{:1.8f}".format(resolved_energy)) + '% resolved)')
        nplt.bar('hist_energy_comp_trunc', np.linspace(0, self.num_snapshots - 1, self.num_snapshots), self.energy,
                 'Mode Number', '\lambda', self.outputfolder + 'hist_energy_comp_trunc', 'original')
        nplt.close('hist_energy_comp_trunc')

    def scale_temp_modes(self):
        print('\n   Scaling temporal modes ...')
        for j in range(0, self.num_valid_modes):
            temporal_mode_mag = sum(self.temporal_modes[:, j].real * self.temporal_modes[:, j].real) / self.num_snapshots
            self.temporal_modes[:, j] = self.temporal_modes[:, j] * np.sqrt(self.energy[j].real / temporal_mode_mag)
        initial_conditions = np.array(np.zeros((self.num_valid_modes), dtype=np.float64))
        initial_conditions = self.temporal_modes[0, :].real

    def spatial_modes(self):
        print('\n   Calculating truncated spatial modes ...')
        energy_trunc_inv = np.array(np.zeros((self.num_modes_trunc, self.num_modes_trunc), dtype=np.float64))
        energy_trunc_inv = np.diag(np.ones(self.num_modes_trunc) / self.energy[0:self.num_modes_trunc].real, 0)

        self.spatial_modes_trunc = np.array(np.zeros((self.num_components * self.num_points, self.num_modes_trunc), dtype=np.float64))
        self.spatial_modes_trunc = np.dot(np.dot(self.A[:, 0:self.num_snapshots], self.temporal_modes[:, 0:self.num_modes_trunc].real), energy_trunc_inv) / self.num_snapshots

    def pod_output(self):

        # write eigenvalue statistics
        self.write_eigenvalues()

    def write_eigenvalues(self):
        filename = self.outputfolder + 'POD.eigenvalues.dat'

        # calculation of comulative energy
        cumulative_energy = np.array(np.zeros((self.num_valid_modes), dtype=np.float64))
        cumulative_energy[0] = self.energy[0].real

        for i in range(1, self.num_valid_modes):
            cumulative_energy[i] = cumulative_energy[i - 1] + self.energy[i].real

        # total energy
        total_energy = cumulative_energy[self.num_valid_modes - 1]

        print('\n   Writing energies to')
        print('      ', filename)

        # write energy and statistics into file
        file = open(filename, 'w')
        file.write('#\n')
        file.write(
            '# mode, energy, cumulative, percenterage energy, percentage cumulative, condition number (absolute value if negative)\n')
        file.write('#		Note: cummulative energies are set to zero after first negative energy')
        file.write('#\n')
        for i in range(0, self.num_valid_modes):
            file.write('%4.1d %18.10e %18.10e %18.10e %18.10e %18.10e\n' % (
            i + 1, self.energy[i].real, cumulative_energy[i], self.energy[i].real / total_energy * 100,
            cumulative_energy[i] / total_energy * 100, math.sqrt(self.energy[i].real / self.energy[0].real)))
        for i in range(self.num_valid_modes, self.num_snapshots):
            file.write('%4.1d %18.10e %18.10e %18.10e %18.10e %18.10e\n' % (
            i + 1, self.energy[i].real, 0, self.energy[i].real / total_energy * 100, 0,
            math.sqrt(abs(self.energy[i].real / self.energy[0].real))))
        file.close()

class FS:
    def __init__(self, digfilter, pod, outputfolder, energy_target, plot_temp_modes):
        self.outputfolder = outputfolder+'/PODFS/'
        self.num_snapshots = pod.num_snapshots
        self.dt = pod.dt
        self.nm = pod.num_modes_trunc
        self.energy_target = energy_target
        self.temporal_modes = pod.temporal_modes
        self.coordinates = digfilter.coordinates
        self.num_fcs = self.num_snapshots
        self.plot_temp_modes = plot_temp_modes

        print('\n==============================================================')
        print('Calculating Fourier series ...')

        # calculate time constants
        self.time_constants()

        # calculate fourier coefficients for every mode
        self.fourier()

    def time_constants(self):
        stride = 1
        first = 1
        last = self.num_snapshots + first
        self.time = np.linspace(0, (self.num_snapshots - 1) * self.dt * stride, self.num_snapshots)
        self.period = self.time[-1] + (self.time[1] - self.time[0])

        print('period = ', self.period, ' sec')
        print('dt check:', (self.time[1] - self.time[0]), '=', self.dt, '?')

    def fourier(self):
        # allocate arrays
        y2 = np.array(np.zeros((self.num_snapshots, self.nm), dtype=np.float64))
        self.c = np.array(np.zeros((self.num_fcs, self.nm), dtype=np.complex64))
        c_count = np.array(np.zeros((self.nm), dtype=np.int64))
        self.c_ind = np.zeros((self.nm, self.num_fcs), dtype=np.int32)
        self.num_fcs_needed = np.array(np.zeros((self.nm), dtype=np.int64))

        ifig = 0
        # for each mode calculate fourier coefficients
        for i in range(0, self.nm):
            y = self.temporal_modes[:, i]
            for n in range(0, self.num_fcs):
                k = n - self.num_fcs / 2
                ctemp = y * np.exp(-1j *2 * k * np.pi * self.time / self.period)
                self.c[n, i] = ctemp.sum() / ctemp.size

            # rank fourier coefficients in order of coefficient value
            cmod = np.abs(self.c[:, i])
            c_copy = np.array(np.zeros((self.num_fcs), dtype=np.complex64))
            c_copy[:] = self.c[:, i]
            for j in range(0, self.num_fcs):
                self.c_ind[i, j] = j
            cmod, self.c_ind[i, :] = zip(*sorted(zip(cmod, self.c_ind[i, :]), reverse=True))

            # decide on how many should be used based on energy critierion
            energy = 0
            energy_sum = np.sum(np.abs(self.c[:, i]))
            c_count[i] = 0
            while energy < energy_sum * self.energy_target:
                energy += np.abs(self.c[self.c_ind[i, c_count[i]], i])
                c_count[i] += 1

            self.num_fcs_needed[i] = c_count[i]
            # for k in range(0, self.num_fcs):
            #     self.c[k,i] = c_copy[c_ind[i,k]]

            print('The number of fourier coeffients used for mode ', i + 1, ' is ', c_count[i])

            # generate approximate solutions
            if self.plot_temp_modes:
                j = 0
                for x in self.time:
                    f = 0
                    for n in range(0, self.num_fcs_needed[i]):
                        k = self.c_ind[i, n] - self.num_fcs / 2
                        f += self.c[self.c_ind[i, n], i] * np.exp(1j * 2 * k * np.pi * x / self.period)
                    # if j == 0:
                    #	print c[c_ind[i,n],i],k
                    y2[j, i] = f.real
                    j += 1

                # plot for comparison
                ifig += 1
                # print y
                # print time
                # print y2[:,i]
                # print ifig
                # print i
                # nplt.timeseries(ifig, y, self.time, '\,', self.outputfolder + 'POD_tmode' + str(i))
                # nplt.timeseries(ifig, y2[:, i], self.time, '\,', self.outputfolder + 'POD_tmode_recon' + str(i))
                nplt.timeseries_comp(ifig, y, y2[:, i], self.time, '\,', self.outputfolder + 'POD_tmode_comp' + str(i))


class HDF5:
    def __init__(self, pod, fourier):
        self.outputfolder = pod.outputfolder
        self.num_modes = pod.num_modes_trunc
        self.period = fourier.period
        self.num_FC = fourier.num_fcs_needed
        self.FC = fourier.c
        self.num_points = pod.num_points
        self.mean_field = pod.mean_field
        self.num_vars = pod.num_components
        self.spatial_modes = pod.spatial_modes_trunc
        self.fluct_comps = pod.fluct_comps
        self.coordinates = np.transpose(pod.coordinates)
        self.coord_comps = pod.coord_comps
        self.num_coord_comps = pod.num_coord_comps
        self.c_ind = fourier.c_ind



        self.write_HDF5()

    def write_HDF5(self):


        f = h5py.File(self.outputfolder+'PODFS_HYDRA.hdf5','w')



        # create 'main' group
        main = f.create_group("main")

        # write number of modes as attribute
        main.attrs.create("number_POD_modes", self.num_modes)
        # write period as attribute
        main.attrs.create("period", self.period)
        main.attrs.create("coordinates", self.coord_comps)
        main.attrs.create("number_points", self.num_points)
        main.attrs.create("number_variables",  self.num_vars)
        main.attrs.create("variables", self.fluct_comps)
        main.attrs.create("number_coord_components", self.num_coord_comps)

        # create dataset of number of fourier coefficients per mode
        N_FC = main.create_dataset("number_FC_per_mode", data = self.num_FC)

        # create dataset of ranking of fourier coefficients per mode
        FC_ind = main.create_dataset("Rank_FC", data=self.c_ind)

        # create dataset of fourier coefficients
        FC = main.create_dataset("FC", data = self.FC)

        # create dataset of coordinates
        coord = main.create_dataset("coordinates", data = self.coordinates)

        # create dataset of mean field
        data = main.create_dataset("mean", data = self.mean_field)

        # create "modes" group
        modes = main.create_group("modes")

        #for each mode
        for i in range (0,self.num_modes):

            counter = '%4.4i'% (i+1)
            # create dataset
            data = modes.create_dataset("mode_"+counter,data = self.spatial_modes[:,i])

        f.close()