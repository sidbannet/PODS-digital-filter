# module for HDF5 saving of PODFS modes

# This module saves the PODFS output as a single HDF5 file

# Format is:


import h5py
import numpy as np

def write_HDF5(i_d):

	f = h5py.File('PODFS/PODFS.hdf5','w')

	# create 'main' group
	main = f.create_group("main") 

	# write number of modes as attribute
	main.attrs["N_POD"] = i_d.nm
	# write period as attribute
	main.attrs["period"] = i_d.period
	
	# create dataset of number of fourier coefficients per mode
	N_FC = main.create_dataset("N_FC",(i_d.nm,), dtype='i')
	N_FC[:] = i_d.N_FC

	# create dataset of fourier coefficients
	FC = main.create_dataset("FC",(np.sum(i_d.N_FC)*3,), dtype=np.float64)
	FC[:] = i_d.FC.reshape(np.sum(i_d.N_FC)*3,order='F')

	# create dataset of mean field
	data = main.create_dataset("mean",(i_d.num_points*6,), dtype=np.float64)
	data[:] = i_d.mean.reshape(i_d.num_points*6,order='F')

	# add number of points as attribute
        data.attrs["Np"] = i_d.num_points
        # add number of variables as  attribute
        data.attrs["Nvar"] = 6
        # add variables as attribute (needs dummy variable at the end for precise to read!)
        data.attrs["Vars"] = np.string_('x,y,z,u,v,w,dummy')
	# add scaling factors as an attribute
        data.attrs["SF"] = [1.,1.,1.,1.,1.,1.]

	# create "modes" group 
	modes = main.create_group("modes") 
	
	#for each mode
	for i in range (0,i_d.nm):
		
		counter = '%4.4i'% (i+1)
		# create dataset
		data = modes.create_dataset("mode_"+counter,(i_d.num_points*6,), dtype=np.float64)
		data[:] = i_d.modes[i,:,:].reshape(i_d.num_points*6,order='F')

		# add number of points as attribute
        	data.attrs["Np"] = i_d.num_points
        	# add number of variables as  attribute
        	data.attrs["Nvar"] = 6
        	# add variables as attribute (needs dummy variable at end)
        	data.attrs["Vars"] = np.string_('x,y,z,u,v,w,dummy')
		# add scaling factors as an attribute
        	data.attrs["SF"] = [1.,1.,1.,1.,1.,1.]

	f.close()
