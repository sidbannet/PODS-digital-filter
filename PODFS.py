# module PODFS.py 

##@package digitalfilters.py
#
# This module is for computing the PODFS of LES inlet data generated
# using the digital filter method. 
#
#DEPENDANCIES:
#vtk,
#numpy,
##
#Version 1.0 
#07/08/2018
#Nicholas Treleaven
#n.c.w.treleaven@lboro.ac.uk

#=======================================================================
# import libraries
import vtk
import vtk.util.numpy_support as VN
import numpy as np
import sys
import nplotlib as plt 
import nsigproclib_no_mpi as sp
import os
import string
import math
from numpy import linalg
#import io_vtk_m as old
#from progressbar import * need to remove this from spray stats module

#=======================================================================
def get_1_var(filename,var_name):
        #reader = vtk.vtkUnstructuredGridReader()
        #reader.ReadAllScalarsOn()
        #reader.ReadAllVectorsOn()
        #reader.SetFileName(filename)
        #reader.Update()
        # extract feature removed, it takes too long
        #if ( (x_min<x_max) and (y_min<y_max) and (z_min<z_max) ):
        #       extract = vtk.vtkExtractUnstructuredGrid()
        #       extract.SetExtent(x_min, x_max, y_min, y_max, z_min, z_max)
        #       extract.SetInputData(reader.GetOutput())
        ##       extract.MergingOn()
        #       extract.Update()
        #       cell_volume = VN.vtk_to_numpy(extract.GetOutput().GetPointData().GetScalars(var_name))
        #else:
        cell_volume = VN.vtk_to_numpy(filename.GetPointData().GetScalars(var_name))
	
	return cell_volume

#======================================================================
def get_1_vector(filename,var_name):
        
        cell_volume = VN.vtk_to_numpy(filename.GetPointData().GetVectors(var_name))
	
        return cell_volume

#=======================================================================
def set_1_var(grid,var_name,var):
        u_vtk = VN.numpy_to_vtk(var)
        u_vtk.SetName(var_name)
        grid.GetPointData().AddArray(u_vtk)

#=======================================================================
#=======================================================================
def get_grid_template(filename,i_d):
        grid = vtk.vtkUnstructuredGrid()
	reader = vtk.vtkUnstructuredGridReader()
	reader.ReadAllScalarsOff()
	reader.ReadAllVectorsOff()
	reader.SetFileName(filename)
	reader.Update()
        # extract feature removed, it takes too long
	if ( (i_d.x_min<i_d.x_max) and (i_d.y_min<i_d.y_max) and (i_d.z_min<i_d.z_max) ):
		extract = vtk.vtkExtractUnstructuredGrid()
		extract.SetExtent(i_d.x_min, i_d.x_max, i_d.y_min, i_d.y_max, i_d.z_min, i_d.z_max)
		extract.SetInputData(reader.GetOutput())
	#	extract.MergingOn()
		extract.Update()
		grid.DeepCopy(extract.GetOutput())
	else:
                grid.DeepCopy(reader.GetOutput())

        return grid

#=======================================================================
def get_grid(filename,i_d):
        grid = vtk.vtkUnstructuredGrid()
        reader = vtk.vtkUnstructuredGridReader()
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.SetFileName(filename)
        reader.Update()
        # extract feature removed, it takes too long
        if ( (i_d.x_min<i_d.x_max) and (i_d.y_min<i_d.y_max) and (i_d.z_min<i_d.z_max) ):
                extract = vtk.vtkExtractUnstructuredGrid()
                extract.SetExtent(i_d.x_min, i_d.x_max, i_d.y_min, i_d.y_max, i_d.z_min, i_d.z_max)
                extract.SetInputData(reader.GetOutput())
        #       extract.MergingOn()
                extract.Update()
                grid.DeepCopy(extract.GetOutput())
        else:
                grid.DeepCopy(reader.GetOutput())

        return grid

#=======================================================================
def get_mean_template(filename,x_min, x_max, y_min, y_max, z_min, z_max):
	reader = vtk.vtkUnstructuredGridReader()
	reader.ReadAllScalarsOff()
	reader.ReadAllVectorsOff()
	reader.SetFileName(filename)
	reader.Update()
        # extract feature removed, it takes too long
	if ( (x_min<x_max) and (y_min<y_max) and (z_min<z_max) ):
		extract = vtk.vtkExtractUnstructuredGrid()
		extract.SetExtent(x_min, x_max, y_min, y_max, z_min, z_max)
		extract.SetInputData(reader.GetOutput())
	#	extract.MergingOn()
		extract.Update()
		grid = extract
	else:
                grid = reader

        return grid

#=======================================================================
def write_field(num_points,num_components,grid,field,filename):
	u_temp = np.array(np.zeros((num_points,3), dtype=np.float64))
	for k in range(0,num_components):
		u_temp[:,k] = field[k*num_points:(k+1)*num_points]
	u_vtk = VN.numpy_to_vtk(u_temp)
	u_vtk.SetName('u')
	grid.GetPointData().AddArray(u_vtk)
	writer = vtk.vtkUnstructuredGridWriter()
	writer.SetFileTypeToBinary()
	writer.SetInputData(grid)
	writer.SetFileName(filename)
	writer.Write()

#=======================================================================
def write_field_cell(num_points,num_components,grid,field,filename):
	u_temp = np.array(np.zeros((num_points,3), dtype=np.float64))
	for k in range(0,num_components):
		u_temp[:,k] = field[k*num_points:(k+1)*num_points]
	u_vtk = VN.numpy_to_vtk(u_temp)
	u_vtk.SetName('u')
	grid.GetCellData().AddArray(u_vtk)
	writer = vtk.vtkUnstructuredGridWriter()
	writer.SetFileTypeToBinary()
	writer.SetInputData(grid)
	writer.SetFileName(filename)
	writer.Write()

#=======================================================================
def write_stats(num_points,num_components,grid,vectors,vector_names, scalars, scalar_names, s_cell, s_cell_names, filename):
        u_temp = np.array(np.zeros((num_points,3), dtype=np.float64))
        scalar_temp = np.array(np.zeros((num_points), dtype=np.float64))
        #print len(vector_names)
	grid.GetPointData().Initialize()
        for i in range (0,len(vector_names)):
            for k in range(0,num_components):
                u_temp[:,k] = vectors[k*num_points:(k+1)*num_points,i]
            u_vtk = VN.numpy_to_vtk(u_temp,1)
            u_vtk.SetName(vector_names[i])
            grid.GetPointData().AddArray(u_vtk)
	    #print 'added ', vector_names[i] ,'of'
        for i in range (0, len(scalar_names)):
            scalar_temp = scalars[i]
            u_vtk = VN.numpy_to_vtk(scalar_temp)
            u_vtk.SetName(scalar_names[i])
            grid.GetPointData().AddArray(u_vtk)
	for i in range (0,len(s_cell_names)):
            #scalar_temp = s_cell[i]
            u_vtk = VN.numpy_to_vtk(s_cell[i])
            u_vtk.SetName(s_cell_names[i])
            grid.GetCellData().AddArray(u_vtk)

	print 'writing'
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileTypeToBinary()
        writer.SetInputData(grid)
        writer.SetFileName(filename)
        writer.Write()
	print 'written'

#=======================================================================
def write_mean_field(num_points,num_components,grid,mean_field):
	print '   Writing mean field to'
	filename = './results/spatial_meanfield.vtk'
	print '      ', filename
	write_field(num_points,num_components,grid,mean_field,filename)
#=======================================================================
def write_mean_field_vor(num_points,num_components,grid,mean_field):
	print '   Writing mean field to'
	filename = './results/spatial_meanfield_vor.vtk'
	print '      ', filename
	write_field(num_points,num_components,grid,mean_field,filename)

#=======================================================================
def write_spatial_POD_modes(num_modes_to_write,num_components,grid,spatial_modes_trunc,var_name,rdir):
	print '\n   Writing spatial POD modes to'
        grid1 = vtk.vtkUnstructuredGrid()
        grid1.DeepCopy(grid)
        grid1.GetPointData().Initialize()
        grid1.GetCellData().Initialize()
        num_cells = grid1.GetNumberOfCells()
	for j in range(0,num_modes_to_write):
		j_index = j+1
		filename = rdir+'POD.spatial_mode_'+var_name+'_' + '%04d'%j_index + '.vtk'
		print '      ', filename
		if len(var_name.split(',')) > 1:
        		var_name_list = var_name.split(',')
        		cc=0
        		for i in range (0,len(var_name_list)):
        			
				if (var_name_list[i] == 'velocity' or var_name_list[i] == 'U' or var_name_list[i] == 'SprayVelocity'):
                        		u_temp = np.array(np.zeros((num_cells,3), dtype=np.float64))
                        		u_temp2 = np.array(np.zeros((num_cells), dtype=np.float64))
                        		u_mag = np.array(np.zeros((num_cells), dtype=np.float64))
                        		for k in range(0,3):
                                		u_temp[:,k] = spatial_modes_trunc[cc*num_cells:(cc+1)*num_cells,j]
                                		u_temp2[:] = spatial_modes_trunc[cc*num_cells:(cc+1)*num_cells,j]
                                                print k, cc
                                                cc+=1
                                		u_mag[:] = u_mag[:]+u_temp2[:]*u_temp2[:]
                                		u_vtk = VN.numpy_to_vtk(u_temp2,1)
                                		u_vtk.SetName(var_name_list[i]+'_'+str(k+1)+'_POD')
                                		grid1.GetCellData().AddArray(u_vtk)

                        		u_mag[:] = np.sqrt(u_mag[:])
                        		u_vtk = VN.numpy_to_vtk(u_mag,1)
                        		u_vtk.SetName(var_name_list[i]+'_magnitude_POD')
                        		grid1.GetCellData().AddArray(u_vtk)
                        

                        		u_vtk = VN.numpy_to_vtk(u_temp,1)
                        		u_vtk.SetName(var_name_list[i]+'_POD')
                        		grid1.GetCellData().AddArray(u_vtk)
                
                		else:
                        		u_temp = np.array(np.zeros((num_cells), dtype=np.float64))
                        		u_temp[:] = spatial_modes_trunc[cc*num_cells:(cc+1)*num_cells,j]
                        		cc+=1
                        		u_vtk = VN.numpy_to_vtk(u_temp,1)
                        		u_vtk.SetName(var_name_list[i]+'_POD')
                        		grid1.GetCellData().AddArray(u_vtk)
                   
                elif (var_name == 'velcity' or 'U'):
                        u_temp = np.array(np.zeros((num_cells,num_components), dtype=np.float64))
                        u_temp2 = np.array(np.zeros((num_cells), dtype=np.float64))
                        u_mag = np.array(np.zeros((num_cells), dtype=np.float64))
                        for k in range(0,num_components):
                                u_temp[:,k] = spatial_modes_trunc[k*num_cells:(k+1)*num_cells,j]
                                u_temp2[:] = spatial_modes_trunc[k*num_cells:(k+1)*num_cells,j]
                                u_mag[:] = u_mag[:]+u_temp2[:]*u_temp2[:]
                                u_vtk = VN.numpy_to_vtk(u_temp2,1)
                                u_vtk.SetName(var_name+'_'+str(k+1)+'_POD')
                                grid1.GetCellData().AddArray(u_vtk)

                        u_mag[:] = np.sqrt(u_mag[:])
                        u_vtk = VN.numpy_to_vtk(u_mag,1)
                        u_vtk.SetName(var_name+'_magnitude_POD')
                        grid1.GetCellData().AddArray(u_vtk)
                        

                        u_vtk = VN.numpy_to_vtk(u_temp,1)
                        u_vtk.SetName(var_name+'_POD')
                        grid1.GetCellData().AddArray(u_vtk)
                
                else:
                        u_temp = np.array(np.zeros((num_cells), dtype=np.float64))
                        u_temp[:] = spatial_modes_trunc[:,j]
                        u_vtk = VN.numpy_to_vtk(u_temp,1)
                        u_vtk.SetName(var_name+'_POD')
                        grid1.GetCellData().AddArray(u_vtk)
                        
                c2p = vtk.vtkCellDataToPointData()
                c2p.SetInputData(grid1)
                c2p.Update()
                writer = vtk.vtkUnstructuredGridWriter()
                writer.SetFileTypeToBinary()
                writer.SetInputData(c2p.GetOutput())
                writer.SetFileName(filename)
                writer.Write()
#=======================================================================
def write_spatial_POD_modes_i_d(num_modes_to_write,num_components,grid,spatial_modes_trunc,var_name,rdir,i_d):
	print '\n   Writing spatial POD modes to'
        grid1 = vtk.vtkUnstructuredGrid()
        grid1.DeepCopy(grid)
        grid1.GetPointData().Initialize()
        grid1.GetCellData().Initialize()
        num_cells = grid1.GetNumberOfCells()
	for j in range(0,num_modes_to_write):
		j_index = j+1
		filename = rdir+'POD.spatial_mode_'+var_name+'_' + '%04d'%j_index + '.vtk'
		print '      ', filename
		if len(var_name.split(',')) > 1:
        		var_name_list = var_name.split(',')
        		cc=0
        		for i in range (0,len(var_name_list)):
                            #print var_name_list[i],i_d.is_POD_var_vec[i]
                            if not i_d.is_POD_var_vec:
				if (var_name_list[i] == 'velocity' or var_name_list[i] == 'U' or var_name_list[i] == 'SprayVelocity'):
                        		u_temp = np.array(np.zeros((num_cells,3), dtype=np.float64))
                        		u_temp2 = np.array(np.zeros((num_cells), dtype=np.float64))
                        		u_mag = np.array(np.zeros((num_cells), dtype=np.float64))
                        		for k in range(0,3):
                                		u_temp[:,k] = spatial_modes_trunc[cc*num_cells:(cc+1)*num_cells,j]
                                		u_temp2[:] = spatial_modes_trunc[cc*num_cells:(cc+1)*num_cells,j]
                                                print k, cc
                                                cc+=1
                                		u_mag[:] = u_mag[:]+u_temp2[:]*u_temp2[:]
                                		u_vtk = VN.numpy_to_vtk(u_temp2,1)
                                		u_vtk.SetName(var_name_list[i]+'_'+str(k+1)+'_POD')
                                		grid1.GetCellData().AddArray(u_vtk)

                        		u_mag[:] = np.sqrt(u_mag[:])
                        		u_vtk = VN.numpy_to_vtk(u_mag,1)
                        		u_vtk.SetName(var_name_list[i]+'_magnitude_POD')
                        		grid1.GetCellData().AddArray(u_vtk)
                        

                        		u_vtk = VN.numpy_to_vtk(u_temp,1)
                        		u_vtk.SetName(var_name_list[i]+'_POD')
                        		grid1.GetCellData().AddArray(u_vtk)
                
                		else:
                        		u_temp = np.array(np.zeros((num_cells), dtype=np.float64))
                        		u_temp[:] = spatial_modes_trunc[cc*num_cells:(cc+1)*num_cells,j]
                        		cc+=1
                        		u_vtk = VN.numpy_to_vtk(u_temp,1)
                        		u_vtk.SetName(var_name_list[i]+'_POD')
                        		grid1.GetCellData().AddArray(u_vtk)
                            else:
                                if (i_d.PODVarVec[i] == 'true'):
                        		u_temp = np.array(np.zeros((num_cells,3), dtype=np.float64))
                        		u_temp2 = np.array(np.zeros((num_cells), dtype=np.float64))
                        		u_mag = np.array(np.zeros((num_cells), dtype=np.float64))
                        		for k in range(0,3):
                                		u_temp[:,k] = spatial_modes_trunc[cc*num_cells:(cc+1)*num_cells,j]
                                		u_temp2[:] = spatial_modes_trunc[cc*num_cells:(cc+1)*num_cells,j]
                                                print k, cc
                                                cc+=1
                                		u_mag[:] = u_mag[:]+u_temp2[:]*u_temp2[:]
                                		u_vtk = VN.numpy_to_vtk(u_temp2,1)
                                		u_vtk.SetName(var_name_list[i]+'_'+str(k+1)+'_POD')
                                		grid1.GetCellData().AddArray(u_vtk)

                        		u_mag[:] = np.sqrt(u_mag[:])
                        		u_vtk = VN.numpy_to_vtk(u_mag,1)
                        		u_vtk.SetName(var_name_list[i]+'_magnitude_POD')
                        		grid1.GetCellData().AddArray(u_vtk)
                        

                        		u_vtk = VN.numpy_to_vtk(u_temp,1)
                        		u_vtk.SetName(var_name_list[i]+'_POD')
                        		grid1.GetCellData().AddArray(u_vtk)
                
                		else:
                        		u_temp = np.array(np.zeros((num_cells), dtype=np.float64))
                        		u_temp[:] = spatial_modes_trunc[cc*num_cells:(cc+1)*num_cells,j]
                        		cc+=1
                        		u_vtk = VN.numpy_to_vtk(u_temp,1)
                        		u_vtk.SetName(var_name_list[i]+'_POD')
                        		grid1.GetCellData().AddArray(u_vtk)
                else:
                 if not i_d.is_POD_var_vec:
                  if (var_name == 'velocity' or 'U'):
                        u_temp = np.array(np.zeros((num_cells,num_components), dtype=np.float64))
                        u_temp2 = np.array(np.zeros((num_cells), dtype=np.float64))
                        u_mag = np.array(np.zeros((num_cells), dtype=np.float64))
                        for k in range(0,num_components):
                                u_temp[:,k] = spatial_modes_trunc[k*num_cells:(k+1)*num_cells,j]
                                u_temp2[:] = spatial_modes_trunc[k*num_cells:(k+1)*num_cells,j]
                                u_mag[:] = u_mag[:]+u_temp2[:]*u_temp2[:]
                                u_vtk = VN.numpy_to_vtk(u_temp2,1)
                                u_vtk.SetName(var_name+'_'+str(k+1)+'_POD')
                                grid1.GetCellData().AddArray(u_vtk)

                        u_mag[:] = np.sqrt(u_mag[:])
                        u_vtk = VN.numpy_to_vtk(u_mag,1)
                        u_vtk.SetName(var_name+'_magnitude_POD')
                        grid1.GetCellData().AddArray(u_vtk)
                        

                        u_vtk = VN.numpy_to_vtk(u_temp,1)
                        u_vtk.SetName(var_name+'_POD')
                        grid1.GetCellData().AddArray(u_vtk)
                
                  else:
                        u_temp = np.array(np.zeros((num_cells), dtype=np.float64))
                        u_temp[:] = spatial_modes_trunc[:,j]
                        u_vtk = VN.numpy_to_vtk(u_temp,1)
                        u_vtk.SetName(var_name+'_POD')
                        grid1.GetCellData().AddArray(u_vtk)
                 else:       
                  if (i_d.PODVarVec == 'true'):
                        u_temp = np.array(np.zeros((num_cells,num_components), dtype=np.float64))
                        u_temp2 = np.array(np.zeros((num_cells), dtype=np.float64))
                        u_mag = np.array(np.zeros((num_cells), dtype=np.float64))
                        for k in range(0,num_components):
                                u_temp[:,k] = spatial_modes_trunc[k*num_cells:(k+1)*num_cells,j]
                                u_temp2[:] = spatial_modes_trunc[k*num_cells:(k+1)*num_cells,j]
                                u_mag[:] = u_mag[:]+u_temp2[:]*u_temp2[:]
                                u_vtk = VN.numpy_to_vtk(u_temp2,1)
                                u_vtk.SetName(var_name+'_'+str(k+1)+'_POD')
                                grid1.GetCellData().AddArray(u_vtk)

                        u_mag[:] = np.sqrt(u_mag[:])
                        u_vtk = VN.numpy_to_vtk(u_mag,1)
                        u_vtk.SetName(var_name+'_magnitude_POD')
                        grid1.GetCellData().AddArray(u_vtk)
                        

                        u_vtk = VN.numpy_to_vtk(u_temp,1)
                        u_vtk.SetName(var_name+'_POD')
                        grid1.GetCellData().AddArray(u_vtk)
                
                  else:
                        u_temp = np.array(np.zeros((num_cells), dtype=np.float64))
                        u_temp[:] = spatial_modes_trunc[:,j]
                        u_vtk = VN.numpy_to_vtk(u_temp,1)
                        u_vtk.SetName(var_name+'_POD')
                        grid1.GetCellData().AddArray(u_vtk)

          #      c2p = vtk.vtkCellDataToPointData()
           #     c2p.SetInputData(grid1)
            #    c2p.Update()
                writer = vtk.vtkUnstructuredGridWriter()
                writer.SetFileTypeToBinary()
                writer.SetInputData(grid1)
                writer.SetFileName(filename)
                writer.Write()

#=======================================================================
def write_mean_field2(grid,spatial_modes_trunc,var_name,num_components,rdir):
	print '\n   Writing spatial mean field to'
        grid1 = vtk.vtkUnstructuredGrid()
        grid1.DeepCopy(grid)
        grid1.GetPointData().Initialize()
        grid1.GetCellData().Initialize()
        num_cells = grid1.GetNumberOfCells()


        filename = rdir+'POD.spatial_mean_field_'+var_name+'.vtk'
        print '      ', filename
        if len(var_name.split(',')) > 1:
        	var_name_list = var_name.split(',')
        	cc=0
        	for i in range (0,len(var_name_list)):
        		temp = var_name_list[i] 
        		if (temp == 'velocity' or temp == 'U'  or temp == 'SprayVelocity' ):
                		u_temp = np.array(np.zeros((num_cells,3), dtype=np.float64))
                		for k in range(0,3):
					print var_name_list[i], cc,k
                        		u_temp[:,k] = spatial_modes_trunc[cc*num_cells:(cc+1)*num_cells]
                        		cc+=1
                                        #                               print k
                                u_vtk = VN.numpy_to_vtk(u_temp,1)
                        	u_vtk.SetName(var_name_list[i]+'_POD')
                        	grid1.GetCellData().AddArray(u_vtk)
        		else:
                		u_temp = np.array(np.zeros((num_cells), dtype=np.float64))
                		u_temp[:] = spatial_modes_trunc[cc*num_cells:(cc+1)*num_cells]
                		cc+=1
                		u_vtk = VN.numpy_to_vtk(u_temp,1)
                		u_vtk.SetName(var_name_list[i]+'_POD')
                		grid1.GetCellData().AddArray(u_vtk)
                        
        elif (var_name == 'velocity' or 'U'):
                u_temp = np.array(np.zeros((num_cells,num_components), dtype=np.float64))
                for k in range(0,num_components):
                        u_temp[:,k] = spatial_modes_trunc[k*num_cells:(k+1)*num_cells]
                        u_vtk = VN.numpy_to_vtk(u_temp,1)
                        u_vtk.SetName(var_name+'_POD')
                        grid1.GetCellData().AddArray(u_vtk)
        else:
                u_temp = np.array(np.zeros((num_cells), dtype=np.float64))
                u_temp[:] = spatial_modes_trunc
                u_vtk = VN.numpy_to_vtk(u_temp)
                u_vtk.SetName(var_name+'_POD')
                grid1.GetCellData().AddArray(u_vtk)

        #print 'ok'
                        
        #c2p = vtk.vtkCellDataToPointData()
        #c2p.SetInputData(grid1)
        #c2p.Update()

        #print 'check'

        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileTypeToBinary()
        writer.SetInputData(grid1)
        writer.SetFileName(filename)
        writer.Write()

def find_axis_points(x,y,num_points,nxtop,nxbottom,nxleft,nxright):
	index1=np.where(x==np.amax(x))
	nxtop[:]=np.product(np.shape(index1))
	
	index1=np.where(x==np.amin(x))
	nxbottom[:]=np.product(np.shape(index1))
	
	index1=np.where(y==np.amin(y))
	nxleft[:]=np.product(np.shape(index1))
	
	index1=np.where(y==np.amax(y))
	nxright[:]=np.product(np.shape(index1))
	

#=======================================================================
def extract_plane(n1,n2,n3,o1,o2,o3,grid,i_d):

    # The cutter uses an implicit function to perform the cutting.
    # Here we define a plane, specifying its center and normal.
    # Then we assign the plane to the cutter.
    pl3d_output = vtk.vtkUnstructuredGrid()
    pl3d_output.DeepCopy(grid)
    #if ( (x_min<x_max) and (y_min<y_max) and (z_min<z_max) ):
    #    extract1 = vtk.vtkExtractUnstructuredGrid()
    #    extract1.SetExtent(x_min, x_max, y_min, y_max, z_min, z_max)
    #    extract1.SetInputData(pl3d_output)
        #extract1.MergingOn()
        #extract1.CellClippingOn()
    #    extract1.Update()
    #    extract = extract1.GetOutput()
    #else:
    extract = pl3d_output
    plane = vtk.vtkPlane()
    planeS = vtk.vtkPlaneSource()
    planeS.SetResolution(i_d.plane_res,i_d.plane_res)
    bounds = extract.GetBounds()
    center = extract.GetCenter()
    s1 = bounds[1]-bounds[0]
    s2 = bounds[3]-bounds[2]
    s3 = bounds[5]-bounds[4]
    trans = vtk.vtkTransform()
    trans1 = vtk.vtkTransform()
    planeS.SetNormal(n1, n2, n3)
    trans1.Scale(s1,s2,s3)
    plane.SetNormal(n1, n2, n3)
 #   print o1,o2,o3 
    if (o1 == i_d.t_o[0]):
	o1 = center[0]
    if (o2 == i_d.t_o[1]):
	o2 = center[1]
    if (o3 == i_d.t_o[2]):
	o3 = center[2]
  #  print center,o1,o2,o3
  #  print bounds,n1,n2,n3
    plane.SetOrigin(o1,o2,o3)
	#planeS.SetCenter(o1,o2,o3)
    trans.Translate(o1,o2,o3)
 
   
    transf1 = vtk.vtkTransformPolyDataFilter()
    transf1.SetInputConnection(planeS.GetOutputPort())
    transf1.SetTransform(trans1)
    transf = vtk.vtkTransformPolyDataFilter()
    transf.SetInputConnection(transf1.GetOutputPort())
    transf.SetTransform(trans)
    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection(transf.GetOutputPort())
    probe.SetSourceData(extract)
    probe.Update()
    extractS = probe

    planeCut = vtk.vtkCutter()
    #planeCut.GenerateTrianglesOff()
    planeCut.SetInputData(extract)
    planeCut.SetCutFunction(plane)
    cutMapper = vtk.vtkPolyDataMapper()
    cutMapper.SetInputConnection(planeCut.GetOutputPort())
    #cutMapper.SetInputConnection(probe.GetOutputPort())

   

    return cutMapper, extract, extractS

#=======================================================================
def write_vi_dat(filename,num_snapshots,time_vec,pprobe,num_probes,T_amb,p_amb,fs,fmax):

	print '\n Number of probes:', num_probes

	target = open(filename+'.dat','w')

	target.write('Point Probes - '+ filename+'\n')
	target.write('\n')
	target.write('# Date : 18 Nov 2015 \n')	
	target.write('# Time : 15:36 \n')
	target.write('# Ambient Pressure (Pa) = '+str(p_amb)+'\n')
	target.write('# Ambient Temperature (C) = '+str(T_amb-273.16)+' \n')
	target.write('# Rig Mass Flow (kg/s) =  0.0000 \n')
	target.write('# Mach number = 0.0000 \n')
	target.write('# Primary Pressure Drop =    0.00 \n')
	target.write('# No. of Sensors = '+str(num_probes)+'\n')
	target.write('# Sampling Frequency (kHz) = '+str(fs)+'\n')
	target.write('# No. of Samples = '+str(num_snapshots)+'\n')
	target.write('# No. of Time Histories =  1 \n')
	target.write('\n')
	target.write('Variables = t')
	for i in range (0, num_probes):
		target.write(',K'+str(i+1))
	target.write('\n')
	target.write('\n')
	for i in range (0, num_probes):
		target.write('# Peak Frequency '+str(i+1)+' (Hz) = '+str(fmax[i])+'\n')
	target.write('\n')
	target.write('ZONE T="R82KU01, H01", i=32768, j=1, F=POINT \n')
	target.write('\n')

	for i in range (0,num_snapshots):
		target.write("%9.9f"% time_vec[i])
		target.write('\t')
		target.write(('\t'.join(map(str,pprobe[:,i])))+'\n') #for j in range (0,num_probes)))
		#target.write('\n'.join(' '.join(str(pprobe[j,i])) for j in range (0,num_probes))

	target.close()
#=======================================================================
def save_Planes(extract,i_d,time,POD_plane_vel,POD_plane_areas,snap):

        i_d.plane_count = np.amin(np.array([len(i_d.save_planes),len(i_d.POD_planes)]))
        if (i_d.make_polar == 'true'):
                Atemp = np.zeros((i_d.plane_cells),dtype=np.float64)
    
        for k in range (0,len(i_d.plane_names)):
                for j in range (0,i_d.plane_count):
                        if (i_d.save_planes[j] == k or i_d.POD_planes[j] == k):
                        
                                if not os.path.exists(i_d.snapshot_dir+'post/'+i_d.plane_names[k]):
                                        os.mkdir(i_d.snapshot_dir+'post/'+i_d.plane_names[k])

                                n1=i_d.plane_normals[k,0]
                                n2=i_d.plane_normals[k,1]
                                n3=i_d.plane_normals[k,2]
                                o1=i_d.plane_centres[k,0]
                                o2=i_d.plane_centres[k,1]
                                o3=i_d.plane_centres[k,2]
                                cutMapper, pl3d_output, src_output1 = extract_plane(n1,n2,n3,o1,o2,o3,extract,i_d)
                                #plt.vtkcontourf(cutMapper,pl3d_output,src_output,'P',100,snapshot_dir+'/snap'+str(k),plane_normals[k,:],(0,0))
                                
                                src_output = src_output1.GetOutput()

                               
                                if (i_d.save_planes[j] == k): 
                                        #print i_d.te_var
                                        u = VN.vtk_to_numpy(src_output.GetPointData().GetVectors(i_d.vel_var))
                                        p = VN.vtk_to_numpy(src_output.GetPointData().GetScalars(i_d.p_var))
                                        tke = VN.vtk_to_numpy(src_output.GetPointData().GetScalars(i_d.te_var))
                                        epsi = VN.vtk_to_numpy(src_output.GetPointData().GetScalars(i_d.ed_var))
                                                
                                        points = VN.vtk_to_numpy(src_output.GetPoints().GetData())
                                
                                        #print np.shape(u_temp_read),np.shape(points)
                                        timestring = "%.5E" % time
                        
                                        filename = i_d.snapshot_dir+'post/'+i_d.plane_names[k]+'/'+timestring+'.prf'

                                        target = open(filename,'w')

                                        target.write('# ' + i_d.plane_names[k] + ' # name of the profile\n')
                                        target.write('# turbulence model, ' + i_d.turbulence_model+'\n')
                                        plane_rhs = (o1+i_d.t_o[0])*n1+(o2+i_d.t_o[1])*n2+(o3+i_d.t_o[2])*n3
                                        target.write('# plane normal and translation ' + str(n1)+'\t'+str(n2)+'\t'+str(n3)+'\t'+str(plane_rhs)+'\n')
                                        #target.write(str(src_output.GetNumberOfPoints())+'\n')
                                
                                        target.write('type, xyz # type of profile (rad or xyz)\n')
                                        target.write('localcs,origin,0,0,0 # origin of local coordinate system\n')
                                        target.write('localcs,xaxis,1,0,0 # x axis direction of local coordinate system\n')
                                        target.write('localcs,yaxis,0,1,0 # y axis direction of local coordinate system\n')
                                        target.write('localcs,zaxis,0,0,1 # z axis direction of local coordinate system\n')
                                        target.write('tolerance, 1.00E-08 # tolerance\n')
                                        target.write('scale,1,1,1,1,1,1,1,1 # scaling factors\n')
                                        if (i_d.turbulence_model == 'k_epsilon'):
                                                target.write('data,x,y,z,u,v,w,k,e\n')
                                        if (i_d.turbulence_model == 'k_omega'):
                                                target.write('data,x,y,z,u,v,w,k,sdr\n')
                                        for i in range (0,src_output.GetNumberOfPoints()):
                                                target.write(str(points[i,0])+','+str(points[i,1])+','+str(points[i,2])+','+str(u[i,0])+','+str(u[i,1])+','+str(u[i,2]) + ',' + str(tke[i])+','+str(epsi[i])+'\n')
                                            #target.write("%9.9f"% time_vec[i])
                                        #target.write('\t')
                                        #target.write(('\t'.join(map(str,pprobe[:,i])))+'\n') #for j in range (0,num_probes)))
                                        #target.write('\n'.join(' '.join(str(pprobe[j,i])) for j in range (0,num_probes))
                                        
                                        target.close()
                                        
                                if (i_d.POD_planes[j] == k):
                                        POD_plane_areas[j,:],cell_centres = calc_cell_areas(src_output1.GetOutput())
                                        p2c = vtk.vtkPointDataToCellData()
                                        p2c.SetInputConnection(src_output1.GetOutputPort())
                                        p2c.Update()
                                        #print np.shape(POD_plane_vel)
                                        # split var names and load all vars
                                        cc=0
                                        bool_ar = np.ones(i_d.plane_cells,dtype=np.float64)
                                        if (i_d.POD_planes_rout != 'false'): # set cells outside taget radii to zero
                                                		r,theta,rc,thetac = GetPolarCoordinates(p2c.GetOutput(),i_d.t_o)
                                                		for ii in range (0,i_d.plane_cells):
                                                        		if (rc[ii]<i_d.POD_planes_rin[j]):
                                                                		bool_ar[ii] = 0
                                                        		if (rc[ii]>i_d.POD_planes_rout[j]):
                                                                		bool_ar[ii] = 0
                                        for jj in range (0,i_d.pp_vcount):
                                                #print i_d.pp_var_name_list[jj], i_d.vel_var
                                        	if i_d.pp_var_name_list[jj]== i_d.vel_var:                                      	
                                        		u_temp_cell = VN.vtk_to_numpy(p2c.GetOutput().GetCellData().GetVectors(i_d.vel_var))
                                                        #print np.mean(u_temp_cell)
                                        		if (i_d.make_planes_polar == 'true'):                                     
                                                	
                                                		POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] = u_temp_cell[:,0]*bool_ar   
                                                		cc+=1
                                                		POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] =  \
                                                        	(u_temp_cell[:,1]*np.cos(thetac) +  u_temp_cell[:,2]*np.sin(thetac))*bool_ar
                                                		cc+=1
                                                		POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] =  \
                                                        	(-u_temp_cell[:,1]*np.sin(thetac) +  u_temp_cell[:,2]*np.cos(thetac))*bool_ar  
                                        			cc+=1
                                        		else:
                                                		for i in range(0,i_d.num_components):
                                                        		POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] = \
                                                        		u_temp_cell[:,i]*bool_ar
                                                        		cc+=1
                                                elif i_d.pp_var_name_list[jj]== i_d.te_var:	
                                          		u_temp_cell = VN.vtk_to_numpy(p2c.GetOutput().GetCellData().GetScalars(i_d.te_var))
                                                        POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] = u_temp_cell*bool_ar
                                          		cc+=1
                                                elif i_d.pp_var_name_list[jj]== i_d.ed_var:
                                          		u_temp_cell = VN.vtk_to_numpy(p2c.GetOutput().GetCellData().GetScalars(i_d.ed_var))
                                          		POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] = u_temp_cell*bool_ar
                                          		cc+=1
                                                
                                          		

                if (len(i_d.save_planes)>len(i_d.POD_planes)):
                        for j in range (i_d.plane_count,len(i_d.save_planes)):
                                if (i_d.save_planes[j] == k):
                                        if not os.path.exists(i_d.snapshot_dir+'post/'+i_d.plane_names[k]):
                                                os.mkdir(i_d.snapshot_dir+'post/'+i_d.plane_names[k])

                                        n1=i_d.plane_normals[k,0]
                                        n2=i_d.plane_normals[k,1]
                                        n3=i_d.plane_normals[k,2]
                                        o1=i_d.plane_centres[k,0]
                                        o2=i_d.plane_centres[k,1]
                                        o3=i_d.plane_centres[k,2]
                                        cutMapper, pl3d_output, src_output1 = extract_plane(n1,n2,n3,o1,o2,o3,extract,i_d)
                                        #plt.vtkcontourf(cutMapper,pl3d_output,src_output,i_d.p_var,100,snapshot_dir+'/snap'+str(k),plane_normals[k,:],(0,0))
                                
                                        src_output = src_output1.GetOutput()

                                
                                        u = VN.vtk_to_numpy(src_output.GetPointData().GetVectors(i_d.vel_var))
                                        p = VN.vtk_to_numpy(src_output.GetPointData().GetScalars(i_d.p_var))
                                        tke = VN.vtk_to_numpy(src_output.GetPointData().GetScalars(i_d.te_var))
                                        epsi = VN.vtk_to_numpy(src_output.GetPointData().GetScalars(i_d.ed_var))
                                                
                                        points = VN.vtk_to_numpy(src_output.GetPoints().GetData())
                                
                                        #print np.shape(u_temp_read),np.shape(points)
                                        timestring = "%.5E" % time
                        
                                        filename = i_d.snapshot_dir+'post/'+i_d.plane_names[k]+'/'+timestring+'.prf'

                                        target = open(filename,'w')

                                        target.write('# ' + i_d.plane_names[k] + ' # name of the profile\n')
                                        target.write('# turbulence model, ' + i_d.turbulence_model+'\n')
                                        plane_rhs = (o1+i_d.t_o[0])*n1+(o2+i_d.t_o[1])*n2+(o3+i_d.t_o[2])*n3
                                        target.write('# plane normal and translation ' + str(n1)+'\t'+str(n2)+'\t'+str(n3)+'\t'+str(plane_rhs)+'\n')
                                        #target.write(str(src_output.GetNumberOfPoints())+'\n')
                                
                                        target.write('type, xyz # type of profile (rad or xyz)\n')
                                        target.write('localcs,origin,0,0,0 # origin of local coordinate system\n')
                                        target.write('localcs,xaxis,1,0,0 # x axis direction of local coordinate system\n')
                                        target.write('localcs,yaxis,0,1,0 # y axis direction of local coordinate system\n')
                                        target.write('localcs,zaxis,0,0,1 # z axis direction of local coordinate system\n')
                                        target.write('tolerance, 1.00E-08 # tolerance\n')
                                        target.write('scale,1,1,1,1,1,1,1,1 # scaling factors\n')
                                        if (i_d.turbulence_model == 'k_epsilon'):
                                                target.write('data,x,y,z,u,v,w,k,e\n')
                                        if (i_d.turbulence_model == 'k_omega'):
                                                target.write('data,x,y,z,u,v,w,k,sdr\n')
                                        for i in range (0,src_output.GetNumberOfPoints()):
                                                target.write(str(points[i,0])+','+str(points[i,1])+','+str(points[i,2])+','+str(u[i,0])+','+str(u[i,1])+','+str(u[i,2]) + ',' + str(tke[i])+','+str(epsi[i])+'\n')
                                        #target.write("%9.9f"% time_vec[i])
                                        #target.write('\t')
                                        #target.write(('\t'.join(map(str,pprobe[:,i])))+'\n') #for j in range (0,num_probes)))
                                        #target.write('\n'.join(' '.join(str(pprobe[j,i])) for j in range (0,num_probes))
                                        
                                        target.close()


                if (len(i_d.save_planes)<len(i_d.POD_planes)):
                        for j in range (i_d.plane_count,len(i_d.POD_planes)):
                                if (i_d.POD_planes[j] == k):


                                        n1=i_d.plane_normals[k,0]
                                        n2=i_d.plane_normals[k,1]
                                        n3=i_d.plane_normals[k,2]
                                        o1=i_d.plane_centres[k,0]
                                        o2=i_d.plane_centres[k,1]
                                        o3=i_d.plane_centres[k,2]
                                        cutMapper, pl3d_output, src_output1 = extract_plane(n1,n2,n3,o1,o2,o3,extract,i_d)
                                        #plt.vtkcontourf(cutMapper,pl3d_output,src_output,i_d.p_var,100,snapshot_dir+'/snap'+str(k),plane_normals[k,:],(0,0))
                                
                                        src_output = src_output1.GetOutput()

                                        POD_plane_areas[j,:],cell_centres = calc_cell_areas(src_output1.GetOutput())
                                        p2c = vtk.vtkPointDataToCellData()
                                        p2c.SetInputConnection(src_output1.GetOutputPort())
                                        p2c.Update()
                                        #print np.shape(POD_plane_vel)
                                        # split var names and load all vars
                                        cc=0
					bool_ar = np.ones(i_d.plane_cells,dtype=np.float64)
                                        if (i_d.POD_planes_rout != 'false'): # set cells outside taget radii to zero
                                                		r,theta,rc,thetac = GetPolarCoordinates(p2c.GetOutput(),i_d.t_o)
                                                		for ii in range (0,i_d.plane_cells):
                                                        		if (rc[ii]<i_d.POD_planes_rin[j]):
                                                                		bool_ar[ii] = 0
                                                        		if (rc[ii]>i_d.POD_planes_rout[j]):
                                                                		bool_ar[ii] = 0
                                        for jj in range (0,i_d.pp_vcount):
                                        	if i_d.pp_var_name_list[jj]== i_d.vel_var:                                      	
                                        		u_temp_cell = VN.vtk_to_numpy(p2c.GetOutput().GetCellData().GetVectors(i_d.vel_var))
                                        		if (i_d.make_planes_polar == 'true'):                                     
                                                	
                                                		POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] = u_temp_cell[:,0]*bool_ar   
                                                		cc+=1
                                                		POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] =  \
                                                        	(u_temp_cell[:,1]*np.cos(thetac) +  u_temp_cell[:,2]*np.sin(thetac))*bool_ar
                                                		cc+=1
                                                		POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] =  \
                                                        	(-u_temp_cell[:,1]*np.sin(thetac) +  u_temp_cell[:,2]*np.cos(thetac))*bool_ar  
                                        			cc+=1
                                        		else:
                                                		for i in range(0,i_d.num_components):
                                                        		POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] = \
                                                        		u_temp_cell[:,i]*bool_ar
                                                        		cc+=1
                                                elif i_d.pp_var_name_list[jj]== i_d.te_var:	
                                          		u_temp_cell = VN.vtk_to_numpy(p2c.GetOutput().GetCellData().GetScalars(i_d.te_var))
                                                        POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] = u_temp_cell*bool_ar
                                          		cc+=1
                                                elif i_d.pp_var_name_list[jj]== i_d.ed_var:
                                          		u_temp_cell = VN.vtk_to_numpy(p2c.GetOutput().GetCellData().GetScalars(i_d.ed_var))
                                          		POD_plane_vel[j,cc*i_d.plane_cells:(cc+1)*i_d.plane_cells,snap] = u_temp_cell*bool_ar
                                          		cc+=1
#=======================================================================
def save_plane(u,i_d,i):
        
        cc = vtk.vtkCellCenters()
        cc.SetInputData(i_d.grid)
        cc.VertexCellsOn()
        cc.Update()
        points = VN.vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
        npt = cc.GetOutput().GetNumberOfPoints()

        #print np.shape(u_temp_read),np.shape(points)
#        timestring = "%.5E" % i_d.time
#        
#        filename = './PODFS/'+timestring+'.prf'      
#        
#        target = open(filename,'w')
#
#        target.write('# Generated using the digital filter method # name of the profile\n')
#        target.write('# turbulence model, ' + 'none\n')
#        plane_rhs = (i_d.t_o[0])*i_d.n[0]+(i_d.t_o[1])*i_d.n[1]+(i_d.t_o[2])*i_d.n[2]
#        target.write('# plane normal and translation ' + str(i_d.n[0])+'\t'+str(i_d.n[1])+'\t'+str(i_d.n[2])+'\t'+str(plane_rhs)+'\n')
#        #target.write(str(src_output.GetNumberOfPoints())+'\n')
#                                
#        target.write('# type, xyz # type of profile (rad or xyz)\n')
#        target.write('# localcs,origin,0,0,0 # origin of local coordinate system\n')
#        target.write('# localcs,xaxis,1,0,0 # x axis direction of local coordinate system\n')
#        target.write('# localcs,yaxis,0,1,0 # y axis direction of local coordinate system\n')
#        target.write('# localcs,zaxis,0,0,1 # z axis direction of local coordinate system\n')
#        target.write('# tolerance, 1.00E-08 # tolerance\n')
#        target.write('# scale,1,1,1,1,1,1 # scaling factors\n')
#        target.write('data,x,y,z,u,v,w\n')
#        for i in range (0,npt):
#                target.write(sp.str(points[i,0])+','+sp.str(points[i,1])+','+sp.str(points[i,2])+','+sp.str(u[i])+','+sp.str(u[i+npt])+','+sp.str(u[i+2*npt]) + '\n')
#                                        
#        target.close()
        
        
        # export for cfx
        u_temp=np.zeros(npt)
        v_temp=np.zeros(npt)
        w_temp=np.zeros(npt)
        
        for l in range (0,npt):
            u_temp[l]=u[l]
            v_temp[l]=u[l+npt]
            w_temp[l]=u[l+2*npt]
            
        Vel_mag=np.sqrt(u_temp**2+v_temp**2+w_temp**2)
        
        u_comp=u_temp/Vel_mag
        v_comp=v_temp/Vel_mag
        w_comp=w_temp/Vel_mag
        
        
        
        timestring = "%.5E" % i_d.time
        timestepnum = str(i+1)
        
        filename = './PODFS/'+timestring+'.csv'      
        
        target = open(filename,'w')

        target.write('[Name]\n')
        target.write('turbInlet\n')
        target.write('[Spatial Fields]\n')
        target.write('x, y, z\n')
        target.write('[Data]\n')
        target.write('x [ m ],y [ m ],z [ m ],u [ ],v [ ],w [ ]\n')
        for i in range (0,npt):
                target.write(sp.str(points[i,0])+','+sp.str(points[i,1])+','+sp.str(points[i,2])+','+sp.str(u_comp[i])+','+sp.str(v_comp[i])+','+sp.str(w_comp[i]) + '\n')
                                        
        target.close()       
                                        
        filename = './PODFS/inlet.dat'+timestepnum      
        
        target = open(filename,'w')
        
        for m in range (0,npt):
                target.write(sp.str(points[m,0])+'\t'+sp.str(points[m,1])+'\t'+sp.str(points[m,2])+'\t'+sp.str(u_comp[m])+'\t'+sp.str(v_comp[m])+'\t'+sp.str(w_comp[m]) + '\t'+sp.str(w_comp[m]) +'\n')
                #target.write(sp.str(points[m,0])+'\t'+sp.str(points[m,1])+'\t'+sp.str(points[m,2])+'\t100000\t300\n')                       
        target.close() 
#=======================================================================
def calc_cell_areas(extract):

        num_cells = extract.GetNumberOfCells()
        cell_areas = np.zeros((num_cells), dtype = np.float64)
        cell_centres = np.zeros((num_cells,3), dtype = np.float64)
        # for each cell find cell type, and cell volume
        num_tri = 0
        num_quad = 0
        num_unk = 0
        for i in range (0,num_cells):
            cell_type=extract.GetCellType(i)
            points = extract.GetCell(i).GetPoints()

            if (cell_type == 5): #then triangle 
                num_tri += 1
                np_pts = np.array([points.GetPoint(j) for j in range(points.GetNumberOfPoints())])
                a = np.linalg.norm(np_pts[0,:]-np_pts[1,:])
                b = np.linalg.norm(np_pts[1,:]-np_pts[2,:])
                c = np.linalg.norm(np_pts[2,:]-np_pts[0,:])
                s = 0.5*(a+b+c)
                cell_areas[i] = np.sqrt(s*(s-a)*(s-b)*(s-c))
                cell_centres[i,:] = (np_pts[0,:]+np_pts[1,:]+np_pts[2,:])/3.

            elif(cell_type == 9): # then quad
                num_quad += 1  
                np_pts = np.array([points.GetPoint(j) for j in range(points.GetNumberOfPoints())])
                a = np.linalg.norm(np_pts[0,:]-np_pts[1,:])
                b = np.linalg.norm(np_pts[1,:]-np_pts[2,:])
                c = np.linalg.norm(np_pts[2,:]-np_pts[0,:])
                s = 0.5*(a+b+c)
                cell_areas[i] = np.sqrt(s*(s-a)*(s-b)*(s-c))

                a = np.linalg.norm(np_pts[0,:]-np_pts[2,:])
                b = np.linalg.norm(np_pts[2,:]-np_pts[3,:])
                c = np.linalg.norm(np_pts[3,:]-np_pts[0,:])
                s = 0.5*(a+b+c)
                cell_areas[i] = cell_areas[i]+np.sqrt(s*(s-a)*(s-b)*(s-c)) 
                cell_centres[i,:] = (np_pts[0,:]+np_pts[1,:]+np_pts[2,:]+np_pts[3,:])/4.

            #else: # then unknown
             #   print 'unknown cell, id: ' cell_type

        cell_vol_filt = cell_areas[np.nonzero(cell_areas)[0]]
        mean_vol = np.mean(cell_vol_filt)
        min_vol = np.amin(cell_vol_filt)
        max_vol = np.amax(cell_vol_filt)

        print 'Plane Area Statistics:'
        print 'Number of cells: ', num_cells
        print 'Number of quads: ', num_quad
        print 'Number of triangles: ', num_tri
        print 'Smallest cell volume: ', min_vol
        print 'Mean cell volume: ', mean_vol
        print 'Largest cell volume: ', max_vol

        return cell_areas,cell_centres

#=======================================================================
def save_instants(grid,num_points,num_snapshots,scalars,vectors,scalar_names,vector_names, \
                  extract_planes,plane_names,plane_normals,plane_centres,x_min,\
                  x_max,y_min,y_max,z_min,z_max, \
                  snapshot_dir,file_list,num_components,t_o,scalar_range, \
                  extract_cylinders,cyl_names,cyl_normals,cyl_radius,cyl_centres,roi,i_d):
        
        if not os.path.exists(snapshot_dir):
                os.mkdir(snapshot_dir)
                                
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileTypeToBinary()

        u_temp = np.array(np.zeros((num_points,3), dtype=np.float64))
        scalar_temp = np.array(np.zeros((num_points), dtype=np.float64))
	grid.GetPointData().Initialize()
        filter1 = grid

        print 'ok'

        for j in range(0,num_snapshots):
		print '      Saving Instantaneous Data to: ', file_list[j].strip(), 'file number ', j+1, ' of ', num_snapshots
                counter = "%9.9d"% j
                # Add pdash and udash
                for i in range(0,len(vector_names)):
                        for k in range(0,num_components):
                                u_temp[:,k] = vectors[k*num_points:(k+1)*num_points,j,i]
                        u_vtk = VN.numpy_to_vtk(u_temp,1)
                        u_vtk.SetName(vector_names[i])
                        filter1.GetPointData().AddArray(u_vtk)
                for i in range(0,len(scalar_names)):
                        scalar_temp[:] = scalars[:,j,i]
                        u_vtk = VN.numpy_to_vtk(scalar_temp,1)
                        u_vtk.SetName(scalar_names[i])
                        filter1.GetPointData().AddArray(u_vtk)

                # Extract planes and save image
                levels = 200
                if (extract_planes == 'true' and i_d.plot_vtk) :

                        print '            Extracting and Saving Planes...'

                        for k in range (0,len(plane_names)):

                                if not os.path.exists(snapshot_dir+plane_names[k]):
                                        os.mkdir(snapshot_dir+plane_names[k])

                                n1=plane_normals[k,0]
                                n2=plane_normals[k,1]
                                n3=plane_normals[k,2]
                                o1=plane_centres[k,0]
                                o2=plane_centres[k,1]
                                o3=plane_centres[k,2]
                                cutMapper, pl3d_output, src_output = extract_plane(n1,n2,n3,o1,o2,o3,filter1,i_d)

                                for i in range (0,len(scalar_names)):
                                        sr = (scalar_range[0,i],scalar_range[1,i])

                                        if not os.path.exists(snapshot_dir+plane_names[k]+'/'+i_d.scalar_save_names[i]):
                                                os.mkdir(snapshot_dir+plane_names[k]+'/'+i_d.scalar_save_names[i])
                                                                
                
                                        plt.vtkcontourf(cutMapper,pl3d_output,src_output,scalar_names[i],levels,snapshot_dir+plane_names[k]+'/'+i_d.scalar_save_names[i]+'/snap'+counter,plane_normals[k,:],sr)

                if (extract_cylinders == 'true' and i_d.plot_vtk):

                        print "\n   Extracting Cylinders..."

                        #var_name=scalar_names[4]
                        levels = 200

                        for k in range (0,len(cyl_names)):

                                if not os.path.exists(snapshot_dir+cyl_names[k]):
                                        os.mkdir(snapshot_dir+cyl_names[k])

                                        # print cyl_normals[k,:],cyl_radius[k],cyl_centres[k,:]

                                cutMapper, pl3d_output, src_output = extract_cylinder(cyl_normals[k,:],cyl_radius[k],cyl_centres[k,:],roi,filter1,t_o)

                                for i in range (0,len(scalar_names)):
               
                                        if not os.path.exists(snapshot_dir+cyl_names[k]+'/'+scalar_save_names[i]):
                                                os.mkdir(snapshot_dir+cyl_names[k]+'/'+scalar_save_names[i])
                                                                
               
                                        plt.vtkcontourf(cutMapper,pl3d_output,src_output,scalar_names[i],levels,snapshot_dir+cyl_names[k]+'/'+i_d.scalar_save_names[i]+'/snap'+counter,[0,0,1],(0,0))


                                        
                # Save vtk
                writer.SetInputData(filter1)
                writer.SetFileName(snapshot_dir + file_list[j].strip())
                writer.Write()

#=======================================================================
def write_text_file(lhs,M,filename):

        if not os.path.exists('./post/textfiles'):
                os.mkdir('./post/textfiles')

	print '\n   Writing text file to:'
        output_file_name = './post/textfiles/'+filename+'.dat'
        print '      ', output_file_name
        n = np.size(lhs)
        m = np.ndim(M)
        target = open(output_file_name,'w')
        for i in range(0,n):
                target.write("%9.9f"% lhs[i])
		target.write('\t')
                if (m > 1):
                        target.write(('\t'.join(map(str,M[:,i])))+'\n')
                else:
                        target.write("%9.9f"% M[i])
                        target.write('\n')
        target.close()                 
#=======================================================================
def write_text_file2(lhs,M,dire,filename):

        if not os.path.exists(dire+'textfiles'):
                os.mkdir(dire+'textfiles')

	print '\n   Writing text file to:'
        output_file_name = dire+'textfiles/'+filename+'.dat'
        print '      ', output_file_name
        n = np.size(lhs)
        m = np.ndim(M)
        target = open(output_file_name,'w')
        for i in range(0,n):
                target.write("%9.9f"% lhs[i])
		target.write('\t')
                if (m > 1):
                        target.write(('\t'.join(map(str,M[:,i])))+'\n')
                else:
                        target.write("%9.9f"% M[i])
                        target.write('\n')
        target.close()                 

#======================================================================
def GetPolarCoordinates(grid,t_o):

        #print t_o

        num_points = grid.GetNumberOfPoints()
        num_cells = grid.GetNumberOfCells()
        #print num_cells

        r= np.zeros((num_points),dtype=np.float64)
        theta= np.zeros((num_points),dtype=np.float64)
        rc= np.zeros((num_cells),dtype=np.float64)
        thetac= np.zeros((num_cells),dtype=np.float64)

        for i in range (0,num_points):
                point = grid.GetPoints().GetPoint(i)
                y = point[1]-t_o[1]
                z = point[2]-t_o[2]
                r[i] = np.sqrt(y**2+z**2)
                if (r[i] < np.finfo(np.float32).eps):
                        r[i] = np.finfo(np.float32).eps
                theta[i] = np.arctan2(z,y)

        for i in range (0,num_cells):
                cell_type=grid.GetCellType(i)
                points = grid.GetCell(i).GetPoints()
                np_pts = np.array([points.GetPoint(j) for j in range(points.GetNumberOfPoints())])
                y = np.mean(np_pts[:,1])-t_o[1]
                z = np.mean(np_pts[:,2]) -t_o[2]
                #print y,z
                rc[i] = np.sqrt(y**2+z**2)
                if (rc[i] < np.finfo(np.float32).eps):
                        rc[i] = np.finfo(np.float32).eps
                thetac[i] = np.arctan2(z,y)

        return r,theta,rc,thetac

#=====================================================================
def calculate_cell_volume(cell):

   cell_type=cell.GetCellType() #extract.GetOutput().GetCellType(i)
   points = cell.GetPoints() #extract.GetOutput().GetCell(i).GetPoints()

   if (cell_type == 10): #then tetrahedral
           #num_tet += 1
           
        np_pts = np.array([points.GetPoint(j) for j in range(points.GetNumberOfPoints())])
        
        l = np.linalg.norm(np_pts[2,:]-np_pts[1,:])
        j = np_pts[3,:]-np_pts[0,:]
        
        r = np_pts[1,:]-np_pts[0,:]
        s = np_pts[2,:]-np_pts[0,:]
        h = np.linalg.norm(s)
        e = np.linalg.norm(r)
        
        s1 = 0.5*(h+l+e)
           
        A1 = np.sqrt(s1*(s1-h)*(s1-l)*(s1-e))
        
        n = np.cross(r,s)
        h = np.linalg.norm(np.dot(j,(n/np.linalg.norm(n))))

        cell_volume = A1*h/3.0


   elif(cell_type == 11): # then voxel
        #num_vox +=  1
        print 'Error: Voxel detected, Cannot Calculate Volume!'
        jnsglo
 
   elif(cell_type == 12): # then hexahedron
        #num_hex += 1
        np_pts = np.array([points.GetPoint(j) for j in range(points.GetNumberOfPoints())])
        # contruct matrices
        M1=np.ones((4,4))
        M2=np.ones((4,4))
        M3=np.ones((4,4))
        M4=np.ones((4,4))
        
        M1[:,1:4]=[np_pts[0,:],np_pts[1,:],np_pts[2,:],np_pts[5,:]]
        M2[:,1:4]=[np_pts[2,:],np_pts[5,:],np_pts[6,:],np_pts[7,:]]
        M3[:,1:4]=[np_pts[0,:],np_pts[2,:],np_pts[3,:],np_pts[7,:]]
        M4[:,1:4]=[np_pts[0,:],np_pts[4,:],np_pts[5,:],np_pts[7,:]]
        
        #find determinants
        v1=np.abs(np.linalg.det(M1))
        v2=np.abs(np.linalg.det(M2))
        v3=np.abs(np.linalg.det(M3))
        v4=np.abs(np.linalg.det(M4))
        
        # calculate volume
        cell_volume = (v1+v2+v3+v4)/6.0
         
   elif(cell_type == 13): # then wedge
        #num_wedge += 1
        np_pts = np.array([points.GetPoint(j) for j in range(points.GetNumberOfPoints())])
        # contruct matrices
        M1=np.ones((4,4))
        M2=np.ones((4,4))
        M3=np.ones((4,4))
        
        M1[:,1:4]=[np_pts[0,:],np_pts[1,:],np_pts[2,:],np_pts[4,:]]
        M2[:,1:4]=[np_pts[0,:],np_pts[2,:],np_pts[3,:],np_pts[4,:]]
        M3[:,1:4]=[np_pts[2,:],np_pts[3,:],np_pts[4,:],np_pts[5,:]]
        
        #find determinants
        v1=np.abs(np.linalg.det(M1))
        v2=np.abs(np.linalg.det(M2))
        v3=np.abs(np.linalg.det(M3))
        
        # calculate volume
        cell_volume = (v1+v2+v3)/6.0
        
   elif(cell_type == 14): # then pyramid
        #num_pyr += 1
        np_pts = np.array([points.GetPoint(j) for j in range(points.GetNumberOfPoints())])
        #plt.scatter_3d(1,np_pts[:,0],np_pts[:,1],np_pts[:,2],['1','2','3','4','5'])
        h = np.linalg.norm(np_pts[3,:]-np_pts[2,:])
        e = np.linalg.norm(np_pts[3,:]-np_pts[0,:])
        g = np.linalg.norm(np_pts[2,:]-np_pts[1,:])
        f = np.linalg.norm(np_pts[1,:]-np_pts[0,:])
        
        r = np_pts[1,:]-np_pts[3,:]
        s = np_pts[2,:]-np_pts[0,:]
        l = np.linalg.norm(s)
        j = np_pts[4,:] - np_pts[2,:]
           
        s1 = 0.5*(h+l+e)
        s2 = 0.5*(g+f+l)
           
        A1 = np.sqrt(s1*(s1-h)*(s1-l)*(s1-e))
        A2 = np.sqrt(s2*(s2-g)*(s2-f)*(s2-l))

        A=A1+A2

        n = np.cross(r,s)
        h = np.linalg.norm(np.dot(j,(n/np.linalg.norm(n))))

        cell_volume = A*h/3.0

   return cell_volume

#=======================================================================
def read_inflow_data(i_d):

 A = np.zeros((i_d.num_cells*i_d.num_components,i_d.num_snapshots),dtype=np.float64)

 for i in range (0,i_d.num_snapshots):
	# load file
	cstring = '%6.6i' % i
	M = np.loadtxt('inflow.'+cstring)
	A[:,i] = M.reshape(i_d.num_cells*i_d.num_components,order='F')
 i_d.A = A
 i_d.cell_volume = np.ones(i_d.num_cells,dtype=np.float64)
	
#=======================================================================
def make_inflow_plane(i_d):
    planeS = vtk.vtkPlaneSource()
    planeS.SetResolution(i_d.kma,i_d.jma)
    nx = i_d.n[0]
    ny = i_d.n[1]
    nz = i_d.n[2]

    # find angles between normal and principle directions
    alpha = np.arccos(nx)*180/np.pi # angle between plane normal and [1,0,0]
    beta = np.arctan2(nz,ny)*180/np.pi # twist correction
    
    # scale plane
    s1 = 0 
    s2 = i_d.res*float(i_d.jma)*float(i_d.jma)/(float(i_d.jma)-1) 
    s3 = i_d.res*float(i_d.kma)*float(i_d.kma)/(float(i_d.kma)-1) 

    # perform transformations
    trans = vtk.vtkTransform()
    trans1 = vtk.vtkTransform()
    trans2 = vtk.vtkTransform()
    trans3 = vtk.vtkTransform()
    planeS.SetNormal(1, 0, 0)
    trans1.Scale(s1,s2,s3)
    trans2.RotateWXYZ(alpha,0,-nz,ny)
    trans3.RotateWXYZ(beta+i_d.rot,nx,ny,nz)
    #trans3.RotateY(beta)
    #trans2.RotateX(alpha)
    print nx,ny,nz
    o1 = i_d.t_o[0]
    o2 = i_d.t_o[1]
    o3 = i_d.t_o[2]
    trans.Translate(o1,o2,o3)

    transf1 = vtk.vtkTransformPolyDataFilter()
    transf1.SetInputConnection(planeS.GetOutputPort())
    transf1.SetTransform(trans1)
    transf2 = vtk.vtkTransformPolyDataFilter()
    transf2.SetInputConnection(transf1.GetOutputPort())
    transf2.SetTransform(trans2)
    transf3 = vtk.vtkTransformPolyDataFilter()
    transf3.SetInputConnection(transf2.GetOutputPort())
    transf3.SetTransform(trans3)
    transf = vtk.vtkTransformPolyDataFilter()
    transf.SetInputConnection(transf3.GetOutputPort())
    transf.SetTransform(trans)
    transf.Update()
    extractS = transf.GetOutput()
    return extractS


#=======================================================================
def POD(A,num_snapshots,num_points,num_components,correct_for_cell_volumes,cell_volume,restart_dir,restart_flag,tol_CN,num_modes_trunc,num_modes_to_write,test_POD_orthogonality,write_matrices,grid,mean_field,dt,var_name,ifig,N,iwindow,stride,i_d):


        print '\n=============================================================='
        print 'Calculating POD modes ...'

        C = np.array(np.zeros((num_snapshots,num_snapshots), dtype=np.float64))
        print '\n   Calculating covariance matrix ...'
        
        calculate_correlation_matrix(num_snapshots,num_points,num_components,correct_for_cell_volumes,cell_volume,A,C)
        
        
        print '\n   Solving eigenvalue problem ...'
        energy = np.array(np.zeros((num_snapshots), dtype=np.complex64))
        temporal_modes = np.array(np.zeros((num_snapshots,num_snapshots), dtype=np.complex64))
        energy,temporal_modes = linalg.eig(C)
        sort_eigenvalues(num_snapshots,energy,temporal_modes)

        num_valid_modes = 0
        while ( (energy[num_valid_modes].real/energy[0].real>pow(tol_CN,2.0)) and (num_valid_modes<num_snapshots-2) \
		and(energy[num_valid_modes].real>0.0) ):
                num_valid_modes += 1
                if ( (energy[num_valid_modes].real/energy[0].real > pow(tol_CN,2.0)) and (energy[num_valid_modes].real > 0.0) ):
                        num_valid_modes += 1
        print '      Number of valid POD modes with positive energies = ', num_valid_modes
        if ( (num_modes_trunc < 0) or (num_modes_trunc > num_valid_modes) ):
                num_modes_trunc = num_valid_modes

        print '\n   Scaling temporal modes ...'
        for j in range(0,num_valid_modes):
                temporal_mode_mag = sum(temporal_modes[:,j].real * temporal_modes[:,j].real)/num_snapshots
                temporal_modes[:,j] = temporal_modes[:,j] * np.sqrt(energy[j].real/temporal_mode_mag)
        initial_conditions = np.array(np.zeros((num_valid_modes), dtype=np.float64))
        initial_conditions = temporal_modes[0,:].real

        print '\n   Calculating truncated spatial modes ...'
        energy_trunc_inv = np.array(np.zeros((num_modes_trunc,num_modes_trunc), dtype=np.float64))
        energy_trunc_inv = np.diag(np.ones(num_modes_trunc)/energy[0:num_modes_trunc].real,0)
        spatial_modes_trunc = np.array(np.zeros((num_components*num_points,num_modes_trunc), dtype=np.float64))
        spatial_modes_trunc = np.dot( np.dot(A[:,0:num_snapshots], temporal_modes[:,0:num_modes_trunc].real), energy_trunc_inv) / num_snapshots

        #=======================================================================
        print '\n=============================================================='
        print 'Writing POD output ...'

        write_eigenvalues(num_valid_modes,num_snapshots,energy,restart_dir+'POD.eigenvalues.dat')

        write_mean_field2(grid,mean_field,var_name,num_components,restart_dir)
        del mean_field


	# scale modes
	energy_trunc = np.sum(energy[0:num_modes_trunc].real)
	energy_total = np.sum(energy.real)
	scale1 = np.sqrt(energy_total/energy_trunc)
	#spatial_modes_trunc = spatial_modes_trunc*scale1

        
        if i_d.verbose:


        
         write_temporal_modes(num_valid_modes,num_snapshots,dt,temporal_modes,restart_dir)

         num_POD_modes_to_write = num_modes_to_write
         if ( (num_modes_to_write < 0) or (num_modes_to_write > num_modes_trunc) ):
                num_POD_modes_to_write = num_modes_trunc
         print '\n   Number of spatial POD modes to write = ', num_POD_modes_to_write
         write_spatial_POD_modes_i_d(num_POD_modes_to_write,num_components,grid,spatial_modes_trunc,var_name,restart_dir,i_d)


        #======================================================================
        if i_d.verbose:
         #Plot things!
         # calculate PSD of temporal mode
         time_vec = dt*stride*np.linspace(0,num_snapshots,num_snapshots)
         modes=np.linspace(0,num_snapshots,num_snapshots)
         Sxx = np.array(np.zeros((N,num_POD_modes_to_write), dtype=np.float64))
         for j in range(0,num_POD_modes_to_write):
                fs,Sxx[:,j],M = sp.fct_welch(temporal_modes[:,j],1./dt/stride,N,iwindow)

         ifig+=1
         #plt.eigs(ifig,energy,'\lambda',restart_dir+'POD_mode_energies')

	        
         #plot modes
         for j in range(0,num_POD_modes_to_write):
                ifig=ifig+1
                #plt.timeseries(ifig,temporal_modes[:,j],time_vec,'\,',restart_dir+'POD_time_mode_'+str(j))
                plt.close(ifig)
                ifig=ifig+1
                plt.PSD(ifig,Sxx[:,j],fs,'Power',restart_dir+'POD_PSD_mode_'+str(j))
                plt.close(ifig)

        print '\nCode completed.'
        print '==============================================================\n'  
	
        i_d.temporal_modes = temporal_modes
        i_d.spatial_modes = spatial_modes_trunc
        i_d.nm = num_modes_trunc
	
#=======================================================================
def write_initial_conditions(num_valid_modes,initial_conditions,rdir):
	print '\n   Writing initial conditions to'
	filename = rdir+'POD.initial_conditions.dat'
	print '      ', filename
	file = open(filename,'w')
	file.write('#\n')
	file.write('# mode, amplitude\n')
	file.write('#\n')
	for i in range(0,num_valid_modes):
		file.write('%4.1d %18.10e\n' % (i+1, initial_conditions[i]) )
	file.close()

#=======================================================================
def write_eigenvalues(num_valid_modes,num_snapshots,energy,filename):
	cumulative_energy = np.array(np.zeros((num_valid_modes), dtype=np.float64))
	cumulative_energy[0] = energy[0].real
	for i in range(1,num_valid_modes):
		cumulative_energy[i] = cumulative_energy[i-1] + energy[i].real
	total_energy = cumulative_energy[num_valid_modes-1]

	print '\n   Writing energies to'
	print '      ', filename
	file = open(filename,'w')
	file.write('#\n')
	file.write('# mode, energy, cumulative, percenterage energy, percentage cumulative, condition number (absolute value if negative)\n')
	file.write('#		Note: cummulative energies are set to zero after first negative energy')
	file.write('#\n')
	for i in range(0,num_valid_modes):
		file.write('%4.1d %18.10e %18.10e %18.10e %18.10e %18.10e\n' % (i+1, energy[i].real, cumulative_energy[i], energy[i].real/total_energy*100.0, cumulative_energy[i]/total_energy*100.0, math.sqrt(energy[i].real / energy[0].real)) )
	for i in range(num_valid_modes,num_snapshots):
		file.write('%4.1d %18.10e %18.10e %18.10e %18.10e %18.10e\n' % (i+1, energy[i].real, 0.0, energy[i].real/total_energy*100.0, 0.0, math.sqrt(abs(energy[i].real / energy[0].real))) )
	file.close()

#=======================================================================
def sort_eigenvalues(num_snapshots,energy,temporal_modes):
        energy_sorted = np.array(np.zeros((num_snapshots), dtype=np.float64))
        mode_index = np.array(np.zeros((num_snapshots), dtype=np.int))
        for k in range(0,num_snapshots):
                mode_index[k] = k
                if ( (math.isnan(energy[k].real)) or (math.isnan(energy[k].imag)) ):
                        energy_sorted[k] = -1.0e10
                        temporal_modes[:,k] = 0.0
                else:
                        energy_sorted[k] = energy[k].real
        energy_sorted[0:num_snapshots], mode_index[0:num_snapshots] = zip(*sorted(zip(energy_sorted[:], mode_index[:]),reverse=True))
        energy[0:num_snapshots] = energy_sorted[0:num_snapshots]

        temporal_modes0 = np.array(np.zeros((num_snapshots,num_snapshots), dtype=np.float64))
        temporal_modes0[0:num_snapshots,0:num_snapshots] = temporal_modes[0:num_snapshots,0:num_snapshots].real.copy()
        for k in range(0,num_snapshots):
                temporal_modes[0:num_snapshots,k] = temporal_modes0[0:num_snapshots,mode_index[k]]
        return


#=======================================================================
def calculate_correlation_matrix(num_snapshots,num_points,num_components,correct_for_cell_volumes,cell_volume,A,C):
#        print np.mean(A[:,:],0)
	C[:,:] = 0.0
	if (correct_for_cell_volumes=="false"):
		C[:,:] = np.dot( A[:,0:num_snapshots].T, A[:,0:num_snapshots]) / num_snapshots
	elif (correct_for_cell_volumes=="true"):
		for i in range(0,num_snapshots):
			print "      Correlations with snapshot ", i+1, " of ", num_snapshots
			for j in range(0,num_snapshots):
				if (j>=i):
					for k in range(0,num_components):
						C[i,j] = C[i,j] + np.dot( A[k*num_points:(k+1)*num_points,i].T * cell_volume[:], A[k*num_points:(k+1)*num_points,j] ) / num_snapshots
				if (j<i):
					C[i,j] = C[j,i]


#=======================================================================
def write_temporal_modes(num_valid_modes,num_snapshots,dt,temporal_modes,rdir):
	print '\n   Writing temporal modes to'
	#print temporal_modes[:,0]
	for j in range(0,num_valid_modes):
		j_index = j+1
		output_file_name = rdir+'POD.temporal_mode_' + '%04d'%j_index + '.dat'
		print '      ', output_file_name
		file = open(output_file_name,'w')
		file.write('#\n')
		file.write('# time, amplitude\n')
		file.write('#\n')
		for i in range(0,num_snapshots):
			file.write('%18.10e %18.10e\n' % (i*dt, temporal_modes[i,j].real) )

                file.close()

#=======================================================================
def read_temporal_modes(num_valid_modes,num_snapshots,dt,temporal_modes):
	print '\n   Reading temporal modes from'
	#print temporal_modes[:,0]
	for j in range(0,num_valid_modes):
		j_index = j+1
		output_file_name = './post/POD/POD.temporal_mode_' + '%04d'%j_index + '.dat'
		print '      ', output_file_name
		file = open(output_file_name,'r')
		dummy = file.readline()
                dummy = file.readline()
                dummy = file.readline()

		for i in range(0,num_snapshots):
                        line1=file.readline()
                        linesplit = line1.split()
                        temporal_modes[i,j] = float(linesplit[1])
		file.close()
#=======================================================================
def read_temporal_modes_dir(num_valid_modes,num_snapshots,dt,temporal_modes,dire):
	print '\n   Reading temporal modes from'
	#print temporal_modes[:,0]
	for j in range(0,num_valid_modes):
		j_index = j+1
		output_file_name = dire+'POD.temporal_mode_' + '%04d'%j_index + '.dat'
		print '      ', output_file_name
		file = open(output_file_name,'r')
		dummy = file.readline()
                dummy = file.readline()
                dummy = file.readline()

		for i in range(0,num_snapshots):
                        line1=file.readline()
                        linesplit = line1.split()
                        temporal_modes[i,j] = float(linesplit[1])
		file.close()


#=======================================================================
def fourier_coefficients(i_d):

 num_snapshots = i_d.ns
 dt = i_d.dt
 stride = 1
 num_fcs = num_snapshots
 first = 1
 num_modes = i_d.nm
 energy_target = i_d.et
 rdir = './PODFS/'
 temporal_modes = i_d.temporal_modes
 ifig = 1


 # make constants
 last = num_snapshots+first
 time= np.linspace(0,(num_snapshots-1)*dt*stride,num_snapshots)
 print num_snapshots
 period = time[-1]+(time[1]-time[0])


 #allocate arrays
 #temporal_modes = np.array(np.zeros((num_snapshots,num_modes), dtype=np.float64))
 y2 = np.array(np.zeros((num_snapshots,num_modes), dtype=np.float64))
 c = np.array(np.zeros((num_fcs,num_modes), dtype=np.complex64))
 c_count = np.array(np.zeros((num_modes), dtype=np.int64))
 c_ind = np.zeros((num_modes,num_fcs),dtype=np.int32)

 print 'dt check:', (time[1]-time[0])

 #if test:
#	for i in range (0,num_modes):
#		temporal_modes[:,i] = -1 + np.sin(i*time/period*2*np.pi)+np.sin(i*time/period*4*np.pi)
# else:
        #load temporal coefficients
#        read_temporal_modes_dir(num_modes,num_snapshots,dt,temporal_modes,rdir)

 ifig = 0
 #for each mode calculate fourier coefficients
 for i in range (0,num_modes):
	y = temporal_modes[:,i]
	for n in range(0,num_fcs):
		k = n - num_fcs/2
		ctemp = y*np.exp(-1j*2*k*np.pi*time/period)
		#if n == 0 or n ==num_fcs-1:
		#if n == num_fcs/2:
                #        c[n,i]=ctemp.sum()/ctemp.size/2	
                #else:
                #        c[n,i]=ctemp.sum()/ctemp.size
                c[n,i]=ctemp.sum()/ctemp.size
		c_ind[i,n]= n
        #print c[0,i],c[1,i],c[2,i],c[3,i]


	#rank in order of coefficient value
	cmod = np.abs(c[:,i])
	c_copy = np.array(np.zeros((num_fcs), dtype=np.complex64))
	c_copy[:] = c[:,i]
	#print np.amax(c[:,i])
	for j in range (0,num_fcs):
		c_ind[i,j]=j
	cmod,c_ind[i,:]=zip(*sorted(zip(cmod,c_ind[i,:]),reverse=True))


	#decide on how many should be used based on energy critierion
	energy = 0.0
	energy_sum = np.sum(np.abs(c[:,i]))
	c_count[i] = 0
	while energy < energy_sum*energy_target:
		energy += np.abs(c[c_ind[i,c_count[i]],i])
		c_count[i] += 1
		#print c_count, cmod[c_count], c_ind[c_count]

	#print c[0,i],c[1,i]
	#stop
	#_count=100
	print 'The number of fourier coeffients used for mode ', i+1,' is ',c_count[i] 


	#generate approximate solutions
	j=0	
	for x in time:
                f=0
		for n in range(0,c_count[i]):
			k = c_ind[i,n]-num_fcs/2
			f += c[c_ind[i,n],i]*np.exp(1j*2*k*np.pi*x/period) 
			#if j == 0:
			#	print c[c_ind[i,n],i],k
                y2[j,i] = f
		j+=1 
	
	if i_d.verbose:
                #plot for comparison
                ifig += 1
                #print y
                #print time
                #print y2[:,i]
                #print ifig
                #print i
                plt.timeseries(ifig,y,time,'\,',rdir+'POD_tmode_recon'+str(i))
                plt.timeseries(ifig,y2[:,i],time,'\,',rdir+'POD_tmode_recon'+str(i))

 #save

 print '\n Saving to PODFS.dat ...'

 filename = rdir+'PODFS.dat'

 target = open(filename,'w')

 target.write(str(num_modes))
 target.write('\n'+str(period))
 for i in range (0,num_modes):
	target.write('\n'+ str(i+1)+'\t'+str(c_count[i]))

 for i in range(0,num_modes):
	for j in range(0,c_count[i]):
		target.write('\n'+ str(c_ind[i,j]-num_fcs/2)+'\t'+str(c[c_ind[i,j],i].real)+'\t'+str(c[c_ind[i,j],i].imag))

 target.close()

#======================================================================
def pod2prf(i_d):

 rdir='./PODFS/'
 convert_from_polar = False #convert from polar to cartesian coordinates?
 var = 'velocity'
 num_modes = i_d.nm
 i_d.turbulence_model = 'none'
 i_d.t_o = np.array([0,0,0])
 grid = i_d.grid
 num_points = i_d.num_points
 n1=i_d.n[0]
 n2=i_d.n[1]
 n3=i_d.n[2]
 o1=0
 o2=0
 o3=0
 #i_d.x_min=10
 #i_d.x_max=-10
 #i_d.y_min=10
 #i_d.y_max=-10
 #i_d.z_min=10
 #i_d.z_max=-10


 var_list = var.split(',')

 # convert mean
 print 'Loading file: ', 'mean'
 u = i_d.mean_field.reshape((num_points,3),order='F')
 for j in range (1,len(var_list)):
	if var_list[j] == 'TE':
		tke = io_vtk.get_1_var(grid,var_list[j]+'_POD')
	if var_list[j] == 'SDR':
		epsi = io_vtk.get_1_var(grid,var_list[j]+'_POD') 
 cc = vtk.vtkCellCenters()
 cc.SetInputData(grid)
 cc.VertexCellsOn()
 cc.Update()
 points = VN.vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
                                          
 filename = rdir +'PODFS_mean.prf'

 target = open(filename,'w')

 target.write('# ' + 'PODFS_mean # name of the profile\n')
 target.write('# turbulence model, ' + i_d.turbulence_model+'\n')
 plane_rhs = (o1+i_d.t_o[0])*n1+(o2+i_d.t_o[1])*n2+(o3+i_d.t_o[2])*n3
 target.write('# plane normal and translation ' + str(n1)+'\t'+str(n2)+'\t'+str(n3)+'\t'+str(plane_rhs)+'\n')
                                        #target.write(str(src_output.GetNumberOfPoints())+'\n')
                                
 target.write('type, xyz # type of profile (rad or xyz)\n')
 target.write('localcs,origin,0,0,0 # origin of local coordinate system\n')
 target.write('localcs,xaxis,1,0,0 # x axis direction of local coordinate system\n')
 target.write('localcs,yaxis,0,1,0 # y axis direction of local coordinate system\n')
 target.write('localcs,zaxis,0,0,1 # z axis direction of local coordinate system\n')
 target.write('tolerance, 1.00E-08 # tolerance\n')
 if (i_d.turbulence_model == 'none'):
	target.write('scale,1,1,1,1,1,1 # scaling factors\n')
 else:
	target.write('scale,1,1,1,1,1,1,1,1 # scaling factors\n')
 if (i_d.turbulence_model == 'k_epsilon'):
                target.write('data,x,y,z,u,v,w,k,e\n')
 if (i_d.turbulence_model == 'k_omega'):
        	target.write('data,x,y,z,u,v,w,k,sdr\n')
 else:
        	target.write('data,x,y,z,u,v,w\n')
 if (len(var_list) > 1):
	for j in range (0,grid.GetNumberOfCells()):
        	target.write(sp.str(points[j,0])+','+sp.str(points[j,1])+','+sp.str(points[j,2])+','+sp.str(u[j,0])+','+sp.str(u[j,1])+','+sp.str(u[j,2]) + ',' + sp.str(tke[j])+','+sp.str(epsi[j])+'\n')
 else:
	for j in range (0,grid.GetNumberOfCells()):
        	target.write(sp.str(points[j,0])+','+sp.str(points[j,1])+','+sp.str(points[j,2])+','+sp.str(u[j,0])+','+sp.str(u[j,1])+','+sp.str(u[j,2]) + '\n')
	                       
                 
 target.close()

 # convert modes 

 ############################################
 # for modes chosen calculate energy, load spatial modes and calculate vars
 t_energy = 0
 for i in range (0,num_modes):
    #t_energy += energy[mode]
    #load vtk file
    ii = i+1
    counter = '%4.4i'% ii
    print 'Saving mode: ', counter
    u = i_d.spatial_modes[:,i].reshape((num_points,3),order='F') #
                                         
    filename = rdir +'PODFS_mode_'+counter+'.prf'

    target = open(filename,'w')

    target.write('# ' + 'PODFS_mode_'+counter + ' # name of the profile\n')
    target.write('# turbulence model, ' + i_d.turbulence_model+'\n')
    plane_rhs = (o1)*n1+(o2)*n2+(o3)*n3
    target.write('# plane normal and translation ' + str(n1)+'\t'+str(n2)+'\t'+str(n3)+'\t'+str(plane_rhs)+'\n')
                                        #target.write(str(src_output.GetNumberOfPoints())+'\n')
                                
    target.write('type, xyz # type of profile (rad or xyz)\n')
    target.write('localcs,origin,0,0,0 # origin of local coordinate system\n')
    target.write('localcs,xaxis,1,0,0 # x axis direction of local coordinate system\n')
    target.write('localcs,yaxis,0,1,0 # y axis direction of local coordinate system\n')
    target.write('localcs,zaxis,0,0,1 # z axis direction of local coordinate system\n')
    target.write('tolerance, 1.00E-08 # tolerance\n')
    if (i_d.turbulence_model == 'none'):
        target.write('scale,1,1,1,1,1,1 # scaling factors\n')
    else:
        target.write('scale,1,1,1,1,1,1,1,1 # scaling factors\n')
    if (i_d.turbulence_model == 'k_epsilon'):
                target.write('data,x,y,z,u,v,w,k,e\n')
    if (i_d.turbulence_model == 'k_omega'):
        	target.write('data,x,y,z,u,v,w,k,sdr\n')
    else:
	target.write('data,x,y,z,u,v,w\n')
    if (len(var_list) > 1):
    	for j in range (0,grid.GetNumberOfCells()):
        	target.write(sp.str(points[j,0])+','+sp.str(points[j,1])+','+sp.str(points[j,2])+','+sp.str(u[j,0])+','+sp.str(u[j,1])+','+sp.str(u[j,2]) + ',' + sp.str(tke[j])+','+sp.str(epsi[j])+'\n')
                                       
    else:                                   
        for j in range (0,grid.GetNumberOfCells()):
                target.write(sp.str(points[j,0])+','+sp.str(points[j,1])+','+sp.str(points[j,2])+','+sp.str(u[j,0])+','+sp.str(u[j,1])+','+sp.str(u[j,2]) + '\n')
 
    target.close()



#=======================================================================
# EOF
#=======================================================================
