# nicks post processing library

import numpy as np
import sys
import vtk.util.numpy_support as VN
#from mpi4py import MPI
#import io_vtk_obj_no_mpi as io_vtk


def fct_welch(x,fs,N,iwindow):

    
    # CHECKING BLOCK SIZE
    if(N>x.size):
        print '*** ERROR: Block size N should not be larger than the signal size.'
        print '*** STOP !!!'
        sys.exit()



    #%%%%% WINDOWING FUNCTION
    if(iwindow==1):
        w=np.array(np.ones((N), dtype=np.float64))
        print 'WINDOW: rectangular'
    elif(iwindow==2):
        w=np.hanning(N)
        print 'WINDOW: Hanning'
    elif(iwindow==3):
        w=np.blackman(N)
        print 'Window: Blackman'

    Cw=N/np.sum(w**2);   #%%% Corrective factor


    #%%%%% CALCULATION OF THE NUMBER OF BLOCKS (M)
    noverlap=int(np.floor(N/2)) #   % 50% overlap
    go_on=1
    n_end_block=N
    M=1
    while (go_on==1):
        n_end_block=n_end_block+noverlap
        if(n_end_block<=x.size):
            M=M+1
            go_on=1
        elif(n_end_block>x.size):
            go_on=0

    print 'Welch: number of blocks: M= ' + str(M)


    #%%%%% FREQUENCY VECTOR
    #f=np.array(np.zeros((1,N), dtpype=np.float64))
    f=np.linspace(-N/2,N/2-1,N)/N*fs

    #%%%%% CALCULATION OF THE PSD (Sxx)
    #%%%%% THIS REQUIRES CALCULATING THE PSD OF EACH BLOCK (Sxx_current_block) AND AVERAGING
    Sxx=np.array(np.zeros((N), dtype=np.complex64))
    Sxxsum=np.array(np.zeros((N), dtype=np.float64))
    for j in range (1,M+1):
        #print N, w.size, noverlap
        Sxx[:]=np.fft.fft(x[(j-1)*noverlap:(j-1)*noverlap+N]*w)
        Sxx[:]=np.fft.fftshift(Sxx)
        #disp_NB([size(Sxxsum),size(Sxx)])
        Sxxsum[:]=Sxxsum+Cw/N/fs*Sxx*np.conj(Sxx)

    Sxx[:]=Sxxsum/M

    return f,Sxx,M

#--------------------------------------------------------------------
def fct_iwelch(x,y,fs,N,iwindow):

    
    # CHECKING BLOCK SIZE
    if(N>x.size):
        print '*** ERROR: Block size N should not be larger than the signal size.'
        print '*** STOP !!!'
        sys.exit()



    #%%%%% WINDOWING FUNCTION
    if(iwindow==1):
        w=np.array(np.ones((N), dtype=np.float64))
        print 'WINDOW: rectangular'
    elif(iwindow==2):
        w=np.hanning(N)
        print 'WINDOW: Hanning'
    elif(iwindow==3):
        w=np.blackman(N)
        print 'Window: Blackman'

    Cw=N/np.sum(w**2);   #%%% Corrective factor


    #%%%%% CALCULATION OF THE NUMBER OF BLOCKS (M)
    noverlap=int(np.floor(N/2)) #   % 50% overlap
    go_on=1
    n_end_block=N
    M=1
    while (go_on==1):
        n_end_block=n_end_block+noverlap
        if(n_end_block<=x.size):
            M=M+1
            go_on=1
        elif(n_end_block>x.size):
            go_on=0

    print 'Welch: number of blocks: M= ' + str(M)


    #%%%%% FREQUENCY VECTOR
    #f=np.array(np.zeros((1,N), dtpype=np.float64))
    f=np.linspace(-N/2,N/2-1,N)/N*fs

    #%%%%% CALCULATION OF THE PSD (Sxx)
    #%%%%% THIS REQUIRES CALCULATING THE PSD OF EACH BLOCK (Sxx_current_block) AND AVERAGING
    Sxx=np.array(np.zeros((N), dtype=np.complex64))
    Syy=np.array(np.zeros((N), dtype=np.complex64))
    Sxxsum=np.array(np.zeros((N), dtype=np.complex64))
    for j in range (1,M+1):
        #print N, w.size, noverlap
        Sxx[:]=np.fft.fft(x[(j-1)*noverlap:(j-1)*noverlap+N]*w)
        Sxx[:]=np.fft.fftshift(Sxx)
        Syy[:]=np.fft.fft(y[(j-1)*noverlap:(j-1)*noverlap+N]*w)
        Syy[:]=np.fft.fftshift(Syy)
        #disp_NB([size(Sxxsum),size(Sxx)])
        Sxxsum[:]=Sxxsum+Cw/N/fs*Sxx*np.conj(Syy)

    Sxx[:]=Sxxsum/M

    return f,Sxx,M
#--------------------------------------------------------------------
def cross_correlation(x,y,fs,N,iwindow):

    
    # CHECKING BLOCK SIZE
    if(N>x.size):
        print '*** ERROR: Block size N should not be larger than the signal size.'
        print '*** STOP !!!'
        sys.exit()



    #%%%%% WINDOWING FUNCTION
    if(iwindow==1):
        w=np.array(np.ones((N), dtype=np.float64))
        print 'WINDOW: rectangular'
    elif(iwindow==2):
        w=np.hanning(N)
        print 'WINDOW: Hanning'
    elif(iwindow==3):
        w=np.blackman(N)
        print 'Window: Blackman'

    Cw=N/np.sum(w**2);   #%%% Corrective factor


    #%%%%% CALCULATION OF THE NUMBER OF BLOCKS (M)
    noverlap=int(np.floor(N/2)) #   % 50% overlap
    go_on=1
    n_end_block=N
    M=1
    while (go_on==1):
        n_end_block=n_end_block+noverlap
        if(n_end_block<=x.size):
            M=M+1
            go_on=1
        elif(n_end_block>x.size):
            go_on=0

    print 'Welch: number of blocks: M= ' + str(M)


    #%%%%% FREQUENCY VECTOR
    #f=np.array(np.zeros((1,N), dtpype=np.float64))
    f=np.linspace(-N/2,N/2-1,N)/fs

    #%%%%% CALCULATION OF THE PSD (Sxx)
    #%%%%% THIS REQUIRES CALCULATING THE PSD OF EACH BLOCK (Sxx_current_block) AND AVERAGING
    Sxx=np.array(np.zeros((N), dtype=np.complex64))
    Syy=np.array(np.zeros((N), dtype=np.complex64))
    Sxxsum=np.array(np.zeros((N), dtype=np.complex64))
    for j in range (1,M+1):
        #print N, w.size, noverlap
        Sxx[:]=np.fft.fft(x[(j-1)*noverlap:(j-1)*noverlap+N]*w)
        #Sxx[:]=np.fft.fftshift(Sxx)
        Syy[:]=np.fft.fft(y[(j-1)*noverlap:(j-1)*noverlap+N]*w)
        #Syy[:]=np.fft.fftshift(Syy)
        #disp_NB([size(Sxxsum),size(Sxx)])
        Sxxsum[:]=Sxxsum+Cw/N/fs*np.fft.fftshift(np.fft.ifft(Sxx*np.conj(Syy)))

    Sxx[:]=Sxxsum/M

    # normalise
    #Sxx[:]=Sxx[:]/np.amax(Sxx[:])

    return f,Sxx,M

#--------------------------------------------------------------------
def fct_transfer_func(x,y,fs,N,iwindow):

    f,Sxy,M=fct_iwelch(x,y,fs,N,iwindow)
    f,Sxx,M=fct_welch(x,fs,N,iwindow)

    H = Sxy/Sxx

    return f,H,M
#--------------------------------------------------------------------
def fct_coherence(x,y,fs,N,iwindow):

    f,Sxy,M=fct_iwelch(x,y,fs,N,iwindow)
    f,Sxx,M=fct_welch(x,fs,N,iwindow)
    f,Syy,M=fct_welch(y,fs,N,iwindow)

    H = np.abs(Sxy)**2/Sxx/Syy

    return f,H,M,Sxy

#======================================================
def mean(mat,dim=0):
	size1 = np.shape(mat)
#	print size1, np.size(size1)
        if (np.size(size1) == 2):
            if (dim == 1):
		size2 = size1[0]
		dim2 = 1
		mean = np.array(np.zeros((size2), dtype=np.float64))
		
		for i in range (0,size1[dim2]):
                	mean[:] = mean[:] + mat[:,i]
            else:
		size2 = size1[1]
		dim2 = 0
		mean = np.array(np.zeros((size2), dtype=np.float64))

		for i in range (0,size1[dim2]):
			mean[:] = mean[:] + mat[i,:]

        else:
            size2 = np.size(mat)
            dim2 = 1
            mean = 0.0
		
            for i in range (0,size2):
                	mean = mean + mat[i]
          
	mean = mean/float(size2)

	return mean

#========================================================
def azimuthal_fourier_series(mesh,var,nc,axiseg,rseg,aziseg):

    #define arrays
    u = np.zeros((mesh.GetNumberOfPoints(),nc+1),dtype=np.float64)
    m_array = np.zeros((nc+1,axiseg,rseg,aziseg),dtype=np.float64)
    nan_array = np.ones((nc+1,axiseg,rseg,aziseg),dtype=np.float64)
    x_array = np.zeros((axiseg,rseg),dtype=np.float64)
    r_array = np.zeros((axiseg,rseg),dtype=np.float64)


    #vtk to numpy
    if nc > 1:
        u[:,0:nc] = VN.vtk_to_numpy(mesh.GetPointData().GetVectors(var))
        for i in range (0,nc):
            u[:,nc] += u[:,i]**2
        u[:,nc] = np.sqrt(u[:,nc])
    else:
        u[:,0:nc] = VN.vtk_to_numpy(mesh.GetPointData().GetScalars(var))

    

    #print np.amax(u)
    r = VN.vtk_to_numpy(mesh.GetPointData().GetScalars('r'))
    theta = VN.vtk_to_numpy(mesh.GetPointData().GetScalars('theta'))
    x = VN.vtk_to_numpy(mesh.GetPointData().GetScalars('x'))
    
   
    #for each component of variable
    for i in range (0,nc+1): #for each component 
        index1 = 0
        for jj in range (0,axiseg):  #for each axial segment      
            for k in range(0,rseg): #for each radial segment
                index1 += 1 #move away from centre
                ttemp = np.zeros(aziseg,dtype=np.float64)
                utemp = np.zeros(aziseg,dtype=np.complex64)
                for n in range(0,aziseg): #for each azimuthal segment
                    for m in range(0,aziseg):
                        ttemp[m] = theta[index1]
                        #print i,j,k,n,m,r[index1]
                        #if jj == axiseg-1:
                        #    print u[index1,i],x[index1],r[index1]
                        utemp[m] = u[index1,i]*np.exp(-1j*n*ttemp[m])
                        x_array[jj,k] = x[index1]
                        r_array[jj,k] = r[index1]
                        #if u[index1,i] == 0.0:
                            #nan_array[i,jj,k,n] = np.nan
                        index1 += 1 #move around circle
                    index1 -= aziseg #move back to start point
                    #print np.amaxutemp
                    temp = np.trapz(utemp,ttemp) 
                   
                    m_array[i,jj,k,n] = temp*temp.conj()
                index1 += aziseg-1
              
            index1 += 1

    #m_array *= nan_array
    #print m_array[-1,0,2,:]
    return m_array,x_array,r_array


#========================================================
def MPI_mean_2D(B,shape,i_d):

    global_num_snapshots = i_d.global_num_snapshots
    comm = i_d.comm
    rank = i_d.rank
    size = i_d.size

    #print np.shape(B)

    mean_v = np.array(np.zeros((shape), dtype=np.float64))
    mean_v[:] = np.sum(B,1)

    if rank == 0:
        rbuff = np.array(np.zeros((size,shape), dtype=np.float64))
    else:
        rbuff = []
    comm.Gather(mean_v[:],rbuff, root =0)
    if rank == 0: #calculate and broadcast mean
        mean_v[:] = np.sum(rbuff,0)/global_num_snapshots
    comm.Bcast(mean_v,root=0)


    return mean_v

#========================================================
def MPI_mean_3D(B,shape1,shape2,i_d):

    global_num_snapshots = i_d.global_num_snapshots
    comm = i_d.comm
    rank = i_d.rank
    size = i_d.size

    mean_v = np.array(np.zeros((shape1,shape2), dtype=np.float64))
    mean_v[:,:] = np.sum(B,2)
    
    if rank == 0:
        rbuff = np.array(np.zeros((size,shape1,shape2), dtype=np.float64))
    else:
        rbuff = []
    comm.Gather(mean_v[:],rbuff, root =0)
    if rank == 0: #calculate and broadcast mean
        mean_v[:,:] = np.sum(rbuff,0)/global_num_snapshots
    comm.Bcast(mean_v,root=0)


    return mean_v


#========================================================
def MPI_local_to_global_2D(B,shape1,i_d):


    #print B, i_d.rank,'err'

    if B ==[]:
        return []

    global_num_snapshots = i_d.global_num_snapshots
    comm = i_d.comm
    rank = i_d.rank
    size = i_d.size

   
   
    
    if rank == 0:
        C = np.array(np.zeros((shape1,global_num_snapshots), dtype=np.float64))
        #print size
        for i in range (1,size):
            buffer1 = np.array(np.zeros((shape1,i_d.npp), dtype=np.float64))
            comm.Recv(buffer1,source = i, tag =1)
            #print i,buffer1
            C[:,i_d.npp*(i-1):i_d.npp*i] = buffer1
        C[:,i_d.npp*(size-1):global_num_snapshots] = B
    else:
        C= []
        comm.Send(B,dest = 0,tag =1)
        del B

    return C

#========================================================
def MPI_local_to_global_3D(B,shape1,shape2,i_d):

    if B ==[]:
        return []

    global_num_snapshots = i_d.global_num_snapshots
    comm = i_d.comm
    rank = i_d.rank
    size = i_d.size

    #print np.shape(B), rank
   
    
    if rank == 0:
        C = np.array(np.zeros((shape1,shape2,global_num_snapshots), dtype=np.float64))
        for i in range (1,size):
            buffer1 = np.array(np.zeros((shape1,shape2,i_d.npp), dtype=np.float64))
            comm.Recv(buffer1,source = i, tag =1)
            #print buffer1
            C[:,:,i_d.npp*(i-1):i_d.npp*i] = buffer1
        C[:,:,i_d.npp*(size-1):global_num_snapshots] = B
    else:
        C= []
        comm.Send(B,dest = 0,tag =1)
        del B

    return C
#========================================================
def MPI_local_to_global_1D(B,i_d):

    if B ==[]:
        return []

    global_num_snapshots = i_d.global_num_snapshots
    comm = i_d.comm
    rank = i_d.rank
    size = i_d.size

    #print np.shape(B), rank
   
    
    if rank == 0:
        C = np.array(np.zeros((global_num_snapshots), dtype=np.float64))
        
        for i in range (1,size):
            buffer1 = np.array(np.zeros((i_d.npp), dtype=np.float64))
            comm.Recv(buffer1,source = i, tag =1)
            C[i_d.npp*(i-1):i_d.npp*i] = buffer1
        C[i_d.npp*(size-1):global_num_snapshots] = B
    else:
        C= []
        comm.Send(B,dest = 0,tag =1)
        del B

    return C

#========================================================
def MPI_local_to_global_1D_npp(B,i_d,npp,global_num_snapshots):

    if B ==[]:
        return []

    comm = i_d.comm
    rank = i_d.rank
    size = i_d.size

    #print npp, rank
    #print global_num_snapshots
    
    if rank == 0:
        C = np.array(np.zeros(global_num_snapshots, dtype=np.float64))
        
        for i in range (1,size):
            buffer1 = np.array(np.zeros((npp), dtype=np.float64))
            comm.Recv(buffer1,source = i, tag =1)
            a = npp*(i-1)
            b = npp*i
            print i,a,b, np.shape(buffer1)
            C[a:b] = buffer1
        C[npp*(size-1):global_num_snapshots] = B
    else:
        C= []
        comm.Send(B,dest = 0,tag =1)
        del B

    return C

#===================================================================
def DivideByNP(startig,endig,i_d):
    

    npp = int(np.floor((endig-startig)/i_d.size))
    #print startig,endig,npp
    
    csi = startig
    cei = csi + npp
    for i in range (1,i_d.size):
        if i_d.rank == i:
            starti = csi
            endi = cei
        csi = cei 
        cei =  csi +npp

    if i_d.rank == 0:
        starti = csi
        endi = endig

    return starti,endi,npp

#======================================================================
def phase_average_MPI(B,scalars_array,POD_plane_velocity,i_d):


    print '\n   Phase Averaging ...'
    num_phase = int(1./i_d.phase_freq/i_d.del_t/i_d.stride)
    if i_d.rank != 0:
    	time_vec = i_d.del_t*i_d.stride*np.linspace((i_d.rank-1)*i_d.npp,(i_d.rank-1)*i_d.npp+i_d.num_snapshots,i_d.num_snapshots)
    else:
	time_vec = i_d.del_t*i_d.stride*np.linspace(0,i_d.global_num_snapshots,i_d.global_num_snapshots+1)
    del_t1 = i_d.stride*i_d.del_t
    period = 1./i_d.phase_freq
    #print num_phase,del_t1,period
    if (num_phase >= i_d.global_num_snapshots):
        num_phase = i_d.global_num_snapshots -1
        if i_d.rank == 0:
            time_vec[-1] = time_vec[i_d.num_snapshots]

    i_d.num_phase=num_phase
    dt_phase = 1./i_d.phase_freq/num_phase

    print '    Number of snapshots per period: ', num_phase
    scalars_phase=np.zeros((i_d.num_scalars,i_d.num_points,num_phase), dtype=np.float64)
    v_phase=np.zeros((i_d.num_points*i_d.num_components,num_phase), dtype=np.float64)
    nrs_phase = np.array(np.zeros((i_d.num_points*i_d.num_components,num_phase), dtype=np.float64))
    crs_phase = np.array(np.zeros((i_d.num_points*i_d.num_components,num_phase), dtype=np.float64))
    phase_factor = np.zeros((num_phase),dtype=np.int32)
    #if (turbulence_model == 'k_epsilon'):
    #    k_phase=np.zeros((num_points,num_phase), dtype=np.float64)
    #    epsi_phase=np.zeros((num_points,num_phase), dtype=np.float64)
    #elif (turbulence_model == 'k_omega'):
    #    k_phase=np.zeros((num_points,num_phase), dtype=np.float64)
    #    epsi_phase=np.zeros((num_points,num_phase), dtype=np.float64)


    # set up windows
    Ainbuff = np.zeros((i_d.num_scalars,i_d.num_points),dtype = np.float64)
    Aoutbuff = []
    Awin = []
    counter1 = 0
    for i in range (0,i_d.npp):
		Aoutbuff.append(np.zeros((i_d.num_scalars,i_d.num_points),dtype = np.float64))
		#for j in range (0,i_d.num_scalars):
                    #print scalars_array
                    #Aoutbuff[i][j*i_d.num_points:(j+1)*i_d.num_points]=0
                if i < i_d.num_snapshots :
                    Aoutbuff[i][:,:] = scalars_array[:,:,i]
		Awin.append(MPI.Win.Create(Aoutbuff[i],comm=i_d.comm))
                #print 'Awin set for', i_d.rank
    #print Awin

    Binbuff = np.zeros((i_d.num_points*i_d.num_components),dtype = np.float64)
    Boutbuff = []
    Bwin = []
    for i in range (0,i_d.npp):
                Boutbuff.append(np.zeros(i_d.num_points*i_d.num_components,dtype = np.float64))
                if i < i_d.num_snapshots :
                    Boutbuff[i][:] = B[:,i]
                Bwin.append(MPI.Win.Create(Boutbuff[i],comm=i_d.comm))
                #print 'Bwin set for', i_d.rank



    if i_d.rank == 0:
      #perform phase average
      counter1 = 0
      counter2 = 0
      time = 0
      time_loc = 0
      while (time <= time_vec[-1]):
        print 'Time = ' , time, ' of ', time_vec[-1]
	#for i in range (0,i_d.num_scalars):
        scalars_phase[:,:,counter2] =  scalars_phase[:,:,counter2] + \
                    temporal_interpolation_MPI(scalars_array[:,:,:],Awin, \
                    Ainbuff,time_vec,time,del_t1,i_d)

        v_phase[:,counter2] =  v_phase[:,counter2] + \
		temporal_interpolation_MPI \
		(B,Bwin,Binbuff,time_vec,time,del_t1,i_d)
        nrs_phase[:,counter2] = nrs_phase[:,counter2]+ \
		temporal_interpolation_mult_MPI \
		(B,Bwin,Binbuff,0,i_d.num_points*i_d.num_components, \
		0,i_d.num_points*i_d.num_components,time_vec,time,del_t1,i_d)
        crs_phase[0:i_d.num_points,counter2] = crs_phase[0:i_d.num_points,counter2]+ \
		temporal_interpolation_mult_MPI \
		(B,Bwin,Binbuff,0,i_d.num_points, \
		i_d.num_points,2*i_d.num_points,time_vec,time,del_t1,i_d)
	crs_phase[i_d.num_points:2*i_d.num_points,counter2]  = \
		crs_phase[i_d.num_points:2*i_d.num_points,counter2] + \
		temporal_interpolation_mult_MPI \
		(B,Bwin,Binbuff,2*i_d.num_points,3*i_d.num_points, \
		i_d.num_points,2*i_d.num_points,time_vec,time,del_t1,i_d)
	crs_phase[2*i_d.num_points:3*i_d.num_points,counter2] = \
		crs_phase[2*i_d.num_points:3*i_d.num_points,counter2] + \
		temporal_interpolation_mult_MPI \
		(B,Bwin,Binbuff,0,i_d.num_points,\
		2*i_d.num_points,3*i_d.num_points,time_vec,time,del_t1,i_d)

        counter2 += 1
        counter1 += 1
        phase_factor[counter2-1] += 1
        time += dt_phase
	time_loc += dt_phase
        #print phase_factor[counter2-1], counter2
        if (counter2 == num_phase):
            counter2 = 0


      scalars_phase = scalars_phase/phase_factor
      v_phase = v_phase/phase_factor
      nrs_phase = nrs_phase/phase_factor
      crs_phase = crs_phase/phase_factor


    phase = np.linspace(0,360,num_phase)


    # broadcast phase averages
    i_d.comm.Bcast(scalars_phase, root = 0)
    i_d.comm.Bcast(v_phase, root = 0)
    i_d.comm.Bcast(nrs_phase, root = 0)
    i_d.comm.Bcast(crs_phase, root = 0)


    print '\n   Subtracting phase averaged quantities from instantaneous fields ...'

    # find double-dash quantities
    if i_d.rank == 0: 
        time_vec = i_d.del_t*i_d.stride*np.linspace((i_d.size-1)*i_d.npp,i_d.global_num_snapshots,i_d.num_snapshots)
    counter2 = 0
    time = 0
    i_loc = 0
    while (time < time_vec[-1]):
        if time >= time_vec[i_loc]:
            #if time >= time_vec[i_loc+1]:
            #    i_loc += 1
            i1 = i_loc
            i2 = counter2+1
            if i2 == num_phase:
                i2 = 0
            #i_loc += 1

            # find weights
            w1 = (time_vec[i_loc+1]-time)/i_d.del_t
            if w1 < 0:
                print 'error in phase average'
            w2 = 1.-w1

            #subtract phase
            #if i_d.rank == 0:
            #print i1,i2,counter2
            B[:,i1] = B[:,i1] - 0.5*(w1*v_phase[:,counter2] + w2*v_phase[:,i2])
            scalars_array[:,:,i1] = scalars_array[:,:,i1] - 0.5*(w1*scalars_phase[:,:,counter2]+w2*scalars_phase[:,:,i2])
            t2  = time + dt_phase
            if t2 > time_vec[i_loc]: # posibility of dt<dt_phase
                i_loc += 1
            if t2 > time_vec[i_loc]: # possibility of dt>dt_phase
                i_loc += 1
        counter2 += 1
        time += dt_phase
        if (counter2 == num_phase):
            counter2 = 0

    print '\n   Calculating non-periodic Resolved Reynolds Stresses and RMS of pressure...'
    rms_scalars = np.array(np.zeros((i_d.num_scalars,i_d.num_points), dtype=np.float64))
    rms_scalars[:,:] = np.sqrt(np.mean(scalars_array[:,:,:]*scalars_array[:,:,:],2))

    nrs = np.array(np.zeros((i_d.num_points*i_d.num_components), dtype=np.float64))
    crs = np.array(np.zeros((i_d.num_points*i_d.num_components), dtype=np.float64))

    nrs[:] = MPI_mean_2D(B[:,:]*B[:,:],i_d.num_points*i_d.num_components,i_d)
    crs[0:i_d.num_points] = MPI_mean_2D(B[0:i_d.num_points,:]*B[i_d.num_points:2*i_d.num_points,:],i_d.num_points,i_d)
    crs[i_d.num_points:2*i_d.num_points] = MPI_mean_2D(B[2*i_d.num_points:3*i_d.num_points,:]*B[i_d.num_points:2*i_d.num_points,:],i_d.num_points,i_d)
    crs[2*i_d.num_points:3*i_d.num_points] = MPI_mean_2D(B[0:i_d.num_points,:]*B[2*i_d.num_points:3*i_d.num_points,:],i_d.num_points,i_d)


    for i in range (0,i_d.npp):
                Awin[i].Fence(0)
                Bwin[i].Fence(0)
                Awin[i].Free()
                Bwin[i].Free()

    return scalars_phase,v_phase,nrs_phase,crs_phase,rms_scalars,nrs,crs

#=======================================================================
def temporal_interpolation_MPI(field,Awin,Ain,time_vec,time,dt,i_d):

        # find time index
        t=0
        i=0
	i_tot = 0
	act_proc = 1
        while t < time:
                i=i+1
                i_tot+=1
                t=t+dt
		if i > i_d.npp-1:
			i = 0
			act_proc += 1
		if act_proc > i_d.size-1:
			act_proc = 0
			

        if (i_tot > (time_vec.size -2)):
                i = i - 1
                #print 'act 1',i,act_proc
        else:
                test = time_vec[i_tot]-time - dt
                #print time_vec[i_tot]-time_vec[i_tot-1],dt
                if test > 1.e-16 :
                    
                    t=t-dt
                    i=i-1
                    i_tot -= 1
                    if i < 0 and act_proc == 0:
                            i = i_d.npp-1
                            #i_tot -= 1
                            act_proc = i_d.size -1
                    if i < 0 and act_proc != 0:
                            i = i_d.npp-1
                            act_proc -= 1
                            #i_tot -= 1
                    #print 'act 2', i,act_proc,time_vec[i_tot+1]-time - dt


        #print i_tot,time_vec.size

        # find weights
        w1 = (time_vec[i_tot]-time)/dt
        #print w1
        w2 = 1.-w1
        shape = np.shape(Ain)
        size = np.size(Ain)
        if size == shape[0]:
            interp_field = np.zeros(size,dtype = np.float64)
            dims = 1
        else:
            interp_field = np.zeros((shape[0],shape[1]),dtype = np.float64)
            dims = 2
	if act_proc != 0:
		#Awin[i].Lock(MPI.LOCK_SHARED,act_proc,0)
                Awin[i].Lock(act_proc,MPI.LOCK_SHARED,0)
                #Awin[i].Get(Ain,act_proc)
		interp_field[:] = w2*Ain
		Awin[i].Unlock(act_proc)
		if i == 0:
			act_proc = act_proc - 1
			i = i_d.npp-1
			#Awin[i].Lock(MPI.LOCK_SHARED,act_proc,0)
                	Awin[i].Lock(act_proc,MPI.LOCK_SHARED,0)
                	Awin[i].Get(Ain,act_proc)
                	interp_field[:] += w1*Ain
                	Awin[i].Unlock(act_proc)
		else:
			i -= 1
                        #Awin[i].Lock(MPI.LOCK_SHARED,act_proc,0)
			Awin[i].Lock(act_proc,MPI.LOCK_SHARED,0)
                        Awin[i].Get(Ain,act_proc)
                        interp_field[:] += w1*Ain
                        Awin[i].Unlock(act_proc)
	else:
                if dims == 1 :
                    interp_field[:] = w2*field[:,i]
                else:
                    interp_field[:] = w2*field[:,:,i]
		if i == 0:
			act_proc = i_d.size-1
                        i = i_d.npp -1
                        #print i
                        #Awin[i].Lock(MPI.LOCK_SHARED,act_proc,0)
                        Awin[i].Lock(act_proc,MPI.LOCK_SHARED,0)
                        Awin[i].Get(Ain,act_proc)
                        interp_field[:] += w1*Ain
                        Awin[i].Unlock(act_proc)
		else:
			i -= 1
                        if dims == 1 :
                            interp_field[:] += w1*field[:,i]
                        else:
                            interp_field[:] += w1*field[:,:,i]

        return interp_field


#=======================================================================
def temporal_interpolation_mult_MPI(field,Awin,Ain,start1,end1,start2,end2,time_vec,time,dt,i_d):

    # find time index
        t=0
        i=0
	i_tot = 0
	act_proc = 1
        while t < time:
                i=i+1
                i_tot+=1
                t=t+dt
		if i > i_d.npp-1:
			i = 0
			act_proc += 1
		if act_proc > i_d.size-1:
			act_proc = 0
			

        if (i_tot > (time_vec.size -2)):
                i = i - 1
                #print 'act 1',i,act_proc
        else:
                test = time_vec[i_tot]-time - dt
                #print time_vec[i_tot]-time_vec[i_tot-1],dt
                if test > 1.e-16 :
                    
                    t=t-dt
                    i=i-1
                    i_tot -= 1
                    if i < 0 and act_proc == 0:
                            i = i_d.npp-1
                            #i_tot -= 1
                            act_proc = i_d.size -1
                    if i < 0 and act_proc != 0:
                            i = i_d.npp-1
                            act_proc -= 1
                            #i_tot -= 1

        # find weights
        w1 = (time_vec[i]-time)/dt
        w2 = 1.-w1
        
        interp_field1 = np.zeros(end1-start1,dtype = np.float64)
        interp_field2 = np.zeros(end2-start2,dtype = np.float64)
        if act_proc != 0:
                Awin[i].Lock(act_proc,MPI.LOCK_SHARED,0)
                #print i
                #Awin[i].Lock(MPI.LOCK_SHARED,act_proc,0)
                Awin[i].Get(Ain,act_proc)
                interp_field1[:] = w2*Ain[start1:end1]
                interp_field2[:] = w2*Ain[start2:end2]
                Awin[i].Unlock(act_proc)
                if i == 0:
                        act_proc = act_proc - 1
                        i = i_d.npp-1
                        Awin[i].Lock(act_proc,MPI.LOCK_SHARED,0)
                        #Awin[i].Lock(MPI.LOCK_SHARED,act_proc,0)
                        Awin[i].Get(Ain,act_proc)
                        interp_field1[:] += w1*Ain[start1:end1]
                        interp_field2[:] += w1*Ain[start2:end2]
                        Awin[i].Unlock(act_proc)
                else:
                        i -= 1
                        Awin[i].Lock(act_proc,MPI.LOCK_SHARED,0)
                        #Awin[i].Lock(MPI.LOCK_SHARED,act_proc,0)
                        Awin[i].Get(Ain,act_proc)
                        interp_field1[:] += w1*Ain[start1:end1]
                        interp_field2[:] += w1*Ain[start2:end2]
                        Awin[i].Unlock(act_proc)
        else:
        	interp_field1[:] = w2*field[start1:end1,i]
        	interp_field2[:] = w2*field[start2:end2,i]
                if i == 0:
                        act_proc = i_d.size-1
                        i = i_d.npp-1
                        Awin[i].Lock(act_proc,MPI.LOCK_SHARED,0)
                        #Awin[i].Lock(MPI.LOCK_SHARED,act_proc,0)
                        Awin[i].Get(Ain,act_proc)
                        interp_field1[:] += w1*Ain[start1:end1]
                        interp_field2[:] += w1*Ain[start2:end2]
                        Awin[i].Unlock(act_proc)
                else:
                        i -= 1
                        interp_field1[:] += w1*field[start1:end1,i]
                        interp_field2[:] += w1*field[start2:end2,i]

        return interp_field1*interp_field2

#=====================================================================
def str(val):
    
    return "%0.12f" %val


