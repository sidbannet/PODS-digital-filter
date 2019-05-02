
import sys
import string
import os
import re
import numpy as np
import matplotlib
#matplotlib.use('Qt4Agg')
from scipy import interpolate
from scipy.interpolate import griddata
from scipy.interpolate import Rbf
from scipy import spatial 
import matplotlib.mlab as mlab
from matplotlib import cm
from optparse import OptionParser
from optparse import Option, OptionValueError
from io import BytesIO
import PODFS as pod
import math
import nplotlib as plt


PROG = 'DigitalFilters'
VERSION = '1.0.0'

Pi = 3.1415927

class obj(object):
	a=0

class MultipleOption(Option):
    ACTIONS = Option.ACTIONS + ("extend",)
    STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
    TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
    ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

    def take_action(self, action, dest, opt, value, values, parser):
        if action == "extend":
            values.ensure_value(dest, []).append(value)
        else:
            Option.take_action(self, action, dest, opt, value, values, parser)

def coeff3D(a,nfx,nfy,nfz,lnx,lny,lnz):
   
   ax = np.zeros(nfx*2+1)
   ay = np.zeros(nfy*2+1)
   az = np.zeros(nfz*2+1)    

   ax=calccoeff(ax,nfx,lnx)   
   ay=calccoeff(ay,nfy,lny)
   az=calccoeff(az,nfz,lnz)      


   for i in range(nfx*2+1):
      for j in range(nfy*2+1):
         for k in range(nfz*2+1):
            a[0,i,j,k] = ax[i]*ay[j]*az[k]
	    
   '''
   norm = 0.0
   for k in range(nfz*2+1):
      for j in range(nfy*2+1):
          for i in range(nfx*2+1):
	     norm=norm+a[0,i,j,k]  	     
	     	     
   a[:,:,:,:] = a[:,:,:,:] / norm 	     
   '''
 

def calccoeff(a,n,ln):

   norm = 0.0
    
   for i in range(n*2+1):
      k = float(i-n)
       
      a[i] = np.exp(-Pi*k*k/(2.0*ln*ln))
      norm = norm+a[i]**2
   
   print np.amax(a)

   norm = np.sqrt(norm)
   b = a/norm

   print np.amax(b)
   return b

def filter3D(x,y,a,jma,kma,nfx,nfy,nfz):

   #  Filtering of the data. Input :x Output y                        *
   #  Filter coefficients: a  Filter length: 2*( nx,ny,nz )+1
   
   for k in range(kma):
      for j in range(jma):
          y[j,k] = np.sum(x[:,j:j+2*nfy+1,k:k+2*nfz+1]*a[0,:,:,:]) #y_filtered
          

def adapt1d(yu,yv,yw,uin,uuin,vvin,wwin,uwin,jma,kma):		   
   
   R = np.zeros((3,3))
   A = np.zeros((3,3))
   

   for k in range(kma):

      R[0,0] = uuin[k]
      R[1,1] = vvin[k]
      R[2,2] = wwin[k]
      #R[0,2] = uwin[k]
      R[2,0] = uwin[k] # changed from negative (doesn't matter) 
      
      A[:,:] = 1.0
      A[0,0] = np.sqrt(R[0,0])	 
      A[0,1] = 0.0 	 
      A[0,2] = 0.0
	      
      A[1,0] = R[1,0]/(A[0,0]+1e-20)      
      #print k, R[1,1],A[1,0],vvin[k],kma
      A[1,1] = np.sqrt(R[1,1]-A[1,0]*A[1,0])      
      A[1,2] = 0.0
      
      A[2,0] = R[2,0]/(A[0,0]+1e-20)      
      A[2,1] = (R[2,1]-A[1,0]*A[2,0])/(A[1,1]+1e-20)      
      A[2,2] = np.sqrt(R[2,2]-A[2,0]*A[2,0]-A[2,1]*A[2,1])      
	      
      for j in range(jma):
         xu = yu[j,k]     
	 xv = yv[j,k]
	 xw = yw[j,k]
	
	 yu[j,k] =  A[0,0]*xu + A[0,1]*xv + A[0,2]*xw + uin[k]
	 yv[j,k] =  A[1,0]*xu + A[1,1]*xv + A[1,2]*xw
	 yw[j,k] =  A[2,0]*xu + A[2,1]*xv + A[2,2]*xw

def adapt2prf(yu,yv,yw,uin,vin,win,uuin,vvin,wwin,uvin,uwin,vwin,jma,kma):		   
   
   R = np.zeros((3,3))
   A = np.zeros((3,3))
   
   #print jma,kma,np.shape(yu),np.shape(uin)
   for k in range(kma):
      for j in range(jma):

        R[0,0] = uuin[j,k]
        R[1,1] = vvin[j,k]
        R[2,2] = wwin[j,k]
        R[0,1] = uvin[j,k]
        R[1,0] = uvin[j,k]  # negative?
        R[0,2] = uwin[j,k]
        R[2,0] = uwin[j,k] # negative?
        R[1,2] = vwin[j,k]
        R[2,1] = vwin[j,k] # negative?
	 
        A[:,:] = 1.0
        A[0,0] = np.sqrt(R[0,0])	 
        A[0,1] = 0.0 	 
        A[0,2] = 0.0
	if A[0,0] > 0.:      
        	A[1,0] = R[1,0]/(A[0,0]+1e-20)
	else:
		A[1,0] = 0.0      
        if A[1,0]**2>R[1,1]:
		A[1,1]=0.0
	else:
        	A[1,1] = np.sqrt(R[1,1]-A[1,0]*A[1,0])      
        A[1,2] = 0.0
        if A[0,0] > 0.0:
        	A[2,0] = R[2,0]/(A[0,0]+1e-20)
	else:
		A[2,0] = 0.0
	if A[1,1]> 0.0:      
        	A[2,1] = (R[2,1]-A[1,0]*A[2,0])/(A[1,1]+1e-20)
	else:
		A[2,1] = 0.0
	if R[2,2]<A[2,0]*A[2,0]+A[2,1]*A[2,1]:
		A[2,2]=0.0
	else:      
        	A[2,2] = np.sqrt(R[2,2]-A[2,0]*A[2,0]-A[2,1]*A[2,1])      
	      
        xu = yu[j,k]     
        xv = yv[j,k]
        xw = yw[j,k]
	
        yu[j,k] =  A[0,0]*xu + A[0,1]*xv + A[0,2]*xw + uin[j,k]
        yv[j,k] =  A[1,0]*xu + A[1,1]*xv + A[1,2]*xw + vin[j,k]
        yw[j,k] =  A[2,0]*xu + A[2,1]*xv + A[2,2]*xw + win[j,k]

def adapt2d(yu,yv,yw,uin,uuin,vvin,wwin,uwin,jma,kma,mean_profile,inner_d):

  R = np.zeros((3,3))
  A = np.zeros((3,3))

  if mean_profile == 'double-hyperbolic-tangent':

   zArray = np.linspace(-1.,1,kma)
   zi = np.linspace(-1.,1,jma)
   uinj  = interpolate.splev(zi, interpolate.splrep(zArray, uin, s=0), der=0)
   uuinj = interpolate.splev(zi, interpolate.splrep(zArray, uuin, s=0), der=0)
   vvinj = interpolate.splev(zi, interpolate.splrep(zArray, vvin, s=0), der=0)
   wwinj = interpolate.splev(zi, interpolate.splrep(zArray, wwin, s=0), der=0)
   uwinj = interpolate.splev(zi, interpolate.splrep(zArray, uwin, s=0), der=0)

   uinj[0] = uin[0]	
   uuinj[0] = uuin[0]	
   vvinj[0] = vvin[0]	
   wwinj[0] = wwin[0]	
   uwinj[0] = uwin[0]	

   uinj[-1] = uin[-1] 
   uuinj[-1] = uuin[-1]
   vvinj[-1] = vvin[-1]
   wwinj[-1] = wwin[-1]
   uwinj[-1] = uwin[-1]

   for j in range(jma):
	
      # check for negative values in interpolation
      if uuinj[j] < 0.:
		uuinj[j] = 0.0
      if vvinj[j] < 0.:
		vvinj[j] = 0.0
      if wwinj[j] < 0.:
		wwinj[j] = 0.0

      for k in range(kma):
      
	R[0,0] = np.sqrt(uuin[k]*uuinj[j])#/2.
      	R[1,1] = np.sqrt(vvin[k]*vvinj[j])#/2.
      	R[2,2] = np.sqrt(wwin[k]*wwinj[j])#/2.
      	R[0,2] = np.sign(uwin[k]+uwinj[j])*np.sqrt(abs(uwin[k]*uwinj[j]))#/2.
      	R[2,0] = np.sign(uwin[k]+uwinj[j])*np.sqrt(abs(uwin[k]*uwinj[j]))#/2. #changed from minus

      	A[:,:] = 1.0
	temp1 = R[0,0]
	if temp1 <0:
		temp1 = 0
      	A[0,0] = np.sqrt(R[0,0])
      	A[0,1] = 0.0
      	A[0,2] = 0.0

      	A[1,0] = R[1,0]/(A[0,0]+1e-20)
      	#print k, R[1,1],A[1,0],vvin[k],kma
	temp1 = R[1,1]-A[1,0]*A[1,0]
	if temp1<0:
		temp1 = 0.0
      	A[1,1] = np.sqrt(temp1)
      	A[1,2] = 0.0

      	A[2,0] = R[2,0]/(A[0,0]+1e-20)
      	A[2,1] = (R[2,1]-A[1,0]*A[2,0])/(A[1,1]+1e-20)
	temp1 = R[2,2]-A[2,0]*A[2,0]-A[2,1]*A[2,1]
        if temp1<0:
                temp1 = 0.0
      	A[2,2] = np.sqrt(temp1)

        xu = yu[j,k]
        xv = yv[j,k]
        xw = yw[j,k]

        yu[j,k] =  A[0,0]*xu + A[0,1]*xv + A[0,2]*xw + np.sqrt(uin[k]*uinj[j])#/2.
        yv[j,k] =  A[1,0]*xu + A[1,1]*xv + A[1,2]*xw
        yw[j,k] =  A[2,0]*xu + A[2,1]*xv + A[2,2]*xw

  elif mean_profile == 'circular-hyperbolic-tangent':

    x = np.linspace(-1.,1.,jma)
    y = np.linspace(-1.,1.,kma)

    # find centre value index
    ci = np.argmax(uin)

    # re-interpolate to r
    zArray = np.linspace(0,1,len(uin)-ci)
    uinj  = interpolate.splrep(zArray, uin[ci:], s=0)
    uuinj = interpolate.splrep(zArray, uuin[ci:], s=0)
    vvinj = interpolate.splrep(zArray, vvin[ci:], s=0)
    wwinj = interpolate.splrep(zArray, wwin[ci:], s=0)
    uwinj = interpolate.splrep(zArray, uwin[ci:], s=0)

    for j in range (jma):
    
      for k in range(kma):
	
    	r = np.sqrt(x[j]**2+y[k]**2)
        uinr  = interpolate.splev(r, uinj, der=0)
        uuinr = interpolate.splev(r, uuinj, der=0)
        vvinr = interpolate.splev(r, vvinj, der=0)
        wwinr = interpolate.splev(r, wwinj, der=0)
        uwinr = interpolate.splev(r, uwinj, der=0)

	# reset boundaries to avoid dodgy values
	if r==0.0:
	 uinr = uin[ci]
	 uuinr = uuin[ci]
	 vvinr = vvin[ci]
	 wwinr = wwin[ci]
	 uwinr = uwin[ci]

        if r==1.0:
         uinr = uin[-1]
         uuinr = uuin[-1]
         vvinr = vvin[-1]
         wwinr = wwin[-1]
         uwinr = uwin[-1]

	# if r>0 then set to zero!
        if r>1.0:
         uinr = 0.0
         uuinr = 0.0
         vvinr = 0.0
         wwinr = 0.0
         uwinr = 0.0

        R[0,0] = uuinr
        R[1,1] = vvinr
        R[2,2] = wwinr
        R[0,2] = uwinr
        R[2,0] = uwinr #changed from minus

        A[:,:] = 1.0
        temp1 = R[0,0]
        if temp1 <0:
                temp1 = 0
        A[0,0] = np.sqrt(R[0,0])
        A[0,1] = 0.0
        A[0,2] = 0.0

        A[1,0] = R[1,0]/(A[0,0]+1e-20)
        temp1 = R[1,1]-A[1,0]*A[1,0]
        if temp1<0:
                temp1 = 0.0
        A[1,1] = np.sqrt(temp1)
        A[1,2] = 0.0

        A[2,0] = R[2,0]/(A[0,0]+1e-20)
        A[2,1] = (R[2,1]-A[1,0]*A[2,0])/(A[1,1]+1e-20)
        temp1 = R[2,2]-A[2,0]*A[2,0]-A[2,1]*A[2,1]
        if temp1<0:
                temp1 = 0.0
        A[2,2] = np.sqrt(temp1)

        xu = yu[j,k]
        xv = yv[j,k]
        xw = yw[j,k]

        yu[j,k] =  A[0,0]*xu + A[0,1]*xv + A[0,2]*xw + uinr
        yv[j,k] =  A[1,0]*xu + A[1,1]*xv + A[1,2]*xw
        yw[j,k] =  A[2,0]*xu + A[2,1]*xv + A[2,2]*xw

  elif mean_profile == 'ring-hyperbolic-tangent':

    x = np.linspace(-1.,1.,jma)
    y = np.linspace(-1.,1.,kma)

    # re-interpolate to r
    zArray = np.linspace(inner_d,1.,kma)
    uinj  = interpolate.splrep(zArray, uin, s=0)
    uuinj = interpolate.splrep(zArray, uuin, s=0)
    vvinj = interpolate.splrep(zArray, vvin, s=0)
    wwinj = interpolate.splrep(zArray, wwin, s=0)
    uwinj = interpolate.splrep(zArray, uwin, s=0)

    for j in range (jma):
    
      for k in range(kma):
	
    	r = np.sqrt(x[j]**2+y[k]**2)
        uinr  = interpolate.splev(r, uinj, der=0)
        uuinr = interpolate.splev(r, uuinj, der=0)
        vvinr = interpolate.splev(r, vvinj, der=0)
        wwinr = interpolate.splev(r, wwinj, der=0)
        uwinr = interpolate.splev(r, uwinj, der=0)

	# reset boundaries to avoid dodgy values
	if r==inner_d:
	 uinr = uin[0]
	 uuinr = uuin[0]
	 vvinr = vvin[0]
	 wwinr = wwin[0]
	 uwinr = uwin[0]

        if r==1.0:
         uinr = uin[-1]
         uuinr = uuin[-1]
         vvinr = vvin[-1]
         wwinr = wwin[-1]
         uwinr = uwin[-1]

	# if r>0 then set to zero!
        if r>1.0:
         uinr = 0.0
         uuinr = 0.0
         vvinr = 0.0
         wwinr = 0.0
         uwinr = 0.0

	# if r<inner_d then set to zero!
        if r<inner_d:
         uinr = 0.0
         uuinr = 0.0
         vvinr = 0.0
         wwinr = 0.0
         uwinr = 0.0


        R[0,0] = uuinr
        R[1,1] = vvinr
        R[2,2] = wwinr
        R[0,2] = uwinr
        R[2,0] = uwinr #changed from minus

        A[:,:] = 1.0
        temp1 = R[0,0]
        if temp1 <0:
                temp1 = 0
        A[0,0] = np.sqrt(R[0,0])
        A[0,1] = 0.0
        A[0,2] = 0.0

        A[1,0] = R[1,0]/(A[0,0]+1e-20)
        temp1 = R[1,1]-A[1,0]*A[1,0]
        if temp1<0:
                temp1 = 0.0
        A[1,1] = np.sqrt(temp1)
        A[1,2] = 0.0

        A[2,0] = R[2,0]/(A[0,0]+1e-20)
        A[2,1] = (R[2,1]-A[1,0]*A[2,0])/(A[1,1]+1e-20)
        temp1 = R[2,2]-A[2,0]*A[2,0]-A[2,1]*A[2,1]
        if temp1<0:
                temp1 = 0.0
        A[2,2] = np.sqrt(temp1)

        xu = yu[j,k]
        xv = yv[j,k]
        xw = yw[j,k]

        yu[j,k] =  A[0,0]*xu + A[0,1]*xv + A[0,2]*xw + uinr
        yv[j,k] =  A[1,0]*xu + A[1,1]*xv + A[1,2]*xw
        yw[j,k] =  A[2,0]*xu + A[2,1]*xv + A[2,2]*xw
	
def read_profile(profilefile,kma):

   profiledata = np.genfromtxt(profilefile, names=True,autostrip=True, comments='#')
   print profiledata.dtype.names
   npoints = profiledata.shape[0]

   for i in reversed(profiledata[0:npoints-2]): 
     profiledata = np.append(profiledata,i)

   profiledata['y'][npoints:] = (-(profiledata['y'][npoints:] - 1.0) + 1)
   profiledata['uv'][npoints:] = -profiledata['uv'][npoints:]
   
   zArray = profiledata['y']
   UArray = profiledata['U']
   uuArray = profiledata['uu']
   vvArray = profiledata['vv']
   wwArray = profiledata['ww']
   uwArray = profiledata['uv']
   
   zArray = (zArray-np.min(zArray ))/(np.max(zArray )-np.min(zArray ))
   zi = np.linspace(np.min(zArray ), np.max(zArray),kma)


   U = interpolate.splev(zi, interpolate.splrep(zArray, UArray, s=0), der=0)
   uu = interpolate.splev(zi, interpolate.splrep(zArray, uuArray, s=0), der=0)
   vv = interpolate.splev(zi, interpolate.splrep(zArray, vvArray, s=0), der=0)
   ww = interpolate.splev(zi, interpolate.splrep(zArray, wwArray, s=0), der=0)
   uw = interpolate.splev(zi, interpolate.splrep(zArray, uwArray, s=0), der=0)
 
   U[0]=U[-1]=0.   
   uu[0]=uu[-1]=0.   
   vv[0]=vv[-1]=0.   
   ww[0]=ww[-1]=0.   
   uw[0]=uw[-1]=0.   
   
   return U,uu,vv,ww,uw

def read_prf(profilefile,res,mdot,den,bulk_velocity,non_dim):

   f = open(profilefile,'r')
   not_data = True
   data_count = 0
   while not_data:
           line = f.readline()
           data_count += 1
           if line.startswith('data'):
                   line=line.strip()
                   not_data = False
                   data = line.split(',')

   #print data_count
   uudata=vvdata=wwdata=kdata=epsdata=sdrdata=-1 # use this as flag later        

   for i in range (1,len(data)):
           if data[i] == 'x':
                   xdata = i-1
           elif data[i] == 'y':
                   ydata = i-1
           elif data[i] == 'z':
                   zdata = i-1
           elif data[i] == 'u':
                   udata = i-1
           elif data[i] == 'v':
                   vdata = i-1
           elif data[i] == 'w':
                   wdata = i-1
           elif data[i] == 'k':
                   kdata = i-1
           elif data[i] == 'e':
                   epsdata = i-1
           elif data[i] == 'sdr':
                   sdrdata = i-1
           elif data[i] == 'uu':
                   uudata = i-1
           elif data[i] == 'vv':
                   vvdata = i-1
           elif data[i] == 'ww':
                   wwdata = i-1


   profiledata = np.loadtxt(profilefile,skiprows = data_count)
   npoints = profiledata.shape[0]

   xArray = profiledata[:,xdata]
   yArray = profiledata[:,ydata]
   zArray = profiledata[:,zdata]
   UArray = profiledata[:,udata]
   VArray = profiledata[:,vdata]
   WArray = profiledata[:,wdata]
   if uudata != -1:
           uuArray = profiledata[:,uudata]
   if vvdata != -1:  
           vvArray = profiledata[:,vvdata]
   if wwdata != -1:
           wwArray = profiledata[:,wwdata]
   if kdata != -1:
           kArray = profiledata[:,kdata]
   if epsdata != -1:
           epsArray = profiledata[:,epsdata]
   if sdrdata != -1:
           sdrArray = profiledata[:,sdrdata]

   # construct 2 vectors that run in plane
   # this may not be robust...
   # if points are not entered in a logical way then vector 1 & 2 COULD
   # be lying in the same direction

   x2 = xArray[1] - xArray[0] 
   y2 = yArray[1] - yArray[0] 
   z2 = zArray[1] - zArray[0] 

   x1 = xArray[-1] - xArray[0] 
   y1 = yArray[-1] - yArray[0] 
   z1 = zArray[-1] - zArray[0]

   # contruct normal vector
   xn = y1*z2-z1*y2
   yn = z1*x2-x1*z2
   zn = x1*y2-y1*x2

   # normalise!
   nnorm = np.sqrt(xn**2+yn**2+zn**2)
   xn = xn/nnorm
   yn = yn/nnorm
   zn = zn/nnorm

   # contruct centre point
   xc = (np.amax(xArray)+np.amin(xArray))/2
   yc = (np.amax(yArray)+np.amin(yArray))/2
   zc = (np.amax(zArray)+np.amin(zArray))/2

 #  print xn,yn,zn,xc,yc,zc
  # print np.amin(xArray)-np.amax(xArray), \
  #      np.amin(yArray)-np.amax(yArray), \
  #      np.amin(zArray)-np.amax(zArray)
   

   # the logic here is to reorient the plane such that it lies in the y,z axis
   # then it is easy to find the four corners and the dimension of the plane
   # but there may be a smarter way of doing this...

   # translate plane to 0,0
   xArray1 = xArray-xc
   yArray1 = yArray-yc
   zArray1 = zArray-zc

   # find (negative) rotation angle betwen 1,0,0 and normal
   theta = -np.arccos(xn)
   # find (negative) twist correction
   beta = -np.arctan2(zn,yn)

   # rotate about -nz,ny axis
   ## find transfomration matrix
   C = np.cos(theta)
   S = np.sin(theta)
   t = 1-C
   nx=0
   ny=-zn
   nz=yn

   T = np.matrix([[t*nx**2+C,t*nx*ny-S*nz,t*nx*nz+S*ny], \
                 [t*nx*ny+S*nz,t*ny**2+C,t*ny*nz-S*nx], \
                 [t*nx*nz-S*ny,t*ny*nz+S*nx,t*nz**2+C]])
   
   points = np.matrix([xArray1.T,yArray1.T,zArray1.T])

   
   # rotate!
   points = T*points

   
   # rotate about normal!
   ## find transfomration matrix
   C = np.cos(beta)
   S = np.sin(beta)
   t = 1-C
   nx=xn
   ny=yn
   nz=zn

   T = np.matrix([[t*nx**2+C,t*nx*ny-S*nz,t*nx*nz+S*ny], \
                 [t*nx*ny+S*nz,t*ny**2+C,t*ny*nz-S*nx], \
                 [t*nx*nz-S*ny,t*ny*nz+S*nx,t*nz**2+C]])
   
   # rotate!
   points = T*points

   # work out kma and jma
   yspan = np.amax(points[1,:]) - np.amin(points[1,:])
   zspan = np.amax(points[2,:]) - np.amin(points[2,:])
   kma = int(math.ceil(zspan/res))
   jma = int(math.ceil(yspan/res))

   print 'jma = ', jma ,' and kma = ', kma

   #tri = spatial.Delaunay(points.T)

   yarray = np.zeros(np.size(points[1,:]),dtype=np.float64)    
   zarray = np.zeros(np.size(points[1,:]),dtype=np.float64)    
   
   yArray = points[1,:]
   yarray[:] = points[1,:]
   zArray = points[2,:]
   zarray[:] = points[2,:]

   #yArray = (yArray-np.min(yArray ))/(np.max(yArray )-np.min(yArray ))
   yi = np.linspace(np.min(yArray ), np.min(yArray)+res*jma,jma)    
   
   #zArray = (zArray-np.min(zArray ))/(np.max(zArray )-np.min(zArray ))
   zi = np.linspace(np.min(zArray ), np.min(zArray)+res*kma,kma)    
   
   y,z = np.meshgrid(yi,zi)

   #matplotlib.pyplot.plot(yArray,zArray,'ro')
   #matplotlib.pyplot.plot(y,z,'bo')
   #matplotlib.pyplot.savefig('ha.png',bbox_inches='tight')
   
   #points = np.array([y,z])
   
   # interpolate to new grid!
   print np.min(points[1:,:].T), np.shape(UArray)

   #f=interpolate.CloughTocher2DInterpolator(points[1:,:].T,UArray,tol=1e-6)
   U=interpolate.griddata(points[1:,:].T,UArray,(y,z),fill_value=0.0,method='linear')
   V=interpolate.griddata(points[1:,:].T,VArray,(y,z),fill_value=0.0,method='linear')
   W=interpolate.griddata(points[1:,:].T,WArray,(y,z),fill_value=0.0,method='linear')
   #print U 
   if uudata != -1:
           uu=interpolate.griddata(points[1:,:].T,uuArray,(y,z),fill_value=0.0,method='linear')
           flag = uu<0
           uu[flag] = 0
   if vvdata != -1:  
           vv=interpolate.griddata(points[1:,:].T,vvArray,(y,z),fill_value=0.0,method='linear')
           flag = vv<0
           vv[flag] = 0
   if wwdata != -1:
           ww=interpolate.griddata(points[1:,:].T,wwArray,(y,z),fill_value=0.0,method='linear')
           flag = ww<0
           ww[flag] = 0
   if kdata != -1:
           k=interpolate.griddata(points[1:,:].T,kArray,(y,z),fill_value=0.0,method='linear')
           flag = k<0
           k[flag] = 0
   if epsdata != -1:
           eps=interpolate.griddata(points[1:,:].T,epsArray,(y,z),fill_value=0.0,method='linear')
           flag = eps<0
           eps[flag] = 0
   if sdrdata != -1:
           sdr=interpolate.griddata(points[1:,:].T,sdrArray,(y,z),fill_value=0.0,method='linear')
           flag = sdr<0
           sdr[flag] = 0
           # approximate epsilon
           eps = 0.09*k*sdr
	   eps[np.where(eps>100000000)]=0

   # scale variables if mass flow rate is different.
   if mdot !=0.0:
	# calc old mdot
	c_area = res**2
	A = c_area*(kma-1)*(jma-1)
	meanu = np.mean(U)
	meanv = np.mean(V)
	meanw = np.mean(W)
	udotn = meanu*xn+meanv*yn+meanw*zn
	mdot_old = udotn*A*den
	# find TI
        flag = (eps)>0 # flag to avoid 0 values from biasing average
	TI = np.sqrt(2./3.*k[flag])/np.sqrt(U[flag]**2+V[flag]**2+W[flag]**2)
	#find length scales
   	L = k[flag]**(3/2)/eps[flag]
	# find new velocities
	scale = mdot/mdot_old
	U=U*scale
	V=V*scale
	W=W*scale
	# find new k
	k[flag] = TI**2*(U[flag]**2+W[flag]**2+V[flag]**2)
	# find new epsilon
	eps[flag] = k[flag]**(3/2)/L
	print 'new maximum velocity is: ', np.amax(U),np.amax(V),np.amax(W)

   elif bulk_velocity != 1.0:
        meanu = np.mean(U)
        meanv = np.mean(V)
        meanw = np.mean(W)
        udotn = meanu*xn+meanv*yn+meanw*zn
        # find TI
        flag = (eps)>0 # flag to avoid 0 values from biasing average
        TI = np.sqrt(2./3.*k[flag])/np.sqrt(U[flag]**2+V[flag]**2+W[flag]**2)
        #find length scales
        L = k[flag]**(3/2)/eps[flag]
        # find new velocities
        scale = bulk_velocity/udotn
        U=U*scale
        V=V*scale
        W=W*scale
        # find new k
        k[flag] = TI**2*(U[flag]**2+W[flag]**2+V[flag]**2)
        # find new epsilon
        eps[flag] = k[flag]**(3/2)/L


   # calculate mean velocity gradients
   #dU = interpolate.CloughTocher2DInterpolator(points[1:,:].T,UArray,tol=1e-6)
   #dudy= dU(y,z)

   flag = np.where(eps==0.0)
   flag1 = np.where(U==0.0)
   U[flag] = 0
   V[flag] = 0
   W[flag] = 0
   k[flag] = 0
   eps[flag1]=0   

   dU = np.gradient(U,res)
   dV = np.gradient(V,res)
   dW = np.gradient(W,res)

   dUdy=dU[0]
   dUdz=dU[1]
   dVdy=dV[0]
   dVdz=dV[1]
   dWdy=dW[0]
   dWdz=dW[1]

   dUdy[flag]=0
   dUdz[flag]=0
   dVdy[flag]=0
   dVdz[flag]=0
   dWdy[flag]=0
   dWdz[flag]=0

   # smooth gradients (may not be necessary)
   dUdy1=dUdy.copy()
   dUdz1=dUdz.copy()
   dVdy1=dVdy.copy()
   dVdz1=dVdz.copy()
   dWdy1=dWdy.copy()
   dWdz1=dWdz.copy()
   for i in range(1,kma-1):
	for j in range (1,jma-1):
		dUdy[i,j]=np.mean(dUdy[i-1:i+1,j-1:j+1])
		dUdz[i,j]=np.mean(dUdz[i-1:i+1,j-1:j+1])
		dVdy[i,j]=np.mean(dVdy[i-1:i+1,j-1:j+1])
		dVdz[i,j]=np.mean(dVdz[i-1:i+1,j-1:j+1])
		dWdy[i,j]=np.mean(dWdy[i-1:i+1,j-1:j+1])
		dWdz[i,j]=np.mean(dWdz[i-1:i+1,j-1:j+1])

   if non_dim:
	y=y/np.amax(z)
	z=z/np.amax(z)

   plt.contourf(1,y,z,dUdy,100,'x','y','dudy','PODFS/dudy',figsize=(8,8*kma/jma))
   plt.close(1)
   plt.contourf(2,y,z,dUdz,100,'x','y','dudz','PODFS/dudz',figsize=(8,8*kma/jma))
   plt.close(2)
   plt.contourf(3,y,z,dVdy,100,'x','y','dvdy','PODFS/dvdy',figsize=(8,8*kma/jma))
   plt.close(3)
   plt.contourf(4,y,z,dVdz,100,'x','y','dvdz','PODFS/dvdz',figsize=(8,8*kma/jma))
   plt.close(4)
   plt.contourf(5,y,z,dWdy,100,'x','y','dwdy','PODFS/dwdy',figsize=(8,8*kma/jma))
   plt.close(5)
   plt.contourf(6,y,z,dWdz,100,'x','y','dwdz','PODFS/dwdz',figsize=(8,8*kma/jma))
   plt.close(6)
   plt.contourf(7,y,z,U,100,'x','y','U','PODFS/U',figsize=(8,8*kma/jma))
   plt.close(7)
   plt.contourf(8,y,z,V,100,'x','y','V','PODFS/V',figsize=(8,8*kma/jma))
   plt.close(8)
   plt.contourf(9,y,z,W,100,'x','y','W','PODFS/W',figsize=(8,8*kma/jma))
   plt.close(9)
   plt.contourf(10,y,z,eps,100,'x','y','eps','PODFS/eps',figsize=(8,8*kma/jma))
   plt.close(10)
   plt.contourf(10,y,z,k,100,'x','y','k','PODFS/k',figsize=(8,8*kma/jma))
   plt.close(10)

   # calculate dudx from incompressibility (only an approximation!)
   dUdx = -dVdy -dWdz

   # assume other gradients in flow direction are zero! 
   # (might be a better option..)

   dVdx = np.zeros((kma,jma),dtype=np.float64)
   dWdx = np.zeros((kma,jma),dtype=np.float64)

   # calculate lengthscale and for now average in space
   # lengthscale from eps/sdr and k is calculated in these 3 lines:
   #flag = (eps)>0 # flag to avoid 0 values from biasing average
   #L = k[flag]**(3/2)/eps[flag]
   #L1 = np.mean(L) # this value lools small, revise??
   
   # lengthscale is approximated as 0.07*hydraulic dimeter
   B=2*np.amax(points[1,:])
   C=2*np.amax(points[2,:])
   L = 0.07*2*B*C/(B+C) 
   # much larger lengthscale! makes things more expensive!
        
   #find lengthscale in terms of res
   lnx = math.ceil(L/res)

   # commented out section is an attempt to use an algebraic
   # stress model to compute the reynolds stresses
   # The method appears to be very sensitive to the calculation
   # of the velocity gradients.

   '''

   # apply algebriac stress model
   uu = np.zeros((kma,jma),dtype=np.float64)
   vv = np.zeros((kma,jma),dtype=np.float64)
   ww = np.zeros((kma,jma),dtype=np.float64)
   uv = np.zeros((kma,jma),dtype=np.float64)
   uw = np.zeros((kma,jma),dtype=np.float64)
   vw = np.zeros((kma,jma),dtype=np.float64)
   for i in range (0,kma):
    for j in range (0,jma):
     if eps[i,j]>0.0:
	dudx=dUdx[i,j]
	dvdx=dVdx[i,j]
	dwdx=dWdx[i,j]
	dudy=dUdy[i,j]
	dvdy=dVdy[i,j]
	dwdy=dWdy[i,j]
	dudz=dUdz[i,j]
	dvdz=dVdz[i,j]
	dwdz=dWdz[i,j]
   	# calculate production tensor matrix
  	Pij = np.matrix([[-2*dudx,-2*dudy,-2*dudz,0,0,0],\
                    [-dvdx,-(dvdy+dudx),-dvdz,-dudy,-dudz,0],\
                    [-dwdx,-dwdy,-(dwdz+dudx),-1,-dudy,-dudz],\
                    [0,-2*dvdx,0,-2*dvdy,-2*dvdz,0],\
                    [0,-dwdx,-dvdx,-dwdy,-(dwdz+dvdy),-dvdz],\
                    [0,0,-2*dwdx,0,-2*dwdy,-2*dwdz]])

   	nu_t = 0.09*k[i,j]**2/eps[i,j]

   	Pk = nu_t*(2*dudx**2+dudy**2+dudz**2 \
	      +dvdx**2+2*dvdy**2+dvdz**2 \
	      +dwdx**2+dwdy**2+2*dwdz**2 \
	      +dudy*dvdx + dudz*dwdx \
	      +dvdx*dudy + dvdz*dwdy \
	      +dwdx*dudz + dwdy*dvdz)

   
   	LHS = Pk-eps[i,j]+1.8*eps[i,j]
   	#print nu_t,Pk
	LHS = LHS - (1.-0.6)*Pij
	#print LHS

   	RHS = np.zeros((kma,jma),dtype=np.float64)
   	RHS = 2./3.*(0.6*Pk+(1.8-1.)*eps[i,j])
	#print RHS
   	RHS = np.matrix([[RHS],[0],[0],[RHS],[0],[RHS]])

	try:
   		Rij = np.linalg.solve(LHS,RHS)
	except: 
		Rij = [0.,0.,0.,0.,0.,0.]

	#print Rij
	#print np.shape(RHS),np.shape(LHS)

	uu[i,j] = Rij[0]
	uv[i,j] = Rij[1]
	uw[i,j] = Rij[2]
	vv[i,j] = Rij[3]
	vw[i,j] = Rij[4]
	ww[i,j] = Rij[5]

   # filter
   val= 2.0
   meanuu = np.abs(np.mean(uu))
   print meanuu
   flag = np.where(np.abs(uu)>100000)#val*meanuu)
   uu[flag] = 0.0
   flag = np.where((uu)<0)
   ww[np.where(ww<0.0)] = 0.0
   vv[np.where(vv<0.0)] = 0.0
   meanuu = np.mean(uu)

   uu1=uu.copy()
   vv1=vv.copy()
   ww1=ww.copy()
   uv1=uv.copy()
   uw1=uw.copy()
   vw1=vw.copy()
   for i in range(1,kma-1):
        for j in range (1,jma-1):
                uu[i,j]=np.mean(uu1[i-1:i+1,j-1:j+1])
                vv[i,j]=np.mean(vv1[i-1:i+1,j-1:j+1])
                ww[i,j]=np.mean(ww1[i-1:i+1,j-1:j+1])
                uv[i,j]=np.mean(uv1[i-1:i+1,j-1:j+1])
                uw[i,j]=np.mean(uw1[i-1:i+1,j-1:j+1])
                vw[i,j]=np.mean(vw1[i-1:i+1,j-1:j+1])

   '''

   # try eddy viscosity model instead:
   nu_t = np.zeros((kma,jma),dtype=np.float64)
   flag = np.where(eps>0)
   nu_t[flag] = 0.09*k[flag]**2/eps[flag]
   uu = -2.*nu_t*dUdx+2./3.*k
   vv = -2.*nu_t*dVdy+2./3.*k
   ww = -2.*nu_t*dWdz+2./3.*k
   uv = -nu_t*(dUdy+dVdx)
   uw = -nu_t*(dUdz+dWdx)
   vw = -nu_t*(dVdz+dWdy)

   # filter for negative normal stress components
   uu[np.where(uu<0.0)]=0.0
   vv[np.where(vv<0.0)]=0.0
   ww[np.where(ww<0.0)]=0.0

   plt.contourf(11,y,z,uu,100,'x','y','uu','PODFS/uu',figsize=(8,8*kma/jma))
   plt.close(11)
   plt.contourf(12,y,z,vv,100,'x','y','vv','PODFS/vv',figsize=(8,8*kma/jma))
   plt.close(12)
   plt.contourf(13,y,z,ww,100,'x','y','ww','PODFS/ww',figsize=(8,8*kma/jma))
   plt.close(13)
   plt.contourf(14,y,z,uv,100,'x','y','uv','PODFS/uv',figsize=(8,8*kma/jma))
   plt.close(14)
   plt.contourf(15,y,z,uw,100,'x','y','uw','PODFS/uw',figsize=(8,8*kma/jma))
   plt.close(15)
   plt.contourf(16,y,z,vw,100,'x','y','vw','PODFS/vw',figsize=(8,8*kma/jma))
   plt.close(16)

   U=np.flip(U,0)
   V=np.flip(V,0)
   W=np.flip(W,0)
   uu=np.flip(uu,0)
   vv=np.flip(vv,0)
   ww=np.flip(ww,0)
   uv=np.flip(uv,0)
   uw=np.flip(uw,0)
   vw=np.flip(vw,0)

   return U.T,V.T,W.T,uu.T,vv.T,ww.T,uv.T,uw.T,vw.T,lnx \
		,kma,jma,xn,yn,zn,xc,yc,zc


def build_profile(mean_profile,turb_profile,bulk_velocity,turbulence_intensity,kma):

   if mean_profile in ['hyperbolic-tangent','double-hyperbolic-tangent', \
			'circular-hyperbolic-tangent','ring-hyperbolic-tangent']:
           
           y = np.linspace(-0.5,0.5,kma)
           U = bulk_velocity/2*(1.+np.tanh(10.*(-np.abs(y)+0.5)))
           #print np.tanh(np.abs(y))
   else:
	raise Exception('Invalid mean_profile chosen, type \'python digitalfilters.py -h\' for help.')

   if turb_profile == 'top-hat':
           uu = (turbulence_intensity*U)**2
           vv = (turbulence_intensity*U)**2
           ww = (turbulence_intensity*U)**2
           uw = 0.0*U
   elif turb_profile == 'none':
	   uu = 0.0
           vv = 0.0
           ww = 0.0
           uw = 0.0
   else:
	raise Exception('Invalid turb_profile chosen, type \'python digitalfilters.py -h\' for help.')

   return U,uu,vv,ww,uw


def main():
   description = """ LES Inflow Generator after Klein et.al. """

   parser = OptionParser(option_class=MultipleOption,
                          usage='usage: %prog [options]',
                          version='%s %s' % (PROG, VERSION),
                          description=description)


   parser.add_option("-i", "--inputfile", dest="profilefile", default='none',
			 help="1d turbulent profile file", metavar="FILE")

   parser.add_option("-p", "--mean_profile", dest="mean_profile", default='hyperbolic-tangent',
			 help="What kind of mean flow profile would you like? Options: hyperbolic-tangent, double-hyperbolic-tangent, ring-hyperbolic-tangent, circular-hyperbolic-tangent. Using the -i option along with this will adapt the -i profile to the shape of the -p profile ie. plane jet, square jet, round jet, annulus jet.", metavar="STRING")
   
   parser.add_option("--turb_profile", dest="turb_profile", default='top-hat',
			 help="What kind of turbulence flow profile would you like? Options: top-hat, none", metavar="STRING")

   parser.add_option("--U0","--bulk_velocity", type="float",dest="bulk_velocity", default=1.0,
			 help="What is the bulk velocity magnitude?, this \
			option can also be used to scale the velocities \
			and turbulent quantities similar to the massflow \
			option is using a .prf profile.", metavar="NUM")

   parser.add_option("--u_dash", type="float",dest="turbulence_intensity", default='0.02',
			 help="What is the desired u'/U0 value with u'=v'=w'?", metavar="NUM")

   parser.add_option("-n", "--nsteps", type="int", dest="nsteps",default=20,
			 help="number of steps", metavar="INT")
			 
   parser.add_option("-l","--lengthscale", type="float", dest="lengthscale", default=3.0,
		       help="turbulent lengthscale in terms of grid spacing ", metavar="NUM")
			 
   parser.add_option("-f", "--fwidth", type="float", dest="fwidth",default=2.0,
			 help="half filter width in lengthscales, should be greater than 2", metavar="NUM")
   			 
   parser.add_option("-k", "--nk", type="int", dest="kma",default=11,
			 help="number of points in k (wall-normal) direction", metavar="INT")

   parser.add_option("-j", "--nj", type="int", dest="jma",default=10,
			 help="number of points in j (spanwise) direction", metavar="INT")
   
   parser.add_option("-t", "--dt", type="float", dest="dt",default=0.0,
			 help="time step (s)", metavar="NUM")

   parser.add_option("-m", "--nm", type="int", dest="nm",default=20,
			 help="number of POD modes", metavar="INT")

   parser.add_option("-e", "--et", type="float", dest="et",default=0.9,
			 help="target energy for Fourier reconstruction", metavar="NUM")

   parser.add_option("-v", "--verbose", dest="verbose",default=False,
			 help="Save the mean flow, POD spatial and temporal modes?",action='store_true')
   
   parser.add_option("--non_dim", dest="non_dim",default=False,
			 help="Non-dimensionalise lengths if using .prf",action='store_true')

   parser.add_option("-r", "--resolution", type="float", dest="res",default=0.1,
			 help="plane resolution in meters per grid point", metavar="NUM")

   parser.add_option("--nx", type="float", dest="nx",default=1.0,
			 help="plane normal direction, x-component", metavar="NUM")

   parser.add_option("--ny", type="float", dest="ny",default=0.0,
			 help="plane normal direction, y-component", metavar="NUM")

   parser.add_option("--nz", type="float", dest="nz",default=0.0,
			 help="plane normal direction, z-component", metavar="NUM")


   parser.add_option("--ox", type="float", dest="ox",default=0.0,
			 help="plane origin, x-component", metavar="NUM")

   parser.add_option("--oy", type="float", dest="oy",default=0.0,
			 help="plane origin, y-component", metavar="NUM")

   parser.add_option("--oz", type="float", dest="oz",default=0.0,
			 help="plane origin, z-component", metavar="NUM")

   parser.add_option("--rotate", type="float", dest="rot",default=0.0,
			 help="rotate plane about its normal (degrees)", metavar="NUM")

   parser.add_option("--ring", type="float", dest="ring",default=0.5,
                         help="inner diameter of ring as a proportion of the outer \
			diameter if using the ring-hyperbolic-tangent option", metavar="NUM")

   parser.add_option("--massflow", type="float", dest="mdot",default=0.0,
                         help="If using a .prf file, the velocities can be  \
			scaled to achieve the desired mass flow rate. Mean \
			velocities are scaled equally. k and epsilon/sdr \
			are scaled to maintain the same turbulence intensity \
			and length scale. Requires density!", metavar="NUM")

   parser.add_option("--density", type="float", dest="den",default=0.0,
                         help="If massflow is specified, a density must  \
                        also be specified.", metavar="NUM")



   if len(sys.argv) == 1:
	       parser.parse_args(['--help'])
    
   (options, args) = parser.parse_args()

   profilefile = options.profilefile   

   if not os.path.exists('PODFS'):
       os.mkdir('PODFS')
   V=W=0
 
   # if profilefile == none then build profile from options...
   if profilefile == 'none':
	mean_profile = options.mean_profile
	turb_profile = options.turb_profile
	bulk_velocity = options.bulk_velocity
	turbulence_intensity = options.turbulence_intensity
   else: # allow adoption of profile to different shapped inlets
	mean_profile = options.mean_profile 	

   nsteps = options.nsteps  

   #kma = 11
   #jma = 10 #For generating 2D inflow data set nfy=0 and jma=1
   
   kma = options.kma
   jma = options.jma

   mdot = options.mdot
   den = options.den

   res = options.res   
   dt = options.dt

   inner_d = options.ring

   lnx = options.lengthscale 
   lny = lnx
   lnz = lnx
   
   # calculate required filter size from lengthscale and fwidth
   nf = int(math.ceil(options.fwidth*options.lengthscale))
   nfx = nf
   nfy = nf
   nfz = nf
  
   nx1 = options.nx
   ny1 = options.ny
   nz1 = options.nz
   #normalise
   nx = nx1/np.sqrt(nx1**2+ny1**2+nz1**2)
   ny = ny1/np.sqrt(nx1**2+ny1**2+nz1**2)
   nz = nz1/np.sqrt(nx1**2+ny1**2+nz1**2)

   ox = options.ox
   oy = options.oy
   oz = options.oz
 
   if (options.profilefile and os.path.isfile(profilefile)):
        if profilefile.endswith('.prf'): # precise 2D profile file
                U,V,W,uu,vv,ww,uv,uw,vw,lnx,kma,jma \
		,nx,ny,nz,ox,oy,oz = \
		read_prf(profilefile,res,mdot,den,options.bulk_velocity \
			,options.non_dim)
		lny=lnz=lnx
        else:
                U,uu,vv,ww,uw = read_profile(profilefile,kma)
   else: # contruct profile
       #exit()
       U,uu,vv,ww,uw = build_profile(mean_profile,turb_profile,bulk_velocity,turbulence_intensity,kma)
   if dt == 0.: #calculate dt from res and U
	   flag = np.where(U**2+V**2+W**2!=0)
           dt = res/np.mean(U[flag])
           print 'timestep set to: ', dt , ' seconds'
   else: #rescale lnx and nfx
	   flag = np.where(U**2+V**2+W**2!=0)
           dt1 = res/np.mean(U[flag])
           factor = dt1/dt
           lnx = factor*lnx
           nfx = int(math.ceil(float(options.fwidth)*lnx))
           print 'Lengthscale in x-direction set to: ' , lnx , 'grid points'
           print 'Filter width in x-direction set to: ' , nfx , 'grid points'
     


   #half filter length nf? should be greater than 2 ln?
   #For generating 2D inflow data set nfy=0 and jma=1
   
   #nfx = 3
   #nfy = 3 
   #nfz = 3

   nm = options.nm
   et = options.et
   verbose = options.verbose
   
   #number of different filters (normally 1)
   na = 1 

   # Variance of random numbers should be equal to 1
   # If using uniform pdf then range should be -sqrt(3) < x < sqrt(3)
   # I think...

   pdfr = np.sqrt(3.0)
 
   # check values of reynolds stress are positive
   if not profilefile.endswith('.prf'):
     for k in range (0,kma):
       if uu[k] < 0.0:
           uu[k] = 0.0
       if vv[k] < 0.0:
           vv[k] = 0.0
       if ww[k] < 0.0:
           ww[k] = 0.0

   a = np.zeros((na,nfx*2+1,nfy*2+1,nfz*2+1))
  

   coeff3D(a,nfx,nfy,nfz,lnx,lny,lnz)

   xu = np.random.uniform(low=-pdfr, high=pdfr, 
                          size=(nfx*2+1,nfy*2+1+jma,nfz*2+1+kma)) #random fields
   xv = np.random.uniform(low=-pdfr, high=pdfr, 
                          size=(nfx*2+1,nfy*2+1+jma,nfz*2+1+kma)) #random fields
   xw = np.random.uniform(low=-pdfr, high=pdfr, 
                          size=(nfx*2+1,nfy*2+1+jma,nfz*2+1+kma)) #random fields

   yu = np.zeros((jma,kma))
   yv = np.zeros((jma,kma))
   yw = np.zeros((jma,kma))
   
   
   #pvar = np.zeros(nsteps)
   
   # define input object
   i_d = obj()
   i_d.kma = kma
   i_d.jma = jma
   i_d.ns = nsteps
   i_d.dt = dt
   i_d.nm = nm
   i_d.et = et
   i_d.rot = options.rot
   i_d.t_o = [ox,oy,oz]
   i_d.res = res
   i_d.n = [nx,ny,nz]
   num_points = jma*kma
   i_d.num_points = num_points
   i_d.is_POD_var_vec = False
   grid = pod.make_inflow_plane(i_d)
   i_d.grid=grid

   A = np.zeros((jma*kma*3,nsteps),dtype=np.float64)
 
   if not os.path.exists('PODFS'):
       os.mkdir('PODFS')


   for i in range(nsteps):
   
     print 'Generating inlet ' ,i+1, ' of ', nsteps
 
    
     filter3D(xu,yu,a,jma,kma,nfx,nfy,nfz) 
     filter3D(xv,yv,a,jma,kma,nfx,nfy,nfz)
     filter3D(xw,yw,a,jma,kma,nfx,nfy,nfz)        
  
     if profilefile.endswith('.prf'):
 	adapt2prf(yu,yv,yw,U,V,W,uu,vv,ww,uv,uw,vw,jma,kma)
     elif mean_profile in ['double-hyperbolic-tangent', \
	'ring-hyperbolic-tangent','circular-hyperbolic-tangent']:
     	adapt2d(yu,yv,yw,U,uu,vv,ww,uw,jma,kma,mean_profile,inner_d)
     else:	
	adapt1d(yu,yv,yw,U,uu,vv,ww,uw,jma,kma)
     

     xu = np.roll(xu, -1, axis=0)
     xv = np.roll(xv, -1, axis=0)
     xw = np.roll(xw, -1, axis=0)

     
   
     xu[nfx*2,:,:] = np.random.uniform(low=-pdfr, high=pdfr, 
                          size=(nfy*2+1+jma,nfz*2+1+kma)) #random fields

     xv[nfx*2,:,:] = np.random.uniform(low=-pdfr, high=pdfr, 
                          size=(nfy*2+1+jma,nfz*2+1+kma)) #random fields
			  
     xw[nfx*2,:,:] = np.random.uniform(low=-pdfr, high=pdfr, 
                          size=(nfy*2+1+jma,nfz*2+1+kma)) #random fields

     #pvar[i] = yu[int(jma/2),int(kma/2)]
     
     A[0:jma*kma,i] = yu.reshape(jma*kma)
     A[jma*kma:2*jma*kma,i] = yv.reshape(jma*kma)
     A[2*jma*kma:3*jma*kma,i] = yw.reshape(jma*kma)
   
     if verbose: # write fields!
         i_d.time = i*dt
         pod.save_plane(A[:,i],i_d)

   # prepare for POD
   if verbose:
       nmw = nm
   else:
       nmw = 0
   i_d.verbose = verbose

   # subtract mean
   mean_field = np.mean(A,1)
   i_d.mean_field = mean_field
   for j in range(0,nsteps):
       A[:,j] = A[:,j] - mean_field[:]



   # Do POD
   pod.POD(A,nsteps,num_points,3,'false',[],'PODFS/','false',1.0e-15,nm,nmw,'false','false',grid,mean_field,dt,'velocity',1,nsteps,1,1,i_d)


   # Compute fourier coefficients
   pod.fourier_coefficients(i_d)

   # Save spatial modes as prfs
   pod.pod2prf(i_d)


if __name__ == '__main__':
  main()
