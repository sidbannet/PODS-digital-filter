# nicks plot library

import matplotlib as mpl
from matplotlib.patches import Arc
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import vtk
import vtk.util.numpy_support as VN
import numpy as np
import os
import socket

def hold(on):
    plt.hold(on)

def show():
    plt.show()

def axis(var):
    plt.axis(var)

def X_is_running():
    # get the display from the environment
    display_env = os.environ['DISPLAY']
    
    # parse the display string
    display_host, display_num = display_env.split(':')
    display_num_major, display_num_minor = display_num.split('.')

    # calculate the port number
    display_port = 6000 + int(display_num_major)

    # attempt a TCP connection
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        sock.connect((display_host, display_port))
    except socket.error:
        return False
    finally:
        sock.close()
    return True

#=======================================================================
def bar(fig,x,y,x_axis,y_axis,filename):

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.figure(fig)
    f=plt.bar(x, y)
    plt.xlabel('$'+x_axis+'$')
    plt.ylabel('$'+y_axis+'$')

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def sbar(fig,x,y,x_axis,y_axis,filename): #stacked bar

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    blah = y.shape
    colours = ('b','g','r','c','m', 'y', 'k', 'w','b','g','r','c','m', 'y', 'k', 'w')
    colours = colours[0:blah[0]]
    sum1 = np.zeros(blah[1])
    plt.figure(fig)
    
    for i in range (0,blah[0]):
        if i ==0:
            f=plt.bar(x, y[i,:],color=colours[i])
	if i<8:
            f=plt.bar(x,y[i,:],color=colours[i],bottom=sum1)
        else:
            f=plt.bar(x,y[i,:],color=colours[i],bottom=sum1,hatch='//')
        sum1 += y[i,:]
    plt.xlabel('$'+x_axis+'$')
    plt.ylabel('$'+y_axis+'$')

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def contourf(fig,x,y,z,levels,x_axis,y_axis,c_bar,filename):

    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')

    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

    plt.figure(fig)
    f=plt.contourf(x, y, z,levels)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    cbar=plt.colorbar(f)
    cbar.ax.set_ylabel('$'+c_bar+'$')
    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def contourfquiver(fig,x,y,z,u,v,levels,x_axis,y_axis,c_bar,filename):

    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')

    plt.figure(fig)
    f=plt.contourf(x, y, z,levels,cmap=plt.cm.RdBu)
    plt.quiver(x,y,u,v)
    plt.xlabel('$'+x_axis+'$')
    plt.ylabel('$'+y_axis+'$')
    cbar=plt.colorbar(f)
    cbar.ax.set_ylabel(c_bar)
    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def contourfcontour(fig,x,y,z,z2,levels,x_axis,y_axis,c_bar,filename):

    plt.rc('text', usetex=True)
    plt.rc('font', family='normal')
    plt.rc('font', size=22)
    origin = 'lower'
    plt.figure(fig)
    f=plt.contourf(x, y, z,levels,cmap=plt.cm.RdBu)
    CS2 = plt.contour(x,y,z2, 7,
                        colors = 'k',
                        hold='on')

    plt.xlabel('$'+x_axis+'$')
    plt.ylabel('$'+y_axis+'$')
    cbar=plt.colorbar(f)
    cbar.ax.set_ylabel('$'+c_bar+'$')
    #plt.contour(x,y,z2,10,colours='green')
    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def timeseries(fig,y,t,y_axis,filename):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    font = {'family' : 'normal',
        'weight' : 'bold',
            'size'   : 12}

    plt.rc('font', **font)

    major_ticks = np.linspace(np.amin(t),np.amax(t),5)

    plt.figure(fig)
    f=plt.plot(t,y)
    plt.xticks(major_ticks)
    plt.xlabel('$Time\, (s)$')
    plt.ylabel('$'+y_axis+'$')
    plt.grid(True)

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

    #plt.close(fig)

def plot(fig,y,t,y_axis,x_axis,filename):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

    plt.rc('font', **font)
    major_ticks = np.linspace(np.amin(t),np.amax(t),3)
    y_ticks = np.linspace(np.amin(y),np.amax(y),5)

    plt.figure(fig)
    f=plt.plot(t,y,linewidth=4)
    plt.xticks(major_ticks)
    plt.yticks(y_ticks)
    plt.xlabel('$'+x_axis+'$')
    plt.ylabel('$'+y_axis+'$')

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

 
def scatter(fig,y,t,y_axis,x_axis,filename):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')

    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}


    plt.figure(fig)
    f=plt.scatter(t,y,s=100)
    plt.xlabel('$'+x_axis+'$',fontsize=22)
    plt.ylabel('$'+y_axis+'$',fontsize=22)

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def scatter2(fig,y,t,y_axis,x_axis,filename):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')

    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}


    plt.figure(fig)
    f=plt.scatter(t,y,s=100,c='g')
    plt.xlabel('$'+x_axis+'$',fontsize=22)
    plt.ylabel('$'+y_axis+'$',fontsize=22)

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')


def scattereq(fig,y,t,y_axis,x_axis,filename,c1):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    plt.figure(fig,figsize=(5,5))
    f=plt.scatter(t,y,c=c1)
    plt.xlabel('$'+x_axis+'$')
    plt.ylabel('$'+y_axis+'$')
    
    xscale = np.amax(t)-np.amin(t)
    yscale = np.amax(y)-np.amin(y)
    if yscale > xscale:
        xmin = np.amin(t)-.1*yscale
        ymin = np.amin(y)-.1*yscale
        xmax = np.amin(t)+1.1*yscale
        ymax = np.amax(y)+.1*yscale
    else:
        xmin = np.amin(t)-.1*xscale
        ymin = np.amin(y)-.1*xscale
        xmax = np.amax(t)+.1*xscale
        ymax = np.amin(y)+1.1*xscale

    axes = plt.gca()
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    major_ticks = np.linspace(xmin, xmax,5)
    plt.xticks(major_ticks)
   # plt.axis('equal')
    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def arrow(fig,y,t,dx,dy,y_axis,x_axis,filename):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    plt.figure(fig,figsize=(5,5))
    for i in range (0,len(dy)):
        f=plt.arrow(t[i],y[i],dy[i],dx[i],length_includes_head=True,head_width=0.001,head_length=0.001,width=0.0001,fc='k',ec='k')
    plt.xlabel('$'+x_axis+'$')
    plt.ylabel('$'+y_axis+'$')
    plt.axis('equal')
    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def arcarrow(fig,y,t,dx,dy,y_axis,x_axis,filename):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')

    ax=plt.axes()
    fig=plt.figure(fig,figsize=(5,5))
    ax=plt.gca() #fig.add_subplot(111)
    #ax = fig.add_axes([0,0,1,1])
    
    l = []
    for i in range (0,len(dy)):
        t2 = t[i]+dy[i]
        x2  = y[i]+dx[i]
        #print t2,x2
        #print np.arctan2(y[i],dx[i])*180/np.pi,np.arctan2(x2,t2)*180/np.pi
        l.append(Arc((0,0),2*t[i],2*x2,theta1=np.arctan2(y[i],dx[i])*180/np.pi,theta2=np.arctan2(x2,t2)*180/np.pi,edgecolor='k'))
        #print arc
        ax.add_patch(l[i])
        f=plt.arrow(t2+0.001,x2,-0.0001*x2,0.0001*t2,length_includes_head=False,head_width=0.001,head_length=0.001,width=0.0001,fc='k',ec='k')
    #ax.plot([0,1],[0,1],'k')
    plt.xlabel('$'+x_axis+'$')
    plt.ylabel('$'+y_axis+'$')

    #plt.axis('equal')
    #plt.show()
    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')


def hist(fig,y,t,x_axis,filename):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 55}

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)



    plt.figure(fig)
    ax = plt.gca()
    plt.hist(y,t)
    plt.xlabel('$'+x_axis+'$',fontsize=22)
    plt.ylabel('$Counts$',fontsize=22)

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')



def logscatter(fig,y,t,y_axis,x_axis,filename):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 55}

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)



    plt.figure(fig)
    ax = plt.gca()
    f=plt.scatter(t,y)
    plt.xlabel('$'+x_axis+'$',fontsize=22)
    plt.ylabel('$'+y_axis+'$',fontsize=22)

    plt.yscale('log')
    plt.xscale('log')

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')




def logscatter2(fig,y,t,y_axis,x_axis,x_range,y_range,filename):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 55}

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.figure(fig)
    ax = plt.gca()
    f=plt.scatter(t,y)
    plt.xlabel('$'+x_axis+'$',fontsize=22)
    plt.ylabel('$'+y_axis+'$',fontsize=22)

    plt.yscale('log')
    plt.xscale('log')

    plt.xlim(x_range[0],x_range[1])
    plt.ylim(y_range[0],y_range[1])

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')



 
def plotwithpoint(fig,y,t,px,py,y_axis,x_axis,filename):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    plt.figure(fig)
    f=plt.scatter(px,py)
    plt.plot(t,y)
    plt.xlabel('$'+x_axis+'$')
    plt.ylabel('$'+y_axis+'$')

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')


def semilogx(fig,y,t,y_axis,x_axis,filename):
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    plt.figure(fig)
    f=plt.semilogx(t,y)
    plt.xlabel('$'+x_axis+'$')
    plt.ylabel('$'+y_axis+'$')

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def PSD(fig,y,t,y_axis,filename):
   # plt.rc('text', usetex=True)
   # plt.rc('font', family='serif')
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

    plt.rc('font', **font)
    plt.figure(fig)
    try:
    	f=plt.loglog(t,y,linewidth=4)
    	plt.xlabel('$Frequency \, (Hz)$')
    	plt.ylabel('$'+y_axis+'$')
    	plt.grid(True)
    except:
	print 'PSD has no positive values'

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

    plt.close(fig)

def WNS(fig,y,t,y_axis,filename):
   # plt.rc('text', usetex=True)
   # plt.rc('font', family='serif')
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

    end1 = np.size(y)
    p1 = int(1./100.*end1)
    p2 = int(1./3.*end1)
    #print end1,p1,p2
    x = np.linspace(t[p1],t[p2],p2-p1+1)
    c = (y[p1]+0.5)*x[p1]**(5./3.) + 0.1
    y2 = c*x**(-5./3.)

    plt.rc('font', **font)
    plt.figure(fig)
    f=plt.loglog(t,y,linewidth=4)
    f=plt.loglog(x,y2,linewidth=4)
    plt.xlabel('$Wavenumber \, (1/m)$')
    plt.ylabel('$'+y_axis+'$')
    plt.grid(True)

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

    plt.close(fig)

def semilogxold(fig,y,t,y_axis,filename):
   # plt.rc('text', usetex=True)
   # plt.rc('font', family='serif')
    plt.figure(fig)
    f=plt.semilogx(t,y)
    plt.xlabel('$Hz$')
    plt.ylabel('$'+y_axis+'$')

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def PSD3(fig,y,y1,y2,t,y_axis,filename):
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.figure(fig)
    f=plt.loglog(t,y,t,y1,t,y2)
    plt.xlabel('$St$')
    plt.ylabel('$'+y_axis+'$')

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def eigs(fig,y,y_axis,filename):
    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
    plt.figure(fig)
    f=plt.loglog(y,'.-')
    plt.xlabel(r"Mode Number")
    plt.ylabel('$'+y_axis+'$')

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def convergence(fig,y,y_axis,filename):
   # plt.rc('text', usetex=True)
   # plt.rc('font', family='serif')
    plt.figure(fig)
    f=plt.plot(y,'.-')
 #   plt.xlabel(r"Timestep")
  #  plt.ylabel('$'+y_axis+'$')

    if (filename != ' '):
        plt.savefig(filename+'.png', bbox_inches='tight')

def close(ifig):
    plt.close(ifig)
        

def scatter_3d(fig,xs,ys,zs,labels):
    fig = plt.figure(fig)
    ax = fig.add_subplot(111, projection = '3d')
    #xs, ys, zs = np.split(points, 3, axis=1)
    sc = ax.scatter(xs,ys,zs)

    # if this code is placed inside a function, then
    # we must use a predefined global variable so that
    # the update function has access to it. I'm not
    # sure why update_positions() doesn't get access
    # to its enclosing scope in this case.
    global labels_and_points
    labels_and_points = []

    for txt, x, y, z in zip(labels, xs, ys, zs):
        x2, y2, _ = proj3d.proj_transform(x,y,z, ax.get_proj())
        label = plt.annotate(
            txt, xy = (x2, y2), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        labels_and_points.append((label, x, y, z))

    def update_position(e):
        for label, x, y, z in labels_and_points:
            x2, y2, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
            label.xy = x2,y2
            label.update_positions(fig.canvas.renderer)
        fig.canvas.draw()

    fig.canvas.mpl_connect('motion_notify_event', update_position)

    plt.show()



#=========================================================================
# VTK PLot functions
#=========================================================================

def vtkcontourf(cutMapper,pl3d_output,src_output,var_name,levels,snapshot_dir,plane_normal,scalarRange):


#    if not X_is_running():
 #       return

    #off-screen rendering
    #gfac = vtk.vtkGraphicsFactory
    #gfac.SetOffScreenOnlyMode(1)
    #gfac.SetUseMesaClasses(1)
    #im_fac = vtk.vtkImagingFactory
    #im_fac.SetUseMesaClasses(1)   

 
    pl3d_output.GetPointData().SetActiveScalars(var_name)
    if (scalarRange == (0,0)):

        #print scalarRange

        # this prevents zero values from setting colourbar
        scalarnump =  VN.vtk_to_numpy(src_output.GetOutput().GetPointData().GetScalars(var_name))
        minval = np.amin(scalarnump)
        try:
            if (minval >= 0.0):
                minval = np.amin(scalarnump[scalarnump > 0.0])
        except Exception:
            pass
        maxval = np.amax(scalarnump)
        scalarRange = (minval,maxval)
        #print scalarRange


        #scalarRange =pl3d_output.GetPointData().GetScalars(var_name).GetRange()

    cutMapper.SetScalarRange(scalarRange)
    #print cutMapper
    print scalarRange
    cutActor = vtk.vtkActor()
    cutActor.SetMapper(cutMapper)
    
    lut = vtk.vtkLookupTable() #MakeLUT()
    lut.SetNumberOfTableValues(levels)
    lut.SetHueRange(0.667, 0.0)
    #lut.SetHueRange(0.801282320, 0.149576098)
    #lut.SetSaturationRange(0.985203081,0.855085320)
    #lut.SetValueRange(0.329415190,0.993247890)
    lut.SetNanColor(1,0,0,1)
    lut.Build()
    cutMapper.SetLookupTable(lut)

    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lut)
    scalarBar.SetTitle(var_name)
    tprop = vtk.vtkTextProperty()
    tprop.SetColor(0,0,0)
    scalarBar.SetTitleTextProperty(tprop)
    scalarBar.SetLabelTextProperty(tprop)
    
    


    # Add the actors to the renderer, set the background and size
    #ren.AddActor(outlineActor)
    #ren.AddActor(planeActor)
      # Create the RenderWindow, Renderer and both Actors
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.SetOffScreenRendering(1) #off-screen render
    renWin.AddRenderer(ren)
    #iren = vtk.vtkRenderWindowInteractor()
    #iren.SetRenderWindow(renWin)
    ren.AddActor(cutActor)
    ren.AddActor2D(scalarBar)

    ren.SetBackground(1, 1, 1)
    renWin.SetSize(1400, 800)

    camera=vtk.vtkCamera()
    camera.SetPosition(plane_normal)
    ren.SetActiveCamera(camera)

    #zoom
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(0.9)
    ren.GetActiveCamera().Dolly(1.4)
    ren.ResetCameraClippingRange()

    #off-screen
    renWin.Render()
    win2im = vtk.vtkWindowToImageFilter()
    win2im.SetInput(renWin)
    win2im.Update()
    writer=vtk.vtkPNGWriter()
    writer.SetFileName(snapshot_dir+'.png')
    writer.SetInputConnection(win2im.GetOutputPort())
    writer.Write()

    #WritePNG(ren,snapshot_dir+'.png', 1)

def vtkcontourf_obj(cutMapper,pl3d_output,src_output,var_name,levels,snapshot_dir,plane_normal,scalarRange,i_d):


#    if not X_is_running():
 #       return

    #off-screen rendering
    #gfac = vtk.vtkGraphicsFactory
    #gfac.SetOffScreenOnlyMode(1)
    #gfac.SetUseMesaClasses(1)
    #im_fac = vtk.vtkImagingFactory
    #im_fac.SetUseMesaClasses(1)   

 
    pl3d_output.GetPointData().SetActiveScalars(var_name)
    if (scalarRange == (0,0)):

        #print scalarRange

        # this prevents zero values from setting colourbar
        scalarnump =  VN.vtk_to_numpy(src_output.GetOutput().GetPointData().GetScalars(var_name))
        minval = np.amin(scalarnump)
        try:
            if (minval >= 0.0):
                minval = np.amin(scalarnump[scalarnump > 0.0])
        except Exception:
            pass
        maxval = np.amax(scalarnump)
        scalarRange = (minval,maxval)
        #print scalarRange


        # set bounds to be symmetric if requested
        if i_d.symmetric:
            if (np.abs(minval) > np.abs(maxval)):
                maxval = -minval
            else:
                minval = -maxval
            scalarRange=(minval,maxval)

        #set minval/maxval to zero if very small
        if (np.abs(minval) < np.finfo(np.float64).eps):
            minval = 0.0
        if (np.abs(maxval) < np.finfo(np.float64).eps):
            minval = 0.0
        scalarRange=(minval,maxval)
       
        #scalarRange =pl3d_output.GetPointData().GetScalars(var_name).GetRange()
    else:
	minval = scalarRange[0]
	maxval = scalarRange[1]

    cutMapper.SetScalarRange(scalarRange)
    cutActor = vtk.vtkActor()
    cutActor.SetMapper(cutMapper)
    #cutActor.SetMapper(src_output)


    lut = vtk.vtkLookupTable() #MakeLUT()
    lut.SetNumberOfTableValues(levels)
    if (i_d.colourmap == 'jet'):
        lut.SetHueRange(0.667, 0.0)
    elif (i_d.colourmap == 'viridis'):
        lut.SetHueRange(0.801282320, 0.149576098)
        lut.SetSaturationRange(0.985203081,0.855085320)
        lut.SetValueRange(0.329415190,0.993247890)
    lut.SetNanColor(1,0,0,1)
    lut.Build()
    cutMapper.SetLookupTable(lut)

    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lut)
    scalarBar.SetTitle(var_name)
    absmax = np.amax([np.abs(minval),np.abs(maxval)])
    if (absmax<1000 and absmax>.1):
        scalarBar.SetLabelFormat("%-#3.2f")
    elif (absmax<.1):
        scalarBar.SetLabelFormat("%-#3.2e")
    tprop = vtk.vtkTextProperty()
    tprop.SetColor(0,0,0)
    #if i_d.non_dimensionalise:
    tprop.SetFontSize(10)
    scalarBar.SetTitleTextProperty(tprop)
    scalarBar.SetLabelTextProperty(tprop)
    scalarBar.SetTextPad(1)
   
    
    


    # Add the actors to the renderer, set the background and size
    #ren.AddActor(outlineActor)
    #ren.AddActor(planeActor)
      # Create the RenderWindow, Renderer and both Actors
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.SetOffScreenRendering(1) #off-screen render
    renWin.AddRenderer(ren)
    #iren = vtk.vtkRenderWindowInteractor()
    #iren.SetRenderWindow(renWin)
    ren.AddActor(cutActor)
    ren.AddActor2D(scalarBar)

    ren.SetBackground(1, 1, 1)
    renWin.SetSize(i_d.render_size[0], i_d.render_size[1])

    camera=vtk.vtkCamera()
    camera.SetPosition(plane_normal)
    ren.SetActiveCamera(camera)

    #zoom
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(0.9)
    ren.GetActiveCamera().Dolly(1.4)
    ren.ResetCameraClippingRange()

    #off-screen
    renWin.Render()
    win2im = vtk.vtkWindowToImageFilter()
    win2im.SetInput(renWin)
    win2im.Update()
    writer=vtk.vtkPNGWriter()
    writer.SetFileName(snapshot_dir+'.png')
    writer.SetInputConnection(win2im.GetOutputPort())
    writer.Write()

    #WritePNG(ren,snapshot_dir+'.png', 1)

def vtkcontourfgrid(cutMapper,pl3d_output,src_output,var_name,levels,snapshot_dir,plane_normal,scalarRange):


#    if not X_is_running():
 #       return

    #off-screen rendering
    #gfac = vtk.vtkGraphicsFactory
    #gfac.SetOffScreenOnlyMode(1)
    #gfac.SetUseMesaClasses(1)
    #im_fac = vtk.vtkImagingFactory
    #im_fac.SetUseMesaClasses(1)   

 
    pl3d_output.GetPointData().SetActiveScalars(var_name)
    if (scalarRange == (0,0)):

        #print scalarRange

        # this prevents zero values from setting colourbar
        scalarnump =  VN.vtk_to_numpy(src_output.GetPointData().GetScalars(var_name))
        minval = np.amin(scalarnump)
        try:
            if (minval >= 0.0):
                minval = np.amin(scalarnump[scalarnump > 0.0])
        except Exception:
            pass
        maxval = np.amax(scalarnump)
        scalarRange = (minval,maxval)
        #print scalarRange


        #scalarRange =pl3d_output.GetPointData().GetScalars(var_name).GetRange()

    cutMapper.SetScalarRange(scalarRange)
    cutActor = vtk.vtkActor()
    cutActor.SetMapper(cutMapper)

    lut = vtk.vtkLookupTable() #MakeLUT()
    lut.SetNumberOfTableValues(levels)
    lut.SetHueRange(0.667, 0.0)
    lut.SetNanColor(1,0,0,1)
    lut.Build()
    cutMapper.SetLookupTable(lut)

    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lut)
    scalarBar.SetTitle(var_name)


    # Add the actors to the renderer, set the background and size
    #ren.AddActor(outlineActor)
    #ren.AddActor(planeActor)
      # Create the RenderWindow, Renderer and both Actors
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.SetOffScreenRendering(1) #off-screen render
    renWin.AddRenderer(ren)
    #iren = vtk.vtkRenderWindowInteractor()
    #iren.SetRenderWindow(renWin)
    ren.AddActor(cutActor)
    ren.AddActor2D(scalarBar)

    ren.SetBackground(0, 0, 0)
    renWin.SetSize(2800, 1600)

    camera=vtk.vtkCamera()
    camera.SetPosition(plane_normal)
    ren.SetActiveCamera(camera)

    #zoom
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(0.9)
    ren.GetActiveCamera().Dolly(1.4)
    ren.ResetCameraClippingRange()

    #off-screen
    renWin.Render()
    win2im = vtk.vtkWindowToImageFilter()
    win2im.SetInput(renWin)
    win2im.Update()
    writer=vtk.vtkPNGWriter()
    writer.SetFileName(snapshot_dir+'.png')
    writer.SetInputConnection(win2im.GetOutputPort())
    writer.Write()

    #WritePNG(ren,snapshot_dir+'.png', 1)


def WritePNG(ren, fn, magnification = 1):
    '''
    Save the image as a PNG
    :param: ren - the renderer.
    :param: fn - the file name.
    :param: magnification - the magnification, usually 1.
    '''
    #renLgeIm = vtk.vtkWindowToImageFilter()
    renLgeIm = vtk.vtkRenderLargeImage()
    renLgeIm.SetInput(ren)
    renLgeIm.SetMagnification(magnification)
    #renLgeIm.SetInputBufferTypeToRGBA()
    #renLgeIm.ReadFrontBufferOff()
    #renLgeIm.Update()
    imgWriter = vtk.vtkPNGWriter()
    imgWriter.SetInputConnection(renLgeIm.GetOutputPort())
    imgWriter.SetFileName(fn)
    imgWriter.Write()

def MakeBands(dR, numberOfBands, nearestInteger):
    '''
    Divide a range into bands
    :param: dR - [min, max] the range that is to be covered by the bands.
    :param: numberOfBands - the number of bands, a positive integer.
    :param: nearestInteger - if True then [floor(min), ceil(max)] is used.
    :return: A List consisting of [min, midpoint, max] for each band.
    '''
    bands = list()
    if (dR[1] < dR[0]) or (numberOfBands <= 0):
        return bands
    x = list(dR)
    if nearestInteger:
        x[0] = math.floor(x[0])
        x[1] = math.ceil(x[1])
    dx = (x[1] - x[0])/float(numberOfBands)
    b = [x[0], x[0] + dx / 2.0, x[0] + dx]
    i = 0
    while i < numberOfBands:
        bands.append(b)
        b = [b[0] + dx, b[1] + dx, b[2] + dx]
        i += 1
    return bands
 
def MakeIntegralBands(dR):
    '''
    Divide a range into integral bands
    :param: dR - [min, max] the range that is to be covered by the bands.
    :return: A List consisting of [min, midpoint, max] for each band.
    '''
    bands = list()
    if (dR[1] < dR[0]):
        return bands
    x = list(dR)
    x[0] = math.floor(x[0])
    x[1] = math.ceil(x[1])
    numberOfBands = int(abs(x[1]) + abs(x[0]))
    return MakeBands(x,numberOfBands, False)
 
def MakeElevations(src):
    '''
    Generate elevations over the surface.
    :param: src - the vtkPolyData source.
    :return: - vtkPolyData source with elevations.
    '''
    bounds = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    src.GetBounds(bounds)
    elevFilter = vtk.vtkElevationFilter()
    elevFilter.SetInputData(src)
    elevFilter.SetLowPoint(0, bounds[2], 0)
    elevFilter.SetHighPoint(0, bounds[3], 0)
    elevFilter.SetScalarRange(bounds[2], bounds[3])
    elevFilter.Update()
    return elevFilter.GetPolyDataOutput()
 
 
def MakePlane():
    '''
    Make a plane as the source.
    :return: vtkPolyData with normal and scalar data.
    '''
    source = vtk.vtkPlaneSource()
    source.SetOrigin(-10.0, -10.0, 0.0)
    source.SetPoint2(-10.0, 10.0, 0.0)
    source.SetPoint1(10.0, -10.0, 0.0)
    source.SetXResolution(20)
    source.SetYResolution(20)
    source.Update()
    return MakeElevations(source.GetOutput())
 
def MakeSphere():
    '''
    Make a sphere as the source.
    :return: vtkPolyData with normal and scalar data.
    '''
    source = vtk.vtkSphereSource()
    source.SetCenter(0.0, 0.0, 0.0)
    source.SetRadius(10.0)
    source.SetThetaResolution(32)
    source.SetPhiResolution(32)
    source.Update()
    return MakeElevations(source.GetOutput())
 
def MakeParametricSource():
    '''
    Make a parametric surface as the source.
    :return: vtkPolyData with normal and scalar data.
    '''
    fn = vtk.vtkParametricRandomHills()
    fn.AllowRandomGenerationOn()
    fn.SetRandomSeed(1)
    fn.SetNumberOfHills(30)
    if fn.GetClassName() == 'vtkParametricRandomHills':
        # Make the normals face out of the surface.
        fn.ClockwiseOrderingOff()
 
    source = vtk.vtkParametricFunctionSource()
    source.SetParametricFunction(fn)
    source.SetUResolution(50)
    source.SetVResolution(50)
    source.SetScalarModeToZ()
    source.Update()
    # Name the arrays (not needed in VTK 6.2+ for vtkParametricFunctionSource)
    source.GetOutput().GetPointData().GetNormals().SetName('Normals')
    source.GetOutput().GetPointData().GetScalars().SetName('Scalars')
    return source.GetOutput()
 
def MakeLUT():
    '''
    Make a lookup table using vtkColorSeries.
    :return: An indexed lookup table.
    '''
    # Make the lookup table.
    colorSeries = vtk.vtkColorSeries()
    # Select a color scheme.
    #colorSeriesEnum = colorSeries.BREWER_DIVERGING_BROWN_BLUE_GREEN_9
    #colorSeriesEnum = colorSeries.BREWER_DIVERGING_SPECTRAL_10
    #colorSeriesEnum = colorSeries.BREWER_DIVERGING_SPECTRAL_3
    #colorSeriesEnum = colorSeries.BREWER_DIVERGING_PURPLE_ORANGE_9
    #colorSeriesEnum = colorSeries.BREWER_SEQUENTIAL_BLUE_PURPLE_9
    #colorSeriesEnum = colorSeries.BREWER_SEQUENTIAL_BLUE_GREEN_9
    colorSeriesEnum = colorSeries.BREWER_QUALITATIVE_SET3
    #colorSeriesEnum = colorSeries.CITRUS
    colorSeries.SetColorScheme(colorSeriesEnum)
    lut = vtk.vtkLookupTable()
    colorSeries.BuildLookupTable(lut)
    lut.SetNanColor(1,0,0,1)
    return lut
 
def ReverseLUT(lut):
    '''
    Create a lookup table with the colors reversed.
    :param: lut - An indexed lookup table.
    :return: The reversed indexed lookup table.
    '''
    lutr = vtk.vtkLookupTable()
    lutr.DeepCopy(lut)
    t = lut.GetNumberOfTableValues() - 1
    revList = reversed(list(range(t + 1)))
    for i in revList:
        rgba = [0,0,0]
        v = float(i)
        lut.GetColor(v,rgba)
        rgba.append(lut.GetOpacity(v))
        lutr.SetTableValue(t - i,rgba)
    t = lut.GetNumberOfAnnotatedValues() - 1
    for i in revList:
        lutr.SetAnnotation(t - i, lut.GetAnnotation(i))
    return lutr
 
def Frequencies(bands, src):
    '''
    Count the number of scalars in each band.
    :param: bands - the bands.
    :param: src - the vtkPolyData source.
    :return: The frequencies of the scalars in each band.
    '''
    freq = dict()
    for i in range(len(bands)):
        freq[i] = 0;
    tuples = src.GetPointData().GetScalars().GetNumberOfTuples()
    for i in range(tuples):
        x = src.GetPointData().GetScalars().GetTuple1(i)
        for j in range(len(bands)):
            if x <= bands[j][2]:
                freq[j] = freq[j] + 1
                break
    return freq
 
def MakeGlyphs(src, reverseNormals):
    '''
    Glyph the normals on the surface.
 
    You may need to adjust the parameters for maskPts, arrow and glyph for a
    nice appearance.
 
    :param: src - the surface to glyph.
    :param: reverseNormals - if True the normals on the surface are reversed.
    :return: The glyph object.
 
    '''
    # Sometimes the contouring algorithm can create a volume whose gradient
    # vector and ordering of polygon (using the right hand rule) are
    # inconsistent. vtkReverseSense cures this problem.
    reverse = vtk.vtkReverseSense()
 
    # Choose a random subset of points.
    maskPts = vtk.vtkMaskPoints()
    maskPts.SetOnRatio(5)
    maskPts.RandomModeOn()
    if reverseNormals:
        reverse.SetInputData(src)
        reverse.ReverseCellsOn()
        reverse.ReverseNormalsOn()
        maskPts.SetInputConnection(reverse.GetOutputPort())
    else:
        maskPts.SetInputData(src)
 
    # Source for the glyph filter
    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(16)
    arrow.SetTipLength(0.3)
    arrow.SetTipRadius(0.1)
 
    glyph = vtk.vtkGlyph3D()
    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.SetInputConnection(maskPts.GetOutputPort())
    glyph.SetVectorModeToUseNormal()
    glyph.SetScaleFactor(1)
    glyph.SetColorModeToColorByVector()
    glyph.SetScaleModeToScaleByVector()
    glyph.OrientOn()
    glyph.Update()
    return glyph
 
def DisplaySurface(st):
    '''
    Make and display the surface.
    :param: st - the surface to display.
    :return The vtkRenderWindowInteractor.
    '''
    surface = st.upper()
    #if  (not(surface in SURFACE_TYPE) ):
    #    print(st, "is not a surface.")
    #    iren = vtk.vtkRenderWindowInteractor()
    #    return iren
    # ------------------------------------------------------------
    # Create the surface, lookup tables, contour filter etc.
    # ------------------------------------------------------------
    src = vtk.vtkPolyData()
    if (surface == "PLANE"):
        src = MakePlane()
    elif (surface == "SPHERE"):
        src = MakeSphere()
    elif (surface == "PARAMETRIC_SURFACE"):
        src = MakeParametricSource()
        # The scalars are named "Scalars"by default
        # in the parametric surfaces, so change the name.
        src.GetPointData().GetScalars().SetName("Elevation");
    scalarRange = src.GetScalarRange()
 
    lut = MakeLUT()
    lut.SetTableRange(scalarRange)
    numberOfBands = lut.GetNumberOfTableValues()
    # bands = MakeIntegralBands(scalarRange)
    bands = MakeBands(scalarRange, numberOfBands, False)
 
    # Let's do a frequency table.
    # The number of scalars in each band.
    #print Frequencies(bands, src)
 
    # We will use the midpoint of the band as the label.
    labels = []
    for i in range(len(bands)):
        labels.append('{:4.2f}'.format(bands[i][1]))
 
    # Annotate
    values = vtk.vtkVariantArray()
    for i in range(len(labels)):
        values.InsertNextValue(vtk.vtkVariant(labels[i]))
    for i in range(values.GetNumberOfTuples()):
        lut.SetAnnotation(i, values.GetValue(i).ToString());
 
    # Create a lookup table with the colors reversed.
    lutr = ReverseLUT(lut)
 
    # Create the contour bands.
    bcf = vtk.vtkBandedPolyDataContourFilter()
    bcf.SetInputData(src)
    # Use either the minimum or maximum value for each band.
    for i in range(0, numberOfBands):
        bcf.SetValue(i, bands[i][2])
    # We will use an indexed lookup table.
    bcf.SetScalarModeToIndex()
    bcf.GenerateContourEdgesOn()
 
    # Generate the glyphs on the original surface.
    glyph = MakeGlyphs(src,False)
 
    # ------------------------------------------------------------
    # Create the mappers and actors
    # ------------------------------------------------------------
    srcMapper = vtk.vtkPolyDataMapper()
    srcMapper.SetInputConnection(bcf.GetOutputPort())
    srcMapper.SetScalarRange(scalarRange)
    srcMapper.SetLookupTable(lut)
    srcMapper.SetScalarModeToUseCellData()
 
    srcActor = vtk.vtkActor()
    srcActor.SetMapper(srcMapper)
    srcActor.RotateX(-45)
    srcActor.RotateZ(45)
 
    # Create contour edges
    edgeMapper = vtk.vtkPolyDataMapper()
    edgeMapper.SetInputData(bcf.GetContourEdgesOutput())
    edgeMapper.SetResolveCoincidentTopologyToPolygonOffset()
 
    edgeActor = vtk.vtkActor()
    edgeActor.SetMapper(edgeMapper)
    edgeActor.GetProperty().SetColor(0, 0, 0)
    edgeActor.RotateX(-45)
    edgeActor.RotateZ(45)
 
    glyphMapper = vtk.vtkPolyDataMapper()
    glyphMapper.SetInputConnection(glyph.GetOutputPort())
    glyphMapper.SetScalarModeToUsePointFieldData()
    glyphMapper.SetColorModeToMapScalars()
    glyphMapper.ScalarVisibilityOn()
    glyphMapper.SelectColorArray('Elevation')
    # Colour by scalars.
    glyphMapper.SetScalarRange(scalarRange)
 
    glyphActor = vtk.vtkActor()
    glyphActor.SetMapper(glyphMapper)
    glyphActor.RotateX(-45)
    glyphActor.RotateZ(45)
 
    # Add a scalar bar.
    scalarBar = vtk.vtkScalarBarActor()
    # scalarBar.SetLookupTable(lut)
    # Use this LUT if you want the highest value at the top.
    scalarBar.SetLookupTable(lutr)
    scalarBar.SetTitle('Elevation (m)')
 
    # ------------------------------------------------------------
    # Create the RenderWindow, Renderer and Interactor
    # ------------------------------------------------------------
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    iren = vtk.vtkRenderWindowInteractor()
 
    renWin.AddRenderer(ren)
    iren.SetRenderWindow(renWin)
 
    # add actors
    ren.AddViewProp(srcActor)
    ren.AddViewProp(edgeActor)
    ren.AddViewProp(glyphActor)
    ren.AddActor2D(scalarBar)
 
    ren.SetBackground(0.7, 0.8, 1.0)
    renWin.SetSize(800, 800)
    renWin.Render()
 
    ren.GetActiveCamera().Zoom(1.5)
 
    return iren
