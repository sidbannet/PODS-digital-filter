
Options
=======

The code is called using::

	python digitalfilters.py

Failure to select any options will activate the help menu which is also accessible using the *-h* or *--help* options tags.

The program once activated will generate a set of 2D planar turbulent fields with the properties specified and then compute the PODFS approximation to those fields. The following options may be selected:

Input file
##########

The desired velocity and Reynolds stress profile can be defined in an input file and loaded using::

	-i [inputfile]
or::

	--inputfile=[inputfile]

This specifies the desired ASCII input profile file. It is also possible to choose from default profiles  (see *--mean\_profile*). The file should contain a series of columns titled `y', the wall normal direction, `U', the mean axial velocity, `uu', the axial component of Reynolds stress, `vv', the wall normal component of Reynolds stress, `ww', the spanwise component of Reynolds stress and `uv', the Reynolds shear stress along the wall. The range of y should be [0,1] where 0 is the wall and 1 is the centreline of the channel. If the input file is of the type `.prf' then a 2D profile will be generated from the input file. In this case the file should have the columns `x',`y',`z',`u',`v',`w' and either `te` (turbulent kinetic energy), and one of `ed' (dissipation rate) or `sdr' (omega), or the Reynolds stress profile `uu', `vv', `ww' (cross stresses are optional). The file should include a header that indicates the order of the columns in the format *data,x,y,z,u,v,w,te...*. If Reynolds stresses are excluded, they are generated using the Boussinesq eddy viscosity approximation. The length scale is assumed to be *0.07Dh* where *Dh* is the hydraulic diameter. But it may be overridden using the *-l* command.


Mean flow profile
#################

The mean flow profile is selected with the option::

	-p [profile_type]
or::

	--mean_profile=[profile_type] 

This option allows the user to select from a range of standard profile shapes. If this option is combined with the *-i* option then the loaded 1D profile will be altered to fit the *[profile_type]* chosen. *default=hyperbolic-tangent*. The options are:



Hyperbolic-tangent
******************

A hyperbolic tangent profile as used in a planar jet.

.. figure:: plane.png
    :width: 496px
    :align: center
    :height: 382px
    :alt: alternate text


    The hyperbolic-tangent mean flow profile

This profile is selected using::

	-p hyperbolic-tangent

or:: 

        --mean_profile=hyperbolic-tangent
 

Double hyperbolic-tangent
*************************

As above but the profile is applied in both principle directions as would be the case for a square jet.


.. figure:: square.png
    :width: 496px
    :align: center
    :height: 382px
    :alt: alternate text


    The double hyperbolic-tangent mean flow profile


This profile is selected using::

        -p double-hyperbolic-tangent

or::

        --mean_profile=double-hyperbolic-tangent


Circular hyperbolic-tangent
***************************

A hyperbolic tangent profile as used in a round jet.

.. figure:: circular.png
    :width: 496px
    :align: center
    :height: 382px
    :alt: alternate text


    The circular hyperbolic-tangent mean flow profile


This profile is selected using::

        -p circular-hyperbolic-tangent

or::

        --mean_profile=circular-hyperbolic-tangent



Ring hyperbolic-tangent
***********************

A hyperbolic tangent profile as used in a ring such as the outermost inlets of a fuel spray nozzle or the inlet of the entire annulus of a jet engine compressor, fan or turbine. 

.. figure:: ring.png
    :width: 496px
    :align: center
    :height: 382px
    :alt: alternate text


    The ring hyperbolic-tangent mean flow profile


This profile is selected using::

        -p ring-hyperbolic-tangent

or::

        --mean_profile=ring-hyperbolic-tangent

This feature should also be used in conjunction with the command::


	--ring=0.5 

which corresponds to the inner diameter of the ring expressed as a proportion of the outer diameter. *Default = 0.5*

Turbulence profile
##################

This option allows the user to select from a range of standard profile shapes (*default=top-hat*) , they are: 

Top-hat
*******

uu=vv=ww over all the inlet. To select this inlet type use::

	--turb_profile=top-hat

None
****

No turbulence, use::

	--turb_profile=none

Bulk velocity
#############

The bulk velocity magnitude in the plane normal direction. *Default = 1.0*

Use::

	--U0=1.0

or::

	--bulk_velocity=1.0 

Turbulence intensity
####################

The ratio of u/U0 can be set with::

	--udash=0.02

It is currently assumed that. u=v=w. *Default = 0.02*

Time stepping
#############

The number of timesteps calculated using the digital filter technique is set with::

	-n 20
or::

	--num_steps=20 

This does not have to be equal to the number of steps in the final simulation but the number of steps multiplied by the time step will give the repeating period of the inlet. *Default = 20*

The time step (in seconds) is set using::

	-t 

or:: 

	--dt=

If this is not specified the time step will be calculated based on the length scale and the mean velocity. If this is specified then the axial length scale will be adjusted based on the time step and mean velocity. The length scales in the j and k directions will not change.

Turbulent length scale
######################

The turbulent length scale in terms of the inlet plane grid spacing (*Default = 3.0*) is set using::

	-l 3.0 

or::
	--lengthscale=3.0

Filter width
############

The half filter width in terms of the number of length scales. This should be greater than 2. *Default = 2.0* 

This is set using::

	-f 

or::

--filter_width=2.0 


Inlet mesh position, rotation, size and resolution
##################################################

The resolution and size of the inlet plane is set with::

	-k 11 

or::

	--nk=11 

Which corresponds to the number of grid points in the k or wall-normal direction. Each cell is assumed to be square, therefore the shape of the inlet can be stretched by adjusting the number of cells in each direction. *Default = 11*

and::

	-j 10 

or:: 
	--nj=10 

which is the number of grid points in the j or span-wise direction. Each cell is assumed to be square, therefore the shape of the inlet can be stretched by adjusting the number of cells in each direction. *Default = 10*

The resolution of the inlet (meters per grid point), is set using::

	-r

or:: 

	--resolution=0.1

*Default = 0.1 m*

The position of the inlet plane is defined using::

	--nx=1.0 

The inlet plane normal direction, x-component. *Default = 1*

and::

	--ny=0.0 

The inlet plane normal direction, y-component. *Default = 0*

and::

	--nz=0.0 

The inlet plane normal direction, z-component. *Default = 0*

and::

	--ox=0.0

The x-coordinate of the centre of the inlet plane. *Default = 0*


and::

	--oy=0.0

The y-coordinate of the centre of the inlet plane. *Default = 0*

and::

	--oz=0.0 

The z-coordinate of the centre of the inlet plane. *Default = 0*


You can rotate the plane about its normal vector using::

	--rotate=0.0 

The entered value should be in degrees. You can use the *POD_spatial_mean_field_velocity.vtk* file to check if the plane is in the correct location. *Default = 0*

PODFS
#####

You can control the PODFS modelling using::

	-m 20

or::

	--num_modes=20

which is the number of POD modes to use in the PODFS model. This should be less than the number of timesteps. A small number will produce a less accurate but faster model, a large number produces a more accurate but slower model. Use the *POD.eigenvalues.dat* file to check how much of the fluctuating energy field is captured using the chosen number of modes. *Default = 20*

and::

	-e 0.9 

or::

	--et=0.9 

which is the energy target for the Fourier reconstruction. 0.9=90\%. A small number will produce a less accurate but faster model, a large number produces a more accurate but slower model. If numbers very close to 1 are used there may be spurious high frequency oscillations in the temporal reconstruction. Use the *-v* option to check the *POD_tmode_recon?.png* files using a small number of *-n*. *Default = 0.9*


Verbosity
#########

The flag::

	-v
or:: 

	--verbose

may be used to save the generated digital filter fields used to produce the PODFS model as .prf files which can be loaded directly into your CFD solver. It also saves the POD temporal coefficients and the Fourier reconstructions of them. The spatial POD modes are also saved in .vtk format.

