Quickstart Guide
================



Prerequisites
#############

You will need a working version of Python 2 on your system. We recommend the Anaconda Miniconda system which can be downloaded from `here <https://docs.conda.io/en/latest/miniconda.html>`_.


Once you have installed miniconda and have a new environment running you will require the following prerequisites::


	Numpy
	Matplotlib
	Scipy
	VTK
	h5py
	sphinx


If you have conda installed, you can install the prerequisites with::

	conda install numpy vtk scipy vtk h5py sphinx

Installation
############

Simply clone the repository using::

	git clone https://ncwt3@bitbucket.org/digital_filters_podfs/digital_filters_podfs.git


You can run the code by typing::

	python digitalfilters.py

to display the help menu.

Running the tests
#################

Currently there are no official tests for the code but a minimum working example can be run with::

	python digitalfilters.py -n 5


which will create 5 snapshots of a turbulent flow field and generate the PODFS model of it. You can view the resulting mean field using ParaView which is generated in::



	PODFS/spatial_mean_field_velocity.vtk


The program writes out the following ASCII files:

* PODFS_mean.prf - The mean field with columns X,Y,Z,U,V,W
* PODFS_mode_????.prf - The POD modes wih columns X,Y,Z,u,v,w
* PODFS.dat - The PODFS control file which includes
	1. The number of POD modes. 
	2. The PODFS period in seconds. 
	3. The number of Fourier Coefficients for each POD mode. 
	4. The Fourier Coefficients.


Contributing
############

Please create a new branch for all modifications and use the apporpriate branch naming conventions:

feature/yourfeature - for new features

bugfix/yourbugfix - for new bufixes

Create a pull request when you have a working feature so we can integrate it into our code.

**PLEASE** remember to write documentation for your features using `sphinx <https://www.sphinx-doc.org/en/master/>`_

Versioning
##########

We use `SemVer <http://semver.org/>`_ for versioning. For the versions available, see the `tags on this repository <https://bitbucket.org/digital_filters_podfs/digital_filters_podfs/downloads/?tab=tags>`_. 

* Version 1.1.0 - For cases that use assumed velocity profile, the mean velocity can now be rotated along with the plane.

* Version 1.0.0 - Initial release.

* Version 0.9.0 - Pre-release version.

Authors
#######

**Max Staufer** - *Initial work* - Rolls-Royce Deutschland

**Nicholas Treleaven** - *Project Lead* - Loughborough University/Rolls-Royce Deutschland/STFS, TU-Darmstadt

**Alessandro Soli** - *Applications to compressor/combustor coupling* - Loughborough University

**Jonathan Gruendler** - *Applications to combustor/turbine coupling and Python 3* - GLR, TU-Darmstadt 

**Vishal Saini** - *Radial combustor inlet and performance improvements* - Loughborough University

For all inquiries, please contact Nick at treleaven@stfs.tu-darmstadt.de

License
#######

* `BSD-3-Clause <https://opensource.org/licenses/BSD-3-Clause>`_
* Copyright 2020 Nicholas Treleaven

Publications
############

The methods used in this code are explained in:

* Treleaven, N. C. W., et al. "Application of the PODFS method to inlet turbulence generated using the digital filter technique." Journal of Computational Physics (2020): 109541. `https://doi.org/10.1016/j.jcp.2020.109541 <https://doi.org/10.1016/j.jcp.2020.109541>`_. 

If you using the PODFS method in your work please reference:

* Treleaven, Nicholas CW, et al. "An efficient method to reproduce the effects of acoustic forcing on gas turbine fuel injectors in incompressible simulations." Flow, Turbulence and Combustion 103.2 (2019): 417-437. `https://doi.org/10.1007/s10494-019-00020-4 <https://doi.org/10.1007/s10494-019-00020-4>`_.

The digital filter method as used is first explained here:

* Klein, Markus, Amsini Sadiki, and Johannes Janicka. "A digital filter based generation of inflow data for spatially developing direct numerical or large eddy simulations." Journal of computational Physics 186.2 (2003): 652-665. `https://doi.org/10.1016/S0021-9991(03)00090-1 <https://doi.org/10.1016/S0021-9991(03)00090-1>`_.




Acknowledgments
###############

* The initail inspiration for the PODFS method comes from work completed with Laurent Cordier and Laborotoire PPRIME in Poitiers, France.




