# PODFS Digital Filter Inlet Turbulence Generator

This repository contains the PODFS/digital-filter code used to generate turbulent fields at the inlet of scale resolved simulations of turbulent fluid flow. It uses the digital filter method (Klein 2003) to generate the fields and the PODFS method (Treleaven 2019) to compress the data into a useful format that can be efficiently loaded into a CFD solver.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

You will need a working version of Python 2 on your system. We recommend the Anaconda Miniconda system which can be downloaded from [here](https://docs.conda.io/en/latest/miniconda.html).

Once you have installed miniconda and have a new environment running you will require the following prerequisites:

```
Numpy
Matplotlib
Scipy
VTK
h5py
```

### Installing

Simply clone the repository using:


```
git clone https://ncwt3@bitbucket.org/digital_filters_podfs/digital_filters_podfs.git
```

You can run the code by typing:


```
python digitalfilters.py
```

to display the help menu.

## Running the tests

Currently there are no official tests for the code but a minimum working example can be run with:


```
python digitalfilters.py -n 5
```

which will create 5 snapshots of a turbulent flow field and generate the PODFS model of it. You can view the resulting mean field using ParaView which is generated in:


```
PODFS/spatial_mean_field_velocity.vtk
```

The program writes out the following ASCII files:

* PODFS_mean.prf - The mean field with columns X,Y,Z,U,V,W
* PODFS_mode_????.prf - The POD modes wih columns X,Y,Z,u,v,w
* PODFS.dat - The PODFS control file which includes 1. The number of POD modes. 2. The PODFS period. 3. The number of Fourier Coefficients for each POD mode. 4. The Fourier Coefficients.

## Contributing

Please create a new branch for all modifications and use the apporpriate branch naming conventions:

feature/yourfeature - for new features

bugfix/yourbugfix - for new bufixes

Create a pull request when you have a working feature so we can integrate it into our code.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://bitbucket.org/digital_filters_podfs/digital_filters_podfs/downloads/?tab=tags). 

* Version 1.1.0 - For cases that use assumed velocity profile, the mean velocity can now be rotated along with the plane.

* Version 1.0.0 - Initial release.

* Version 0.9.0 - Pre-release version.

## Authors

**Max Staufer** - *Initial work* - Rolls-Royce Deutschland

**Nicholas Treleaven** - *Initial work* - Loughborough University/Rolls-Royce Deutschland/STFS, TU-Darmstadt

**Alessandro Soli** - *Applications to spray break-up simulations* - Loughborough University

For all inquiries, please contact Nick at treleaven@stfs.tu-darmstadt.de

## License

We need a License?

## Acknowledgments

* The initail inspiration for the PODFS method comes from work completed with Laurent Cordier and Laborotoire PPRIME in Poitiers, France.
