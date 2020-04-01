import digitalfilter_oop as df
import PODFS as podfs

# main input Digital Filter
profilname = '../Input/inlet_centaur_grid.csv'                                  # input profile
# profilname = '../Input/profile_3_new.prf'
resolution = 0.001                                                              # resolution of filter plane
# resolution = 0.002892995
outputfolder = '../Output/test_polar/HYDRA'                                    # path to output folder
# outputfolder = '../Output/777/'
dt_split = 1                                                                    # split filter time step and interpolate
output_timesteps = 200                                                           # number of time steps
verbose = False
periodic = True                                                                 # use first set of fields for last set
polar = True                                                                    # For Turbine inlet: use polar coordinates
prof2d = True                                                                   # Use 2D profile as input
prof1d = False                                                                   # Use 1D profile as input
hdf5 = True                                                                     # save snapshots in hdf5 file (for HYDRA)
fluct_comps = ['u', 'v', 'w', 'Tt']                                             # variables to create fluctuations for
temp_fluct = 200                                                                # standard deviation of temperature field for the calculation of temperature fluctuations
PODFS = True                                                                    # PODFS methode

# special input scaling
mdot = 0                                                                        # target mass flow for scaling
bulk_velocity = 0                                                               # target bulk velocity for scaling
density = 1.71066                                                               # density for massflow scaling

# special input profile creation
profile = 'hyperbolic-tangent'                                                  # shape of mean velocity profile
turb_prof = 'top-hat'                                                           # shape of turbulence profile
tu_intensity = 0.06                                                             # turbulence intensity
res_y = 20                                                                      # Number of grid points in y direction
res_z = 40                                                                      # Number of grid points in z direction
oy = 0                                                                          # coordinate origin y
oz = 0                                                                          # coordinate origin z
lengthscale = 3                                                                 # length scale in terms of grid spacing
y_range = (0.0102242, 0.2702242)                                                # dimension y direction
z_range = (-0.0232496, 0.0232496)                                               # dimension z direction

# settings for plot output
plot_general = False                                                             # general plot flag
plot_mean = True                                                                # plot mean fields
plot_rand = False                                                                # plot random fields
plot_corr = False                                                                # plot correlated fields
plot_final = False                                                               # plot final fields
plot_flag = [plot_general, plot_mean, plot_rand, plot_corr, plot_final]

# create object filter
digfilter = df.Filter(filename=profilname,
                      res=resolution,
                      outdir=outputfolder,
                      dt_split=dt_split,
                      bulk_velocity=bulk_velocity,
                      plot_flag=plot_flag,
                      density=density,
                      mdot=mdot,
                      polar=polar,
                      prof2d=prof2d,
                      prof1d=prof1d,
                      profile=profile,
                      turb_prof=turb_prof,
                      tu_intensity=tu_intensity,
                      res_y=res_y, res_z=res_z, oy=oy, oz=oz,
                      lengthscale=lengthscale,
                      y_range=y_range, z_range=z_range,
                      fluct_comps=fluct_comps,
                      temp_fluct=temp_fluct,
                      verbose=verbose)
print('Filterzeitschritt: {:8e} | Zeitschritt Inletdateien:  {:8e}'.format(digfilter.dt_filter, digfilter.dt_filter /
                                                                           digfilter.dt_split))

# apply filter and create snapshots
digfilter.filtering(output_timesteps, periodic, hdf5, PODFS)

# main input PODFS

num_modes = 20                                                                  # number of modes to write
correct_for_cell_volumes = False
energy_target = 0.9                                                             # target of procentual energy resolved by fourier series

# settings for plot output
plot_temp_modes = True

if PODFS:
    pod = podfs.POD(digfilter.prepare_POD, num_modes, outputfolder, correct_for_cell_volumes)

    fourier = podfs.FS(digfilter.prepare_POD, pod, outputfolder, energy_target, plot_temp_modes)

    save_hdf5 = podfs.HDF5(pod, fourier)
