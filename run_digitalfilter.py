import digitalfilter_oop as df

# main input
profilname = '../Input/profile_3_new.prf'    # input profile
resolution = 0.001                                                              # resolution of filter plane
outputfolder = '../Output/test'                                                     # path to output folder
dt_split = 1                                                                    # split filter time step and interpolate
output_timesteps = 20                                                           # number of time steps
periodic = True                                                                 # use first set of fields for last set

# special input
mdot = 1                                                                        # target mass flow for scaling
bulk_velocity = 1                                                               # target bulk velocity for scaling
density = 1.5                                                                   # density

# settings for plot output
plot_general = False                                                             # general plot flag
plot_mean = True                                                                # plot mean fields
plot_rand = True                                                                # plot random fields
plot_corr = True                                                                # plot correlated fields
plot_final = True                                                               # plot final fields
plot_flag = [plot_general, plot_mean, plot_rand, plot_corr, plot_final]

# create object filter
digfilter = df.Filter(filename=profilname,
                      res=resolution,
                      outdir=outputfolder,
                      dt_split=dt_split,
                      bulk_velocity=bulk_velocity,
                      plot_flag=plot_flag,
                      density=density,
                      mdot=mdot)
print('Filterzeitschritt: {:8e} | Zeitschritt Inletdateien:  {:8e}'.format(digfilter.dt_filter, digfilter.dt_filter /
                                                                           digfilter.dt_split))

# apply filter and create snapshots
digfilter.filtering(output_timesteps, periodic)
