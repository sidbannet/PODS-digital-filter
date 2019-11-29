import digitalfilter_oop as df
from optparse import OptionParser
from optparse import Option, OptionValueError
import sys

PROG = 'DigitalFilters'
VERSION = '2.0'


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

description = """ LES Inflow Generator after Klein et.al. """

parser = OptionParser(option_class=MultipleOption,
                      usage='usage: %prog [options]',
                      version='%s %s' % (PROG, VERSION),
                      description=description)

parser.add_option("-i", "--inputfile", dest="profilefile", default='none',
                  help="1d or 2d turbulent profile file", metavar="FILE")

parser.add_option("-o", "--outputfolder", dest="outputfolder", default='cwd',
                  help="1d or 2d turbulent profile file", metavar="FILE")

parser.add_option("-p", "--mean_profile", dest="mean_profile", default='hyperbolic-tangent',
                  help="What kind of mean flow profile would you like? Options: hyperbolic-tangent, \
                  double-hyperbolic-tangent, ring-hyperbolic-tangent, circular-hyperbolic-tangent.  \
                  Using the -i option along with this will adapt the -i profile to the shape of the -p profile \
                  ie. plane jet, square jet, round jet, annulus jet.",
                  metavar="STRING")

parser.add_option("--turb_profile", dest="turb_profile", default='top-hat',
                  help="What kind of turbulence flow profile would you like? Options: top-hat, none",
                  metavar="STRING")

parser.add_option("--U0", "--bulk_velocity", type="float", dest="bulk_velocity", default=0.0,
                  help="What is the bulk velocity magnitude?, this \
			      option can also be used to scale the velocities \
			      and turbulent quantities similar to the massflow \
			      option is using a .prf profile.", metavar="NUM")

parser.add_option("--u_dash", type="float", dest="turbulence_intensity", default='0.02',
                  help="What is the desired u'/U0 value with u'=v'=w'?", metavar="NUM")

parser.add_option("-n", "--nsteps", type="int", dest="nsteps", default=20,
                  help="number of steps", metavar="INT")

parser.add_option("-l", "--lengthscale", type="float", dest="lengthscale", default=3.0,
                  help="turbulent lengthscale in terms of grid spacing ", metavar="NUM")

parser.add_option("-f", "--fwidth", type="float", dest="fwidth", default=2.0,
                  help="half filter width in lengthscales, should be greater than 2", metavar="NUM")

parser.add_option("-k", "--nk", type="int", dest="kma", default=11,
                  help="number of points in k (wall-normal) direction", metavar="INT")

parser.add_option("-j", "--nj", type="int", dest="jma", default=10,
                  help="number of points in j (spanwise) direction", metavar="INT")

parser.add_option("-t", "--dt", type="float", dest="dt", default=0.0,
                  help="time step (s)", metavar="NUM")

parser.add_option("-s", "--dtsplit", type="float", dest="dt_split", default=1.0,
                  help="number of timestep splits for interploation between filtered timesteps", metavar="NUM")

parser.add_option("--periodic", type="float", dest="periodic", default=False,
                  help="decide whether create periodic signal or not", metavar="store_true")

parser.add_option("-m", "--nm", type="int", dest="nm", default=20,
                  help="number of POD modes", metavar="INT")

parser.add_option("-e", "--et", type="float", dest="et", default=0.9,
                  help="target energy for Fourier reconstruction", metavar="NUM")

parser.add_option("-v", "--verbose", dest="verbose", default=False,
                  help="Save the mean flow, POD spatial and temporal modes?", action='store_true')

parser.add_option("--non_dim", dest="non_dim", default=False,
                  help="Non-dimensionalise lengths if using .prf", action='store_true')

parser.add_option("-r", "--resolution", type="float", dest="res", default=0.1,
                  help="plane resolution in meters per grid point", metavar="NUM")

parser.add_option("--nx", type="float", dest="nx", default=1.0,
                  help="plane normal direction, x-component", metavar="NUM")

parser.add_option("--ny", type="float", dest="ny", default=0.0,
                  help="plane normal direction, y-component", metavar="NUM")

parser.add_option("--nz", type="float", dest="nz", default=0.0,
                  help="plane normal direction, z-component", metavar="NUM")

parser.add_option("--ox", type="float", dest="ox", default=0.0,
                  help="plane origin, x-component", metavar="NUM")

parser.add_option("--oy", type="float", dest="oy", default=0.0,
                  help="plane origin, y-component", metavar="NUM")

parser.add_option("--oz", type="float", dest="oz", default=0.0,
                  help="plane origin, z-component", metavar="NUM")

parser.add_option("--rotate", type="float", dest="rot", default=0.0,
                  help="rotate plane about its normal (degrees)", metavar="NUM")

parser.add_option("--ring", type="float", dest="ring", default=0.5,
                  help="inner diameter of ring as a proportion of the outer \
			diameter if using the ring-hyperbolic-tangent option", metavar="NUM")

parser.add_option("--massflow", type="float", dest="mdot", default=0.0,
                  help="If using a .prf file, the velocities can be  \
			scaled to achieve the desired mass flow rate. Mean \
			velocities are scaled equally. k and epsilon/sdr \
			are scaled to maintain the same turbulence intensity \
			and length scale. Requires density!", metavar="NUM")

parser.add_option("--density", type="float", dest="den", default=0.0,
                  help="If massflow is specified, a density must  \
                        also be specified.", metavar="NUM")

if len(sys.argv) == 1:
    parser.parse_args(['--help'])

(options, args) = parser.parse_args()

# settings for plot output
plot_general = False  # general plot flag
plot_mean = True  # plot mean fields
plot_rand = True  # plot random fields
plot_corr = True  # plot correlated fields
plot_final = True  # plot final fields
plot_flag = [plot_general, plot_mean, plot_rand, plot_corr, plot_final]

# create object filter
digfilter = df.Filter(filename=options.profilefile,
                      res=options.res,
                      outdir=options.outputfolder,
                      dt_split=options.dt_split,
                      bulk_velocity=options.bulk_velocity,
                      plot_flag=plot_flag,
                      density=options.den,
                      mdot=options.mdot,
                      profile=options.mean_profile,
                      fwidth=options.fwidth,
                      verbose=options.verbose,
                      # keep_all=options.keep_all,
                      res_y=options.jma,
                      res_z=options.kma,
                      oy=options.oy,
                      oz=options.oz,
                      turb_prof=options.turb_profile,
                      tu_intensity=options.turbulence_intensity,
                      lengthscale=options.lengthscale,
                      inner_diameter_norm=options.ring)
print('Filterzeitschritt: {:8e} | Zeitschritt Inletdateien:  {:8e}'.format(digfilter.dt_filter, digfilter.dt_filter /
                                                                           digfilter.dt_split))

# apply filter and create snapshots
digfilter.filtering(options.dt, options.periodic)
