#
#   plot_contours
#                               
#   Script to plot a converged flowfield from the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_contours.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filenames to read the run data
    inname = 'input_' + sys.argv[-1] + '.txt'
    outname = 'out_final_' + sys.argv[-1] + '.bin'

    # Read the settings and the case from file
    av = read_settings(inname)
    g = read_case(outname)

    # When presenting results all values should be non-dimensionalised. Two
    # variables of interest might be:
    #    1. Static pressure coefficient, (p - p_ref) / (pstag_ref - p_ref)
    #    2. Mach number, v / (ga * rgas * t)**0.5

    # First complete the "calc_secondary" function within "routines.py" to
    # calculate static pressure and Mach number, and any others you want!
    g = calc_secondary(av,g)

    # Loop over every i-station, mass average T0 and sum the mass flow
    ni = g['ni']
    mass_i = np.zeros(ni)
    tstag_i = np.zeros(ni)
    for i in range(ni):
        c = cut_i(g,i)
        tstag_i[i], mass = mass_av(c,'tstag')
        mass_i[i] = np.sum(mass)

    # Normalise
    mass_ratio = mass_i / mass_i[0]
    tratio = tstag_i / av['tstag']

    vals_list = [tratio, mass_ratio]
    ylabels = ['Stagnation temperature ratio',
               r'Mass flow ratio']
    i_index = np.arange(1, ni + 1)

    for vals, label in zip(vals_list, ylabels):
        fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes()
        ax.plot(i_index, vals, marker='o', markersize=3, linewidth=1)
        ax.set_xlabel('Cell index')
        ax.set_ylabel(label)
        ax.grid(True)

    plt.show()


    
main()


