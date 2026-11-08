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

    # Use the "cut_i", "mass_av" AND "area_av" functions to calculate the
    # reference pressures at the inlet plane and therefore the static pressure
    # coefficient
    # INSERT
        # Mach number under the name used by the plotting loop
    g['mach'] = g['M']

    # Local stagnation pressure at every node (varies in the flow field, so
    # it must be computed before it can be mass averaged)
    g['pstag'] = g['p'] * (1 + 0.5 * (av['gam'] - 1) * g['mach']**2) \
        ** (av['gam'] / (av['gam'] - 1))

    # Take the inlet plane (i = 0 in Python indexing) and calculate the
    # reference pressures. area_av and mass_av return (value, weights), so
    # keep only the first output.
    c_in = cut_i(g,0)
    p_ref,_ = area_av(c_in,'p')
    pstag_ref,_ = mass_av(c_in,'pstag')

    # Static pressure coefficient, stored in g so the plotting loop can find it
    g['cp'] = (g['p'] - p_ref) / (pstag_ref - p_ref)

    # Specify the parameters to plot
    fieldnames = ['cp', 'mach']; 
    colnames = ['Static pressure coefficient','Mach number']

    # Plot the calculated non-dimensional parameters to show the flow solution
    for n,name in enumerate(fieldnames):

        # Open figure window
        fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    
        # Set aspect ratio as equal and remove axes labels
        ax.set_aspect('equal',adjustable='box'); ax.axis('off')
 
        # Plot filled contour levels
        hc = ax.pcolormesh(g['x'],g['y'],g[name],shading='gouraud')

        # Add colorbar with variable name
        colorbar(hc,colnames[n])

        # Add Mach = 1 contours
        if name == 'mach':
            ax.contour(g['x'],g['y'],g['mach'],[1.0],colors='w',
                linewidths=0.5) 

        # Draw the walls of the block
        plot_wall(ax,g)

    # Show all the plots
    plt.show()

    
main()


