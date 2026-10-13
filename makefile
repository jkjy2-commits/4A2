# This file manages the compilation of the solver code. 
#
# All lines starting with "#" are comment lines not used

SaveName = $(USER)SaveSrc
SaveName = SaveSrc

OBS = types.o routines.o flux_stencil.o smooth_stencil.o solver.o \
    read_settings.o allocate_arrays.o read_geom.o generate_mesh.o read_mesh.o \
    calc_areas.o write_output.o check_mesh.o flow_guess.o set_timestep.o \
    set_secondary.o apply_bconds.o euler_iteration.o check_conv.o check_stop.o

# Choice of comipler
FC = gfortran

# Option of Compilation

# Debug option
FFLAGS = -fcheck=all -g

# Select the optimising option by deleting # in the following line
#FFLAGS = -O2

# Compile the executable and name it "solver.x" using "-o"
solver : $(OBS)
	$(FC)  $(FFLAGS) -o solver.x  $(OBS)  

# remove all .o files and Euler for a new start
clean:
	rm -fr *.o solver.x Save*

# save *.f90 and *.py code in directory SaveSrc and compress it
save:
	if [ ! -d $(SaveName) ]; then mkdir $(SaveName); fi
	cp *.f90 $(SaveName)/
	cp *.py $(SaveName)/
	cp makefile $(SaveName)/
	tar -cf $(SaveName).tar $(SaveName)
	gzip $(SaveName).tar
	rm -r -f $(SaveName)

# Extract files from SaveSrc.tar.gz and save them to SaveSrc directory
extract:
	gunzip $(SaveName).tar.gz
	tar -xf $(SaveName).tar

# Compile .f90 to .o
%.o:%.f90
	$(FC) $(FFLAGS) -c $^ -o $@
