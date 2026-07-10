SETUP_FLAGS := -d=3
NUM_PROCS ?= 2

.PHONY: all clean setup init convert-topdown-thick convert-corner45-top-thick

all: convert-topdown-thick convert-corner45-top-thick

setup: makefile amrvac
	@mkdir -p output

init: setup
	@$(RM) output/radsyn_sph_prom_0000.dat
	@mpirun -np $(NUM_PROCS) ./amrvac -i init_spherical.par > radsyn_sph_init.log

convert-topdown-thick: init
	@$(RM) output/AIA171_sph_native_topdown_thick_0000.vti
	@mpirun -np $(NUM_PROCS) ./amrvac -i convert_171_sph_native_topdown_thick.par > radsyn_sph_topdown_thick.log

convert-corner45-top-thick: init
	@$(RM) output/AIA171_sph_native_corner45_top_thick_0000.vti
	@mpirun -np $(NUM_PROCS) ./amrvac -i convert_171_sph_native_corner45_top_thick.par > radsyn_sph_corner45_top_thick.log

clean:
	$(RM) amrvac makefile *.f *.mod *.o *.log amrvac.h
	$(RM) -r output

makefile: $(AMRVAC_DIR)/arch/amrvac.make
	@$(RM) $@
	@$(AMRVAC_DIR)/setup.pl $(SETUP_FLAGS) > setup.log

amrvac: makefile
	@$(MAKE)
