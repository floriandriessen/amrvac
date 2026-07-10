SETUP_FLAGS := -d=3
NUM_PROCS ?= 2

.PHONY: all clean setup init convert-top-thick convert-oblique-thick

all: convert-top-thick convert-oblique-thick

setup: makefile amrvac
	@mkdir -p output

init: setup
	@$(RM) output/radsyn_tdm_prom_amr_0000.dat
	@mpirun -np $(NUM_PROCS) ./amrvac -i init_amr.par > radsyn_init.log

convert-top-thick: init
	@$(RM) output/AIA171_thick_0000.vti
	@mpirun -np $(NUM_PROCS) ./amrvac -i convert_171_thick.par > radsyn_thick.log

convert-oblique-thick: init
	@$(RM) output/AIA171_cart_dda_thick_oblique_0000.vti
	@mpirun -np $(NUM_PROCS) ./amrvac -i convert_171_cart_dda_thick_oblique.par > radsyn_cart_dda_thick_oblique.log

clean:
	$(RM) amrvac makefile *.f *.mod *.o *.log amrvac.h
	$(RM) -r output

makefile: $(AMRVAC_DIR)/arch/amrvac.make
	@$(RM) $@
	@$(AMRVAC_DIR)/setup.pl $(SETUP_FLAGS) > setup.log

amrvac: makefile
	@$(MAKE)
