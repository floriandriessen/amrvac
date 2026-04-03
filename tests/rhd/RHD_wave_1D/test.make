SETUP_FLAGS := -d=1 -v=1
SCHEME_DIR := ../../schemes

SCHEMES := IMEX_Euler_hll_w5  IMEX_SP_hll_ko \
 IMEX_SP_hll_mm IMEX_SP_tvdlf_ko  \
 IMEX_222_hll_w5 IMEX_Trap_hll_w5 IMEX_ARS3_hll_cada3 IMEX_232_hll_cada3 IMEX_SP_hll_w5

SCHEMES_FIXED := 

SCHEMES_BLOW := IMEX_Midp_hll_w5

TESTS := $(SCHEMES:%=rhd_wave_%.log) $(SCHEMES_FIXED:%=rhd_wave_%_fixes_on.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=rhd_wave_%.log): rhd_wave.par $(SCHEME_DIR)/$(s).par))

$(foreach s, $(SCHEMES_FIXED),\
	$(eval $(s:%=rhd_wave_%_fixes_on.log): rhd_wave.par $(SCHEME_DIR)/$(s).par $(SCHEME_DIR)/fixes_on.par))
