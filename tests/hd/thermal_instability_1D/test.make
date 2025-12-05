SETUP_FLAGS := -d=1 -v=1

SCHEME_DIR := ../../schemes

SCHEMES := 2step_tvdmu_al 3step_hll_cada 4step_hll_mc	\
4step_hllc_ko rk4_tvdlf_cada stsmethod1 stsmethod2

BASE_NAME := TIautotest_wc

PAR_FILE := TIautotest.par

TESTS := $(SCHEMES:%=$(BASE_NAME)_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests, which are used to run them
$(foreach s, $(SCHEMES),\
    $(eval $(s:%=$(BASE_NAME)_%.log): $(PAR_FILE) $(SCHEME_DIR)/$(s).par))
