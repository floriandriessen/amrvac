TEST_DIRS := RHD_heatcool

##TEST_DIRS := RHD_wave_1D RHD_shock_1D RHD_pulse_1D RHD_heatcool

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: all clean $(TEST_DIRS) $(CLEAN_TEST)

all: $(TEST_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) clean
