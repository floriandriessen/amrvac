SETUP_FLAGS := -d=2 -v=2
TESTS := ext_bruselator_2d.log 

include ../../test_rules.make

ext_bruselator_2d.log: spirals.par

