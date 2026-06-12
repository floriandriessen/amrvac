FOBJECTS += mod_eos_container.t mod_eos_shared_functions.t mod_eos_interp.t mod_eos_LTE_tables.t mod_eos_LTE_saha.t mod_eos_LTE_entropy.t mod_eos_LTE_state.t mod_eos_PI_tables.t mod_eos_FI.t mod_eos_PI.t mod_eos_LTE.t mod_eos.t

# Explicit dependencies (otherwise parallel make can race on .mod files)
mod_eos_shared_functions.o: mod_eos_container.mod
mod_eos_interp.o: mod_eos_container.mod
mod_eos_LTE_tables.o: mod_eos_container.mod mod_eos_interp.mod
mod_eos_LTE_saha.o: mod_eos_container.mod mod_eos_interp.mod mod_eos_LTE_tables.mod
mod_eos_LTE_entropy.o: mod_eos_container.mod mod_eos_LTE_tables.mod
mod_eos_PI_tables.o: mod_comm_lib.mod mod_constants.mod mod_global_parameters.mod
mod_eos_FI.o: mod_eos_container.mod
mod_eos_PI.o: mod_eos_container.mod mod_eos_PI_tables.mod mod_eos_FI.mod
mod_eos_LTE.o: mod_eos_container.mod mod_eos_shared_functions.mod mod_eos_interp.mod mod_eos_LTE_saha.mod mod_eos_LTE_entropy.mod mod_eos_LTE_tables.mod
mod_eos.o: mod_eos_container.mod mod_eos_shared_functions.mod mod_eos_interp.mod mod_eos_LTE_saha.mod mod_eos_LTE_entropy.mod mod_eos_LTE_tables.mod mod_eos_FI.mod mod_eos_PI.mod mod_eos_LTE.mod
