FOBJECTS += mod_eos_container.t mod_eos_interp.t mod_eos_saha.t mod_eos_entropy.t mod_eos_tables.t mod_eos.t

# Explicit dependencies (otherwise parallel make can race on .mod files)
mod_eos_interp.o: mod_eos_container.mod
mod_eos_saha.o: mod_eos_container.mod mod_eos_interp.mod
mod_eos_entropy.o: mod_eos_container.mod
mod_eos_tables.o: mod_eos_container.mod mod_eos_interp.mod
mod_eos.o: mod_eos_container.mod mod_eos_interp.mod mod_eos_saha.mod mod_eos_entropy.mod mod_eos_tables.mod
