# Minimal `mod_usr.t` Hook for Topology/QSL Namelist Tasks

This is a snippet, not a complete `mod_usr.t` file. Add the same pattern to
the problem-specific `mod_usr.t` that already initializes the physics and
magnetic field.

The first namelist runner does not add a dedicated global `convert_type`.
Use AMRVAC's user conversion hook, typically with `convert_type='user'`.
If your local problem uses a different hook name, adapt the pointer assignment
accordingly.

```fortran
subroutine usr_init()
  use mod_global_parameters, only: par_files
  use mod_magnetic_topology, only: mt_params_read
  use mod_usr_methods, only: usr_special_convert

  ! Register the user conversion hook and read &magnetic_topology_list.
  usr_special_convert => topology_special_convert
  call mt_params_read(par_files)

  ! Continue with the usual problem-specific initialization here.
end subroutine usr_init

subroutine topology_special_convert(qunitconvert)
  use mod_magnetic_topology, only: mt_run_topology_task
  integer, intent(in) :: qunitconvert

  ! Dispatches the one task configured in &magnetic_topology_list.
  call mt_run_topology_task()
end subroutine topology_special_convert
```
