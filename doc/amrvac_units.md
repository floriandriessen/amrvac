# Units

Describes the
- `hd_physical_units` subroutine in `hd/mod_hd_phys.t`.
- `mhd_physical_units` subroutine in `mhd/mod_mhd_phys.t`.

Consider two physical quantities $A$ and $B$ related by $c$, so that $A=cB$. Write $A$ in units of $A_0$ as $A=A_0\widetilde A$, with $A$ variable and dimensional, $A_0$ constant and dimensional, and $\widetilde A$ variable and dimensionless. If we choose $\widetilde A=\widetilde B$, we have $A_0=cB_0$. The user specifies initial/boundary conditions in $\widetilde X$-quantities or scaled $\widetilde X$-quantities, and MPI-AMRVAC solves $\widetilde X$-quantities.  
More generally of course $\widetilde X$-quantities are related, e.g. if $A=cBD^2$, then $\widetilde A=\widetilde B\widetilde D^2$ and $A_0=cB_0D_0^2$.

MPI-AMRVAC constructs a *dependent* set of units $\{X_0\}$ for the set of physical quantities $\{X\}$ based on a set of *independent but complete* units $\{Y_0\}$ the user provides, various combinations are possible. If a unit $X_0$ cannot be determined since the user does not provide a fully complete set of independent units $\{Y_0\}$, $X_0=1$ is used.  
Constructing a set of units isn't always necessary since the user may simply provide $\widetilde X$-data, but the user may actually want their data to represent $a\widetilde A$ with some scale factor $a$, or MPI-AMRVAC may need to use data from tables which must be converted appropriately to the correct set of units, e.g. for optically thin radiative cooling from `physics/mod_radiative_cooling.t` used in e.g. `hd/mod_hd_phys.t`, specified in tables of $\log_{10}\text{ K}$. In such cases, the user should specify a complete set of units $\{Y_0\}$ in the `usr_init` subroutine of the `mod_usr.t` file.  
The constants $c$ are in general dimensional. Their dimensions are set either in SI units or Gaussian CGS units by MPI-AMRVAC (the user chooses), and these constants may depend also on additional dimensionless quantities the user provides, e.g. ionization levels.

MPI-AMRVAC's internal set of units for hydrodynamics/magnetohydrodynamics are:
- mass ($m_0$)
- length ($L_0$)
- time ($t_0$)
- velocity ($v_0$)
- number density ($n_0$)
- density ($\rho_0$)
- temperature ($T_0$)
- pressure ($P_0$)
- magnetic field ($B_0$)
- electric charge ($q_0$)

They are related by:
- $v_0=L_0/t_0$
- $m_0=\rho_0L_0^3$
- $P_0=\rho_0v_0^2$
- $\rho_0=am_\text p n_0$
- $P_0=bk_\text Bn_0T_0$
- $P_0=B_0^2/\mu_0$
- $q_0=(B_0L_0^2)/(v_0\mu_0)$ (SI)
- $q_0=(B_0cL_0^2)/(v_0^2\mu_0)$ (Gaussian CGS, $c$ is speed of light)

Some information on the constants used.  
The proton mass `mp` is needed to relate number density and mass density, and the Boltzmann constant `kB` is needed to relate temperature, density, and pressure. They may be specified in either SI or Gaussian CGS units by MPI-AMRVAC.  
Helium abundance and ionization levels also affect these relationships.  
Electrons are assumed to be massless, neutrons are assumed to have proton mass, and binding energies (atomic and nuclear) are ignored.  
The constant $\mu_0$ in Gaussian CGS units is $1$, but electromagnetic equations in Gaussian CGS units are off by a factor of $4\pi$. MPI-AMRVAC defines $\mu_0$ to be $4\pi$ to ensure consistent conversions.  
The electric charge units are necessary for particle motion only.

Option 1: `eq_state_units`  
With `eq_state_units`, two molecular scaling factors are used, `a` for mass and `b` for the 'degrees of freedom per mass'.  
MPI-AMRVAC understands number density units in this case as **H number density units** (ionized or not). So `a` is the ratio of fluid mass to H component mass, `a = 1+4*He_abundance`. And this is used as `unit_density = a*mp*unit_numberdensity`.  
So we have $\widetilde\rho=\widetilde n$ and $n_0=am_pn_0$.  
`b` similarly associates number density units to pressure units. It is the number of degrees of freedom that come with one H proton, the unit pressure per H proton per unit temperature, using the equipartition theorem. Each H proton contributes `H_ion_fr` electrons, as well as `He_abundance` He atoms, which each contribute one alpha particle as well as `He_ion_fr` 1st electrons and `He_ion_fr*He_ion_fr2` 2nd electrons. So `b = 1+H_ion_fr+He_abundance*(He_ion_fr*(He_ion_fr2+1)+1)` in general (if `hd_partial_ionization`). If fully ionized (not `hd_partial_ionization`) this simplifies to `b = 2+3*He_abundance`. Now `unit_pressure = b*kB*unit_numberdensity*unit_temperature`, so $P_0=(b/a)(k_B/m_p)\rho_0T_0$ and $\widetilde P=\widetilde n\widetilde T=\widetilde\rho\widetilde T$.

Option 2: not `eq_state_units`  
In this mode `a` and `b` are kept fixed to `1`, so that number density units count **baryon number** (any form of mass). Instead of the units being affected by ionization, ionization affects the $\widetilde X$-equations directly. `b/a` is replaced by `RR`, and $\widetilde P=\text{RR}\widetilde\rho\widetilde T$. `kB` and `mp` are still used to relate $\rho_0$ to $n_0$, and $P_0$ to $n_0$ and $T_0$.
