# Magnetic Topology/QSL Regression

This regression covers the current public topology/QSL products:

- Qperp on coordinate planes, arbitrary planes, seeds, and volume products;
- Scott q0 `logQ` in axis-plane VTU output;
- length, twist, mapping CSV consistency where those products remain public;
- public VTU/VTI schema checks.

The endpoint finite-difference endpoint-FD Q path has been removed from the user-facing module and is not part of this regression. The checker also verifies that old endpoint-FD, mixed-source, and Q-specific visual-alias fields are absent from product files.

Run:

```sh
AMRVAC_DIR=/home/nanami/codes/amrvac_qsl make ARCH=debug -j4
./amrvac -i amrvac.par
python3 check_topology_qsl.py qsl_regression_summary.csv
```
