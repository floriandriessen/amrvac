# Magnetic Topology and QSL Products

This module provides field-line topology products for Cartesian MHD data:

- field-line length;
- twist number `Tw`;
- ordinary-like axis-plane `Q` computed as Scott q0 by tangent-vector transport;
- perpendicular squashing factor `Qperp`.

The user-facing `q` / `logQ` product is Scott q0. It is not the former endpoint finite-difference Q, and it is not the Qperp scalar renamed as Q. Qperp remains an independent perpendicular squashing diagnostic.

## Product Policy

| Workflow | Public products | Notes |
| --- | --- | --- |
| `axis_plane_full_vtu` | length, `Tw`, `logQ`, `connection_type_Q`, Qperp | `length_total`, `logQ`, and Qperp come from the Scott/Qperp tangent trace. `connection_type_Q` is the minimal mask field. |
| `axis_plane_csv` | length, optional `Tw`, optional Qperp | CSV mode does not write ordinary Q. |
| `seed_products` | length, optional `Tw`, optional Qperp | No ordinary Q or mapping product. |
| `arbitrary_plane_products` | length, optional `Tw`, optional Qperp | No ordinary Q or mapping product. |
| `volume_vti` | length, optional `Tw`, optional Qperp | No volume ordinary Q or mapping product. |

Former endpoint finite-difference Q, robust stencil modes, endpoint-stencil quality diagnostics, mixed-source Q fields, and Q-specific `*_vis` aliases are not user-facing products.

## Axis-Plane VTU Fields

Axis-plane QSL 1.0 minimal VTU output contains:

- `length_total`
- `twist_total`
- `logQ`
- `logQperp`
- `connection_type_Q`

Full-detail axis-plane VTU may include raw and diagnostic endpoint/tangent fields, validity/status metadata, mapping metadata, and twist component diagnostics, but it does not expose old endpoint-FD, mixed-source, Scott-alias, or Q-visual-alias fields. If raw `q` is written, it is the same Scott q0 product represented by `logQ`.

Use `connection_type_Q` for minimal-output visualization masks in ParaView. No separate `Q visual alias` field is written.

## OpenMP Post-Processing

The tested OpenMP route is the gfortran `ARCH=openmp` build, which enables
`-fopenmp`. The default `ARCH=default` build is serial/non-OpenMP. No Intel
`ifx`/`ifort` OpenMP architecture is added here.

Build the AMRVAC4 magnetic-topology QSL demo with OpenMP, for example:

```sh
cd /home/nanami/codes/amrvac_qsl/tests/demo/AMRVAC4_MagneticTopologyQSL_3D
AMRVAC_DIR=/home/nanami/codes/amrvac_qsl make ARCH=openmp -j4
```

Run one process with OpenMP threads for local post-processing:

```sh
OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores \
  ./amrvac -i par_showcase/showcase_ref_amr_l3_yz_vertical_512.par
```

The AMRVAC4 demo organizes reproducible TDm/RBSL snapshot setup and QSL
post-processing examples under
`tests/demo/AMRVAC4_MagneticTopologyQSL_3D/`: snapshot-generation par files are
in the demo root, showcase examples in `par_showcase/`, seed-resolution examples
in `par_resolution/`, and timing/efficiency examples in `par_benchmark/`.

OpenMP currently accelerates seed-level tracing for `axis_plane_full_vtu` with
both `mt_vtk_detail = 'minimal'` and `mt_vtk_detail = 'full'`. It also covers
the seed-level product loops for `arbitrary_plane_products`, `seed_products`,
and `volume_vti`. VTU/VTI/CSV writing remains serial. Start with
`OMP_NUM_THREADS=4` or `8` and benchmark locally; avoid oversubscription such
as running many MPI ranks while each rank also uses many OpenMP threads. Very
large volume products can still be expensive because snapshot loading, final
product arrays, and VTI writing are not distributed.

## Namelist Runner

Typical axis-plane VTU use:

```fortran
&magnetic_topology_list
  mt_enable = .true.
  mt_mode = 'axis_plane_full_vtu'
  mt_output_file = 'qsl_axis_xy.vtu'
  mt_plane = 'xy'
  mt_xmin = -1.d0
  mt_xmax =  1.d0
  mt_nx = 128
  mt_ymin = -1.d0
  mt_ymax =  1.d0
  mt_ny = 128
  mt_z0 = 0.d0
  mt_dL = 0.05d0
  mt_max_steps = 1000
  mt_vtk_detail = 'minimal'
  mt_compute_qperp = .true.
/
```

`axis_plane_full_vtu` always computes the Scott/Qperp trace used by `length_total`, `logQ`, and Qperp. It also writes `twist_total` as the axis-plane twist integral product. `mt_compute_twist` and `mt_compute_qperp` are ignored in this mode because the coordinate-plane VTU product set is fixed by `mt_vtk_detail`.

Typical axis-plane CSV use:

```fortran
&magnetic_topology_list
  mt_enable = .true.
  mt_mode = 'axis_plane_csv'
  mt_output_prefix = 'qsl_axis_csv'
  mt_plane = 'xy'
  mt_xmin = -1.d0
  mt_xmax =  1.d0
  mt_nx = 128
  mt_ymin = -1.d0
  mt_ymax =  1.d0
  mt_ny = 128
  mt_z0 = 0.d0
  mt_dL = 0.05d0
  mt_max_steps = 1000
  mt_compute_twist = .true.
  mt_compute_qperp = .true.
/
```

The CSV runner writes length by default, plus twist and Qperp when enabled. It does not write ordinary Q.

## Direct Wrappers

Coordinate-plane wrappers:

- `mt_qsl_plane_vtu_xy/xz/yz`: ParaView-ready axis-plane VTU with Scott q0 `Q` and optional Qperp.
- `mt_topology_plane_xy/xz/yz`: shared length/twist/mapping CSV trace, without ordinary Q.
- `mt_qperp_plane_xy/xz/yz`: Qperp-only CSV.

Other product wrappers:

- `mt_fieldline_products_seeds`: seed-set length, optional twist, optional Qperp.
- `mt_fieldline_products_plane_arbitrary`: arbitrary orthonormal-plane length, optional twist, optional Qperp.
- `mt_fieldline_products_volume_vti`: Cartesian volume length, optional twist, optional Qperp.

## Interpretation

`logQ` and `logQperp` are related but distinct diagnostics. `logQ` is an ordinary-like Scott q0 quantity using endpoint boundary normals and tangent-vector transport. `logQperp` is the perpendicular squashing factor. For vertical and arbitrary cuts, Qperp remains the primary squashing diagnostic; axis-plane `logQ` is provided as an ordinary-like companion product. In minimal output, `connection_type_Q` is the public mask field; full-detail output may include additional validity/status metadata.

Unsupported or intentionally absent products:

- endpoint finite-difference endpoint-FD Q;
- endpoint-stencil diagnostics;
- mixed-source Q fields;
- Q-specific visual aliases;
- arbitrary-plane ordinary Q;
- seed-set ordinary Q;
- volume ordinary Q or `Qlocal`;
- AMR-native topology output.

MPI distributed tracing, RK45/adaptive tracing, and volume/AMR scaling
performance are outside the axis-plane QSL 1.0 scope unless validated and
documented separately.
