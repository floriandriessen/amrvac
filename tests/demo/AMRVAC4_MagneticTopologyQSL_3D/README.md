# AMRVAC4 Magnetic Topology QSL 3D Demo

This demo prepares TDm/RBSL-like snapshots for magnetic topology/QSL
post-processing. The files in this directory initialize snapshots only; they do
not run field-line tracing and do not write QSL VTU products.

The flux-rope contribution uses a direct RBSL magnetic-field formulation based
on Titov et al. (2021) Appendix A/B. This demo-local initialization code
accumulates the RBSL magnetic field directly; it does not use the older
`AI`/`AF` vector-potential accumulation plus finite-difference curl route. This
does not change AMRVAC core QSL tracing or topology post-processing.

Build from this directory, for example:

```sh
AMRVAC_DIR=/home/nanami/codes/amrvac make ARCH=default -j4
```

Run one snapshot initialization:

```sh
mkdir -p output
/usr/bin/time -v ./amrvac -i init_u128_snapshot.par
```

Generated `.dat`, `.log`, and `output/` files are ignored and should not be
committed.

MPI can be used for snapshot initialization. OpenMP acceleration mainly applies
to later QSL seed tracing/post-processing, not necessarily to this initialization
path.

QSL post-processing examples are organized into three par groups:

- `par_showcase/`: five AMR l3 showcase products using the direct-RBSL AMR
  snapshot: point-set, xy bottom plane, yz vertical plane, oblique plane, and
  a high-resolution volume cloud.
- `par_resolution/`: a yz-plane resolution series on the same AMR l3 snapshot:
  128^2, 256^2, 512^2, and 1024^2 seeds.
- `par_benchmark/`: yz-plane and volume-cloud benchmark pars for the uniform
  u128/u256/u512 snapshots and AMR l3/l4 snapshots.

Use one process with OpenMP threads for QSL post-processing runs, for example:

```sh
OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores \
  ./amrvac -i par_showcase/showcase_ref_amr_l3_yz_vertical_512.par
```

Seed-level OpenMP covers the axis-plane VTU path, arbitrary/oblique plane
products, point-set products, and volume VTI products. Output writing is still
serial, and large volume clouds remain expensive because final product arrays
and VTI output scale with the full seed count. The 512^3 and 1024^3 volume
pars are opt-in/manual examples; do not run them as a default demo benchmark.

Benchmark timing is end-to-end: snapshot loading, setup, seed tracing, and
serial output are all included. For repeated thread-count runs, keep generated
per-thread pars and summaries under `output/benchmark/`; `output/` is ignored
and should not be committed.

## Snapshot par files

| file | grid | AMR | output |
| --- | --- | --- | --- |
| `init_u128_snapshot.par` | uniform 128^3 | no | `output/qsl_tdm_direct_u128_snap0000.dat` |
| `init_u256_snapshot.par` | uniform 256^3 | no | `output/qsl_tdm_direct_u256_snap0000.dat` |
| `init_u512_snapshot.par` | uniform 512^3 | no | `output/qsl_tdm_direct_u512_snap0000.dat` |
| `init_amr_base128_l3_snapshot.par` | base 128^3 | max level 3 | `output/qsl_tdm_direct_amr_base128_l3_snap0000.dat` |
| `init_amr_base128_l4_snapshot.par` | base 128^3 | max level 4 | `output/qsl_tdm_direct_amr_base128_l4_snap0000.dat` |

Recommended inputs are:

- default demo: uniform 128^3
- medium demo: uniform 256^3
- high-resolution demo: uniform 512^3
- AMR smoke/demo: base 128^3 with max level 3
- AMR high-resolution smoke: base 128^3 with max level 4

The AMR case uses a base 128^3 mesh with `refine_max_level=3`, so the maximum
local effective resolution is about uniform 512^3 under the current AMRVAC level
convention. It is an AMR demo/smoke case, not a uniform 512^3 replacement across
the whole domain.

The uniform 512^3 case is a high-resolution demo input. On the local benchmark
machine it produced an about 8.1 GB snapshot, used about 57 GiB MaxRSS, and took
about 23 minutes. Prefer reusing an existing u512 snapshot instead of repeatedly
regenerating it.

QSL conversion examples live in the magnetic topology documentation and TDm
benchmark notes. The stable minimal axis-plane QSL fields are `length_total`,
`twist_total`, `logQ`, `logQperp`, and `connection_type_Q`.
