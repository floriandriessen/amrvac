# Radiation Synthesis Validation Notes

This note records the validation scope for the synthetic-emission extensions in
`mod_thermal_emission`.

## Compatibility baseline

The default `&emissionlist` behavior remains unchanged:

- `radiation_transfer='thin'`
- `ray_method='legacy'`
- `emission_model='auto'`
- `output_tau=.false.`
- `output_absorption_fraction=.false.`
- `instrument_postprocess=.false.`
- `radio_frequency=17.d9`
- `radio_beam_fwhm=0.d0`
- `radio_beam_pixel_size=0.d0`

Existing `convert_type='EIvtiCCmpi'`, `EIvtuCCmpi`, `SI...`, and `WI...`
inputs therefore continue to use the historical optically thin path unless a
new option is explicitly selected.

## Current product matrix

The current implementation is intended to keep the historical workflow
(`convert_type` plus `&emissionlist`) while adding explicit transfer and ray
choices for Cartesian post-processing.

| Product | Legacy instrument grid | Cartesian dat-resolution | Cartesian `cart_dda` arbitrary LOS | Thick transfer | Post-process to observational grid | Spherical |
| --- | --- | --- | --- | --- | --- | --- |
| EUV AIA/IRIS/EIS image | yes, thin | yes, axis-aligned thin/thick | yes, thin/thick | yes, H/He absorption | yes, AIA-style PSF for EUV images | legacy thin plus native `sph_intersection` thin/thick EUV instrument-grid |
| SXR image | yes, thin | yes, axis-aligned thin | no | no | no DDA post-process | legacy thin instrument-grid only |
| White light | spherical legacy | no | no | no | legacy LASCO-like PSF | yes, primary supported path |
| `radio_ff` | limited through EUV image path | yes | yes, thin/thick | yes | yes, independent Gaussian radio beam | no native support |
| `pseudo_current` | limited through EUV image path | yes | yes, thin only | no | no | no native support |

Practical recommendations:

- For cold-material EUV absorption in Cartesian data, use
  `radiation_transfer='thick'`; for non-axis-aligned views, also use
  `ray_method='cart_dda'`.
- For Cartesian stretched dat-resolution data, use VTU output. VTI is only safe
  for uniform image spacing.
- For radio quick-look products, the current `radio_ff` path is sufficient for
  multi-view brightness-temperature images. Use `radio_beam_fwhm` only when an
  observational beam is needed.
- For spherical data, use the legacy instrument-resolution projection for
  optically thin SXR and white light. For native EUV ray tracing, use
  `ray_method='sph_intersection'` with `dat_resolution=.false.` and
  `radiation_transfer='thin'` or `'thick'`.

## Local demo guards

The radiation-synthesis examples documented here are intended as local demo and
developer guards. They are not wired into the top-level autotest/CI suites by
default, because the full 3D convert cases are heavier than the regular
regression matrix.

## 3D demo case

The Cartesian demo case is `tests/demo4/RadiationSynthesis_3D`. It builds a
small TDm/prominence AMR cube and keeps two thick EUV reference views:

- axis-aligned optically thick AIA 171 image with `tau` and
  `absorption_fraction`, using `convert_171_thick.par`.
- oblique Cartesian DDA optically thick AIA 171 image with `tau` and
  `absorption_fraction`, using `convert_171_cart_dda_thick_oblique.par`.

Run the 3D Cartesian demo manually with:

```sh
make -C tests/demo4/RadiationSynthesis_3D -f test.make
```

The generated VTI files are intended for visual inspection in ParaView rather
than as bit-for-bit regression baselines.

## Performance baseline

The current performance-sensitive paths are:

- axis-aligned thick transfer, which streams LOS layers in small batches instead
  of allocating a full `Npix * Nlos` column buffer.
- `cart_dda`, which preselects block footprints on the image plane and routes
  per-pixel ray segments to owner MPI ranks before sorting. It emits lightweight
  profile counters for ray tests, ray hits, segment counts, MPI exchange volume,
  and sorting work in the per-conversion logs.

For quick performance comparisons, use the same TDm/prominence demo case and
record wall time for the retained reference conversions:

- `convert_171_thick.par`
- `convert_171_cart_dda_thick_oblique.par`

The most useful next profile points are ray-box intersection counts, collected
segment counts, `MPI_ALLTOALLV` volume, and same-pixel segment sorting time.

## Spherical-coordinate support

The default spherical path is the historical instrument-resolution projector.
It initializes a LOS/image-plane basis from `LOS_theta`, `LOS_phi`, and
`image_rotate`, sub-samples spherical cells in `r/theta/phi`, maps each
sub-sample to the image plane, and distributes the contribution through the
instrument PSF. A native ray/mesh-intersection path is also available for EUV
images with `ray_method='sph_intersection'`.

Supported today:

- optically thin EUV instrument-resolution images on spherical grids.
- native spherical thin and thick EUV images with
  `ray_method='sph_intersection'`, `dat_resolution=.false.`, and
  `radiation_transfer='thin'` or `'thick'`.
- spherical thick EUV diagnostic maps `tau` and `absorption_fraction`.
- optically thin SXR instrument-resolution images on spherical grids.
- LASCO-like white-light `B` and `pB` images on spherical grids.
- EUV spectra on spherical grids through the legacy instrument-resolution
  spectral path.
- occultation/opaque-inner-region handling through `R_opt_thick` for EUV/SXR
  visibility checks and the white-light occultor parameters for LASCO products.

Not supported yet:

- spherical .dat-resolution EUV/SXR images.
- `cart_dda` on spherical data.
- `instrument_postprocess` on spherical data.
- `radio_ff` or `pseudo_current` as validated spherical products.
- `sph_intersection` for SXR, white light, radio, pseudo-current, polar-axis
  crossing domains, or phi-wrapping domains.

The native spherical path constructs ray intersections with spherical grid
faces, sorts path points by LOS distance, and in thick EUV mode reuses the
`j/kappa` transfer layer for strict front-to-back integration. It still rejects
polar-axis crossing, phi-wrapping, non-EUV products, and dat-resolution output.

The spherical demo case is `tests/demo4/RadiationSynthesisSphericalTDm_3D`. It
builds a lightweight spherical TDm/prominence-style AMR shell and keeps two
native `sph_intersection` thick EUV reference views:

- top-down thick AIA 171, using `convert_171_sph_native_topdown_thick.par`.
- 45-degree top-corner thick AIA 171, using
  `convert_171_sph_native_corner45_top_thick.par`.

Run the spherical demo manually with:

```sh
make -C tests/demo4/RadiationSynthesisSphericalTDm_3D -f test.make
```

The generated VTI files contain `AIA171_thick`, `tau`, and
`absorption_fraction` and are intended for ParaView inspection.

## Remaining known scope limits

- `cart_dda` currently requires Cartesian slab dat-resolution grids.
- Stretched Cartesian dat-resolution images should use VTU output; VTI remains
  restricted to uniform image spacing.
- `sph_intersection` currently supports only spherical EUV
  instrument-resolution images.
- `pseudo_current` is optically thin only.
- `radio_ff` is a first thermal free-free implementation with a compact
  Gaunt-factor approximation, brightness-temperature output, and optional
  Gaussian beam post-processing.
