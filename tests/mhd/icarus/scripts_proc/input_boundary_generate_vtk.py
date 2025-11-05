#!/usr/bin/env python3
"""
Generate a VTK (STRUCTURED_POINTS, ASCII) inner boundary file for AMRVAC from EUHFORIA.dat file

- Outputs a single VTK with DIMENSIONS (nclt, nlon, ntime), SPACING (dtheta, dphi, dt),
  ORIGIN (theta0, phi0, 0). Scalars are written as: density, vr, pressure, Br (in code units).

Usage:
  python generate_inner_boundary_condition.py \
    --inputs "../solar_wind_bc_used_in_paper.in" \
    --output "../solar_wind_2015_single_gong.vtk" \
    --images  --images-dir "../previews" \
"""

import argparse
from pathlib import Path
from datetime import datetime
import glob
import numpy as np
import scipy.constants as spc
import matplotlib.pyplot as plt

carrington_synodic_rotation_period     = 27.2753*24.0  # rotation at a latitude of ~26° north or south,
carrington_sidereal_rotation_period    = 25.38*24.0 
solar_equator_synodic_rotation_period  = 26.24*24.0 
solar_equator_sidereal_rotation_period = 24.47*24.0 


# -------------------------- Defaults (can be overridden) --------------------------
DEFAULT_UNIT_LENGTH   = 6.957e8         # m (Rsun)
DEFAULT_UNIT_TIME     = 3600.0          # s (1 hour)
DEFAULT_UNIT_DENSITY  = 1.6726e-19      # kg/m^3 (so that rho~1 nondim around 1 AU in some setups)
# ----------------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Make VTK inner boundary for AMRVAC/EUHFORIA.")
    p.add_argument("--inputs", required=True, 
                   help="Glob pattern for input files (one or many).")
    p.add_argument("--output", default="./solar_wind_inner_boundary.vtk",  # The default should contain the date??
                   help="Output VTK path.")
    p.add_argument("--images", action="store_true", 
                   help="Write preview images per time slice.")
    p.add_argument("--images-dir", default="./",
                   help="Directory for preview images (created if missing).")
    p.add_argument("--rotation-period", type=float, default=carrington_synodic_rotation_period,  
                   help="Rigid rotation period in hours for corotation shift.")
    p.add_argument("--unit-length", type=float, default=DEFAULT_UNIT_LENGTH,
                   help="Code length unit in meters.")
    p.add_argument("--unit-time", type=float, default=DEFAULT_UNIT_TIME,
                   help="Code time unit in seconds.")
    p.add_argument("--unit-density", type=float, default=DEFAULT_UNIT_DENSITY,
                   help="Code density unit in kg/m^3.")
    p.add_argument("--verbose", action="store_true", help="More printout.")
    return p.parse_args()

def wrap_lon(lon, origin):
    return lon  - 2 * np.pi * np.floor( (lon - origin)/(2*np.pi))

def read_block_indices(lines):
    """Find start indices for field sections by their labels."""
    idx = {}
    for j, line in enumerate(lines):
        if line.strip() == "vr":
            idx["vr"] = j
        elif line.strip() == "number_density":
            idx["n"] = j
        elif line.strip() == "temperature":
            idx["T"] = j
        elif line.strip() == "Br":
            idx["Br"] = j
    return idx


def parse_one_file(path: Path, verbose=False):
    """Parse a single input file; returns dict with arrays and metadata."""
    with path.open() as f:
        lines = f.readlines()

    # Metadata lines (by your original structure)
    t_str = lines[1].strip()
    time_stamp = datetime.strptime(t_str, "%Y-%m-%dT%H:%M:%S")

    r_inner = float(lines[3].strip()) 

    n_clt = int(lines[5].strip())
    lon_num_index = 4 + 3 + n_clt + 1 #header: 4 lines, 3 lines clt-specs, clt points,  3 lines lon-specs, lon points, etc 
    n_lon = int(lines[lon_num_index].strip())

    # The lists:
    clt_list  = np.array(lines[7:lon_num_index-1], dtype=float)          # theta (colat)
    lons_list = np.array(lines[lon_num_index+2: lon_num_index+2+n_lon], dtype=float)

    # Field block indices:
    ids = read_block_indices(lines)
    if not {"vr", "n", "T", "Br"} <= ids.keys():
        raise ValueError(f"{path.name}: could not locate all field labels (vr, number_density, temperature, Br).")

    # Figure slice ranges between labels (robust to order)
    labels = sorted(ids.items(), key=lambda kv: kv[1])  # list of (name, index) sorted by appearance
    ranges = {}
    for k, (name, idx0) in enumerate(labels):
        idx1 = labels[k+1][1] if k+1 < len(labels) else len(lines)
        ranges[name] = (idx0+1, idx1)

    # Read raw arrays:
    vr_arr = np.array(lines[ranges["vr"][0]:ranges["vr"][1]], dtype=float)
    n_arr  = np.array(lines[ranges["n"][0] :ranges["n"][1] ], dtype=float)
    T_arr  = np.array(lines[ranges["T"][0] :ranges["T"][1] ], dtype=float)
    Br_arr = np.array(lines[ranges["Br"][0]:ranges["Br"][1]], dtype=float)

    if verbose:
        print(f"[{path.name}] n_clt={n_clt}, n_lon={n_lon}, len(vr)={vr_arr.size}")

    # Sanity checks
    expected = n_clt * n_lon
    for name, arr in [("vr", vr_arr), ("n", n_arr), ("T", T_arr), ("Br", Br_arr)]:
        if arr.size != expected:
            raise ValueError(f"{path.name}: {name} has {arr.size} entries, expected {expected} (n_clt*n_lon).")

    return {
        "time": time_stamp,
        "theta": clt_list,    # length n_clt
        "phi":   lons_list,   # length n_lon
        "vr":    vr_arr,
        "n":     n_arr,
        "T":     T_arr,
        "Br":    Br_arr,
        "n_clt": n_clt,
        "n_lon": n_lon,
        "r_inner": r_inner,
    }

def lon_shift_linear_uniform(A, delta_lon, dphi, axis=-1):
    """
    Shift A by delta_lon (radians) on a uniform-phi grid (spacing dphi),
    using integer roll + first-order linear interpolation for the residual.
    Positive delta_lon moves data to larger phi.
    """
    if dphi == 0.0 or delta_lon == 0.0:
        return A

    # keep delta in (-pi, pi] so residual stays small
    delta_lon = wrap_lon(delta_lon, - np.pi)

    delta_cols = delta_lon / dphi
    s = int(np.floor(delta_cols))       # integer roll
    a = float(delta_cols - s)           # residual in [0, 1)

    base = np.roll(A, s, axis=axis)
    if a == 0.0:
        return base

    # always blend toward the +1 neighbor after the base roll
    neighbor = np.roll(base, +1, axis=axis)
    return (1.0 - a) * base + a * neighbor


def make_grid(theta_list, phi_list):
    """Return 2D mesh for plotting; """
    clt2d, lon2d = np.meshgrid(theta_list, phi_list, indexing='ij')
    return clt2d, lon2d


def extend_periodic_phi(A: np.ndarray) -> np.ndarray:
    # A shape: (nclt, nlon). Append first φ column at the end.
    return np.concatenate([A, A[:, :1]], axis=1)


def main():
    args = parse_args()

    # Gather files
    files = sorted([Path(p) for p in glob.glob(args.inputs)])
    if not files:
        raise SystemExit(f"No input files match: {args.inputs}")

    # Units
    unit_length   = args.unit_length
    unit_time     = args.unit_time
    unit_velocity = unit_length / unit_time
    unit_density  = args.unit_density
    unit_b        = np.sqrt(spc.mu_0 * unit_velocity**2 * unit_density)
    unit_pressure = unit_velocity**2 * unit_density

    print("Information about the dataset and input boundary file for Icarus\n")
    print("Total number of magnetograms: ", len(files))

    # Parse all inputs
    recs = [parse_one_file(f, args.verbose) for f in files]
    t0 = recs[0]["time"]
    tN = recs[-1]["time"]
    print("Start time: ", t0)
    print("End time  : ", tN)
    total_hours = (tN - t0).total_seconds()/3600.0 if len(recs) > 1 else 0.0
    print("Total span [hours]: ", f"{total_hours:.3f}")
    if len(recs) > 1:
        avg_cadence = total_hours / (len(recs)-1)
    else:
        avg_cadence = 0.0
    print("Average timestep between magnetograms [hour]: ", f"{avg_cadence:.6f}\n")

    print("Units for generation (in SI)")
    print(f"Velocity  : {unit_velocity:.6e}")
    print(f"Density   : {unit_density:.6e}")
    print(f"Mag field : {unit_b:.6e}")
    print(f"Pressure  : {unit_pressure:.6e}")

    # Reference phi grid from first record
    nclt = recs[0]["n_clt"]
    nlon = recs[0]["n_lon"]
    theta_vals = recs[0]["theta"]
    phi_vals   = recs[0]["phi"]
    
    if args.verbose:
        print(f"Grid: nclt={nclt}, nlon={nlon}")

    # Derived longitude info
    phi0 = float(phi_vals.min())
    phi1 = float(phi_vals.max())
    Lphi = phi1 - phi0
    if Lphi <= 0:
        raise SystemExit("Longitude list must be strictly increasing.")
    #dphi  = Lphi / (nlon - 1)
    dphi  = 2*np.pi / nlon # only if grid  covers 2π uniformly

    assert np.isclose(Lphi / (nlon - 1), dphi), "This should be equal"
    #phi_extended =  np.concatenate([phi_vals, [phi_vals[-1] + dphi]])

    # theta spacing and origin (not required to be uniform; we just write an average spacing)
    # enforce uniformized theta for header.
    theta0 = float(theta_vals.min())
    theta1 = float(theta_vals.max())
    dtheta = (theta1 - theta0) / (nclt - 1)

    vr_all   = []
    dens_all = []
    pres_all = []
    Br_all   = []

    # For previews
    if args.images:
        Path(args.images_dir).mkdir(parents=True, exist_ok=True)

    # Time loop
    shift_accum = 0.0
    prev_time = None

    # Shift lons such that they are all in [0, 2*pi]
    phi_vals = wrap_lon((phi_vals - phi0),0)  # now in [0, 2π) if your list spans full range
    shift_accum  = phi0
    phi0 = float(phi_vals.min())
    phi1 = float(phi_vals.max())

    for i, R in enumerate(recs):
        # cadence
        if prev_time is None:
            dt_hours = 0.0
            prev_time = R["time"]
        else:
            dt_hours = (R["time"] - prev_time).total_seconds()/3600.0
            prev_time = R["time"]

        # Optional rigid rotation shift
        if args.rotation_period > 0.0:
            delta_lon    = - dt_hours * (2.0*np.pi / args.rotation_period) 
            shift_accum += delta_lon
            shift_accum  = wrap_lon(shift_accum, -np.pi) 
        else:
            delta_lon = 0.0

        # Scalars in SI (from file), convert to code units
        mp = spc.m_p
        # number density -> mass density
        rho_SI = R["n"] * mp * 0.5                       # kg/m^3
        p_SI   = spc.k * R["n"] * R["T"]                 # Pa
        vr_SI  = R["vr"]                                 # m/s
        Br_SI  = R["Br"]                                 # Tesla (assuming input is SI)

        # Non-dimensionalize
        vr  = (vr_SI  / unit_velocity).reshape((nlon, nclt)).T
        rho = (rho_SI / unit_density).reshape((nlon, nclt)).T
        pre = (p_SI   / unit_pressure).reshape((nlon, nclt)).T
        Br  = (Br_SI  / unit_b).reshape((nlon, nclt)).T

        if shift_accum != 0.0:
            vr  = lon_shift_linear_uniform(vr,  shift_accum, dphi, axis=1)
            rho = lon_shift_linear_uniform(rho, shift_accum, dphi, axis=1)
            pre = lon_shift_linear_uniform(pre, shift_accum, dphi, axis=1)
            Br  = lon_shift_linear_uniform(Br,  shift_accum, dphi, axis=1)

        # Enforce periodic seam by adding an extra column (last col == first col)
        # NOT NEEDED?
        #vr, rho, pre, Br = [extend_periodic_phi(A) for A in (vr, rho, pre, Br)] 
        #nlon+=1

        vr_all.append(vr.T.flatten())
        dens_all.append(rho.T.flatten())
        pres_all.append(pre.T.flatten())
        Br_all.append(Br.T.flatten())
        # Optional previews
        if args.images:
            cmap = 'turbo'

            # For plotting only:  1D axes in degrees 
            lon_deg = np.rad2deg(phi_vals)                      # (nlon,)
            lat_deg = np.rad2deg((np.pi/2.0) - theta_vals)      # (nclt,)

            # For plotting only: make latitude ascend (south→north) and reorder fields accordingly
            order = np.argsort(lat_deg)                         # no-op if already ascending
            lat_plot = lat_deg[order]
            vr_plot  = vr[order, :]
            rho_plot = rho[order, :]
            pre_plot = pre[order, :]
            Br_plot  = Br[order, :]

            # Temperature for plotting (use reordered arrays)
            mp = spc.m_p
            n_SI_plot = rho_plot*unit_density/(0.5*mp)
            T_SI_plot = (pre_plot*unit_pressure)/(spc.k*n_SI_plot)

            fig, axs = plt.subplots(4, 1, figsize=(8, 10), sharex=True)

            im1 = axs[0].pcolormesh(lon_deg, lat_plot, vr_plot*unit_velocity * 1e-3, shading="auto", cmap=cmap)
            c1 = fig.colorbar(im1, ax=axs[0]); c1.set_label(r"$V_r$ [km s$^{-1}$]")
            axs[0].set_ylabel("lat [deg]")

            im2 = axs[1].pcolormesh(lon_deg, lat_plot, n_SI_plot * 1e-6, shading="auto", cmap=cmap)
            c2 = fig.colorbar(im2, ax=axs[1]); c2.set_label(r"$n$ [cm$^{-3}$]")
            axs[1].set_ylabel("lat [deg]")

            im3 = axs[2].pcolormesh(lon_deg, lat_plot, T_SI_plot * 1e-6, shading="auto", cmap=cmap)
            c3 = fig.colorbar(im3, ax=axs[2]); c3.set_label(r"$T$ [MK]")
            axs[2].set_ylabel("lat [deg]")

            im4 = axs[3].pcolormesh(lon_deg, lat_plot, Br_plot*unit_b * 1e9, shading="auto", cmap='seismic')
            c4 = fig.colorbar(im4, ax=axs[3]); c4.set_label(r"$B_r$ [nT]")
            axs[3].set_ylabel("lat [deg]")
            axs[3].set_xlabel("lon [deg]")

            fstem = f"BC_{(args.output).split('/')[-1].split('.vtk')[0]}_{i:05d}.png"
            fig.tight_layout()
            fig.savefig(Path(args.images_dir)/fstem, dpi=120)
            plt.close(fig)

    # Concatenate time
    vr_all   = np.concatenate(vr_all)
    dens_all = np.concatenate(dens_all)
    pres_all = np.concatenate(pres_all)
    Br_all   = np.concatenate(Br_all)

    ntime = len(recs)
    dt_out = avg_cadence  # hours between frames
    print("\nVTK metadata:")
    print(f" DIMENSIONS: ({nclt}, {nlon}, {ntime})")
    print(f" SPACING   : (dθ={dtheta:.6e}, dφ={dphi:.6e}, dt={dt_out:.6f} h)")
    print(f" ORIGIN    : (θ0={theta0:.6e}, φ0={phi0:.6e}, t0=0)")

    # Write VTK
    out = Path(args.output)
    with out.open("w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("File generated by make_icarus_vtk.py\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nclt} {nlon} {ntime}\n")
        f.write(f"SPACING {dtheta:.16e} {dphi:.16e} {dt_out:.16e}\n")
        f.write(f"ORIGIN {theta0:.16e} {phi0:.16e} 0\n")
        f.write(f"POINT_DATA {nclt*(nlon)*ntime}\n")

        def write_scalar_block(name, arr):
            f.write(f"SCALARS {name} float 1\n")
            f.write("LOOKUP_TABLE default\n")
            # VTU readers tolerate line breaks anywhere; write in chunks for readability
            for val in arr:
                f.write(f"{float(val):.8e} ")
            f.write("\n")

        write_scalar_block("density", dens_all)
        write_scalar_block("vr",      vr_all)
        write_scalar_block("pressure", pres_all)
        write_scalar_block("Br",      Br_all)

    # Report delta_phi suggestion (starting longitude)
    delta_phi_hint = float(phi_vals[0])
    print("\nSuggested delta_phi for the .par file (radians):", delta_phi_hint)
    print("Generated VTK:", out.resolve())
    if args.images:
        print("Preview images in:", Path(args.images_dir).resolve())


if __name__ == "__main__":
    main()

