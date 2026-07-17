#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  tests/demo4/RadiationSynthesisSphericalTDm_3D/run_u512_np8.sh

Environment:
  NP          MPI size to run. Default: 8
  NP_LIST     Optional MPI sizes to run sequentially, e.g. "16 32 64".
              If set, NP is ignored. Commas are also accepted.
  MPIRUN      MPI launcher. Default: mpirun
  MAKE        Make command. Default: make
  RUN_MAKE    yes | no. Default: no
  LOG_FILE    Benchmark log path. Default: timestamped file under output/benchmark/logs
  TIMEOUT     Timeout in seconds. Default: 3600. Set to 0 to disable.
  BATCH_SIZE  Optional ghost-cell communication batch size override, e.g. 64.
  PIXEL_BATCH Optional radsyn_pixel_batch override, e.g. 64, 256, or 512.
  SEG_BATCH_FACTOR Optional radsyn_segment_batch_factor override, e.g. 0 or 256.
                   0 uses automatic memory-budget sizing.
  SEG_MEMORY_MB Optional radsyn_segment_memory_mb override for auto sizing, e.g. 8 or 32.
  SEG_COMM_FACTOR  Optional radsyn_segment_comm_factor override, e.g. 32.
  RADSYN_VERBOSE yes | no. Default: no. Enables radiation-synthesis profiling output.
  AMRVAC_DIR  AMRVAC repository root. Default: inferred from this script.
  FORCE       yes | no. Default: no. Set yes to run even if another u512 job is found.

This script runs only the spherical u512 corner45 thick benchmark.
It expects output/radsyn_sph_prom_u512_0000.dat to already exist.
EOF
}

case "${1:-}" in
  -h|--help) usage; exit 0 ;;
  "") ;;
  *) usage; exit 2 ;;
esac

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="${AMRVAC_DIR:-$(cd "${script_dir}/../../.." && pwd)}"
export AMRVAC_DIR="${repo_root}"

np_list_raw="${NP_LIST:-${NP:-8}}"
mpirun_cmd="${MPIRUN:-mpirun}"
make_cmd="${MAKE:-make}"
run_make="${RUN_MAKE:-no}"
force="${FORCE:-no}"
timeout_seconds="${TIMEOUT:-3600}"
batch_size="${BATCH_SIZE:-}"
pixel_batch="${PIXEL_BATCH:-}"
seg_batch_factor="${SEG_BATCH_FACTOR:-}"
seg_memory_mb="${SEG_MEMORY_MB:-}"
seg_comm_factor="${SEG_COMM_FACTOR:-}"
radsyn_verbose="${RADSYN_VERBOSE:-no}"
np_list_raw="${np_list_raw//,/ }"
read -r -a np_values <<< "${np_list_raw}"

case "${run_make}" in
  yes|no) ;;
  *) echo "RUN_MAKE must be yes or no" >&2; exit 2 ;;
esac

case "${force}" in
  yes|no) ;;
  *) echo "FORCE must be yes or no" >&2; exit 2 ;;
esac

if [[ "${#np_values[@]}" -eq 0 ]]; then
  echo "NP_LIST must contain at least one MPI size" >&2
  exit 2
fi

for np_value in "${np_values[@]}"; do
  if ! [[ "${np_value}" =~ ^[1-9][0-9]*$ ]]; then
    echo "MPI size must be a positive integer: ${np_value}" >&2
    exit 2
  fi
done

if [[ "${#np_values[@]}" -gt 1 && -n "${LOG_FILE:-}" ]]; then
  echo "LOG_FILE can only be used with a single MPI size" >&2
  exit 2
fi

if ! [[ "${timeout_seconds}" =~ ^[0-9]+$ ]]; then
  echo "TIMEOUT must be a non-negative integer number of seconds" >&2
  exit 2
fi

if [[ -n "${batch_size}" ]] && ! [[ "${batch_size}" =~ ^[1-9][0-9]*$ ]]; then
  echo "BATCH_SIZE must be a positive integer" >&2
  exit 2
fi

if [[ -n "${pixel_batch}" ]] && ! [[ "${pixel_batch}" =~ ^[1-9][0-9]*$ ]]; then
  echo "PIXEL_BATCH must be a positive integer" >&2
  exit 2
fi

if [[ -n "${seg_batch_factor}" ]] && ! [[ "${seg_batch_factor}" =~ ^[0-9]+$ ]]; then
  echo "SEG_BATCH_FACTOR must be a non-negative integer" >&2
  exit 2
fi

if [[ -n "${seg_memory_mb}" ]] && ! [[ "${seg_memory_mb}" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; then
  echo "SEG_MEMORY_MB must be a positive number" >&2
  exit 2
fi

if [[ -n "${seg_memory_mb}" && "${seg_memory_mb}" =~ ^0*([.]0*)?$ ]]; then
  echo "SEG_MEMORY_MB must be a positive number" >&2
  exit 2
fi

if [[ -n "${seg_comm_factor}" ]] && ! [[ "${seg_comm_factor}" =~ ^[1-9][0-9]*$ ]]; then
  echo "SEG_COMM_FACTOR must be a positive integer" >&2
  exit 2
fi

case "${radsyn_verbose}" in
  yes|no) ;;
  *) echo "RADSYN_VERBOSE must be yes or no" >&2; exit 2 ;;
esac

cd "${script_dir}"

snapshot="output/radsyn_sph_prom_u512_0000.dat"
bench_par="par_benchmark/bench_u512_corner45_thick.par"
mkdir -p output/benchmark/logs

if [[ ! -f "${snapshot}" ]]; then
  echo "Missing ${snapshot}; create it first with par_init/u512.par." >&2
  exit 1
fi

if [[ "${run_make}" == "yes" ]]; then
  "${make_cmd}" -f test.make setup
fi

run_one_np() {
  local np="$1"
  local timestamp log_file override_par start_ts end_ts elapsed exit_code
  local run_cmd
  local child_pid

  if [[ "${force}" == "no" ]]; then
    if pgrep -af "prterun -np ${np} ./amrvac -i ${bench_par}|mpirun -np ${np} ./amrvac -i ${bench_par}" >/dev/null; then
      echo "An existing spherical u512 np=${np} benchmark appears to be running." >&2
      echo "Stop it first, or rerun with FORCE=yes if this is intentional." >&2
      return 1
    fi
  fi

  timestamp="$(date +%Y%m%d_%H%M%S)"
  log_file="${LOG_FILE:-output/benchmark/logs/spherical_u512_np${np}_${timestamp}.log}"
  override_par=""
  if [[ -n "${batch_size}" || -n "${pixel_batch}" || -n "${seg_batch_factor}" || \
        -n "${seg_memory_mb}" || \
        -n "${seg_comm_factor}" || "${radsyn_verbose}" == "yes" ]]; then
    override_par="output/benchmark/logs/spherical_u512_np${np}_${timestamp}_override.par"
    {
      if [[ -n "${batch_size}" ]]; then
        cat <<EOF
&methodlist
  ghostcell_comm_batched=.true.
  ghostcell_comm_batch_size=${batch_size}
/
EOF
      fi
      if [[ -n "${pixel_batch}" || -n "${seg_batch_factor}" || -n "${seg_memory_mb}" || \
            -n "${seg_comm_factor}" || "${radsyn_verbose}" == "yes" ]]; then
        printf '&emissionlist\n'
        if [[ -n "${pixel_batch}" ]]; then
          printf '  radsyn_pixel_batch=%s\n' "${pixel_batch}"
        fi
        if [[ -n "${seg_batch_factor}" ]]; then
          printf '  radsyn_segment_batch_factor=%s\n' "${seg_batch_factor}"
        fi
        if [[ -n "${seg_memory_mb}" ]]; then
          printf '  radsyn_segment_memory_mb=%s\n' "${seg_memory_mb}"
        fi
        if [[ -n "${seg_comm_factor}" ]]; then
          printf '  radsyn_segment_comm_factor=%s\n' "${seg_comm_factor}"
        fi
        if [[ "${radsyn_verbose}" == "yes" ]]; then
          printf '  radsyn_verbose=.true.\n'
        fi
        printf '/\n'
      fi
    } > "${override_par}"
  fi

  echo "-- spherical/u512: convert np=${np}"
  echo "-- log: ${log_file}"
  if [[ -n "${override_par}" ]]; then
    echo "-- override: ${override_par}"
  fi

  if [[ "${timeout_seconds}" == "0" ]]; then
    echo "-- timeout: disabled"
    if [[ -n "${override_par}" ]]; then
      run_cmd=("${mpirun_cmd}" -np "${np}" ./amrvac -i "${bench_par}" "${override_par}")
    else
      run_cmd=("${mpirun_cmd}" -np "${np}" ./amrvac -i "${bench_par}")
    fi
  else
    echo "-- timeout: ${timeout_seconds}s"
    if [[ -n "${override_par}" ]]; then
      run_cmd=(timeout "${timeout_seconds}" "${mpirun_cmd}" -np "${np}" ./amrvac -i "${bench_par}" "${override_par}")
    else
      run_cmd=(timeout "${timeout_seconds}" "${mpirun_cmd}" -np "${np}" ./amrvac -i "${bench_par}")
    fi
  fi

  start_ts="$(date +%s)"
  child_pid=""
  cleanup_child() {
    local signal_name="$1"
    trap - INT TERM
    if [[ -n "${child_pid}" ]]; then
      echo "-- spherical/u512: stopping np=${np} after ${signal_name}" >&2
      kill -TERM "-${child_pid}" 2>/dev/null || kill -TERM "${child_pid}" 2>/dev/null || true
      sleep 2
      kill -KILL "-${child_pid}" 2>/dev/null || kill -KILL "${child_pid}" 2>/dev/null || true
      wait "${child_pid}" 2>/dev/null || true
      child_pid=""
    fi
    case "${signal_name}" in
      INT) exit 130 ;;
      TERM) exit 143 ;;
      *) exit 1 ;;
    esac
  }
  trap 'cleanup_child INT' INT
  trap 'cleanup_child TERM' TERM
  set +e
  setsid "${run_cmd[@]}" > "${log_file}" 2>&1 &
  child_pid="$!"
  wait "${child_pid}"
  exit_code="$?"
  set -e
  if [[ "${exit_code}" -ne 0 ]]; then
    kill -TERM "-${child_pid}" 2>/dev/null || true
    sleep 1
    kill -KILL "-${child_pid}" 2>/dev/null || true
  fi
  child_pid=""
  trap - INT TERM

  if [[ "${exit_code}" -eq 0 ]]; then
    end_ts="$(date +%s)"
    elapsed="$((end_ts - start_ts))"
    echo "-- spherical/u512: np=${np} completed in ${elapsed}s"
    if [[ ! -f output/benchmark/radsyn_benchmark_times.csv ]]; then
      printf 'suite,case,np,seconds,parfile,logfile\n' > output/benchmark/radsyn_benchmark_times.csv
    fi
    printf '%s,%s,%s,%s,%s,%s\n' \
      "spherical" "u512" "${np}" "${elapsed}" "${bench_par}" "${log_file}" \
      >> output/benchmark/radsyn_benchmark_times.csv
    return 0
  fi

  end_ts="$(date +%s)"
  elapsed="$((end_ts - start_ts))"
  echo "-- spherical/u512: np=${np} failed with exit code ${exit_code}; see ${log_file}" >&2
  if [[ ! -f output/benchmark/radsyn_benchmark_failures.csv ]]; then
    printf 'suite,case,np,seconds,parfile,logfile,exit_code\n' > output/benchmark/radsyn_benchmark_failures.csv
  fi
  printf '%s,%s,%s,%s,%s,%s,%s\n' \
    "spherical" "u512" "${np}" "${elapsed}" "${bench_par}" "${log_file}" "${exit_code}" \
    >> output/benchmark/radsyn_benchmark_failures.csv
  return "${exit_code}"
}

overall_exit=0
for np_value in "${np_values[@]}"; do
  if ! run_one_np "${np_value}"; then
    overall_exit=1
  fi
done

exit "${overall_exit}"
