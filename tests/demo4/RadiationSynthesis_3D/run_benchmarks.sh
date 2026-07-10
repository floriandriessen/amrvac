#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  tests/demo4/RadiationSynthesis_3D/run_benchmarks.sh

Environment:
  CORE_LIST   MPI sizes to test. Default: "1 2 4 8 16 32 64"
  INIT_NP     MPI size used to create missing snapshots. Default: 1
  RUN_INIT    if-missing | always | skip. Default: if-missing
  MPIRUN      MPI launcher. Default: mpirun
  MAKE        Make command. Default: make
  AMRVAC_DIR  AMRVAC repository root. Default: inferred from this script.

Examples:
  CORE_LIST="1 2 4 8" tests/demo4/RadiationSynthesis_3D/run_benchmarks.sh
  RUN_INIT=always INIT_NP=8 tests/demo4/RadiationSynthesis_3D/run_benchmarks.sh
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

core_list="${CORE_LIST:-1 2 4 8 16 32 64}"
init_np="${INIT_NP:-1}"
run_init="${RUN_INIT:-if-missing}"
mpirun_cmd="${MPIRUN:-mpirun}"
make_cmd="${MAKE:-make}"

case "${run_init}" in
  if-missing|always|skip) ;;
  *) echo "RUN_INIT must be one of: if-missing, always, skip" >&2; exit 2 ;;
esac

case_names=(u128 u256 u512 amr_l3 amr_l4 amr_l3_zstretch)

cd "${script_dir}"
"${make_cmd}" -f test.make setup
mkdir -p output/benchmark/logs
if [[ ! -f output/benchmark/radsyn_benchmark_times.csv ]]; then
  printf 'suite,case,np,seconds,parfile,logfile\n' > output/benchmark/radsyn_benchmark_times.csv
fi

for case_name in "${case_names[@]}"; do
  init_par="par_init/${case_name}.par"
  bench_par="par_benchmark/bench_${case_name}_corner45_thick.par"
  snapshot="output/radsyn_tdm_prom_${case_name}_0000.dat"

  if [[ "${run_init}" == "always" || ( "${run_init}" == "if-missing" && ! -f "${snapshot}" ) ]]; then
    init_log="output/benchmark/logs/init_${case_name}_np${init_np}.log"
    echo "-- cart/${case_name}: creating ${snapshot} with np=${init_np}"
    "${mpirun_cmd}" -np "${init_np}" ./amrvac -i "${init_par}" > "${init_log}" 2>&1
  elif [[ "${run_init}" == "skip" && ! -f "${snapshot}" ]]; then
    echo "Missing ${snapshot}; run without RUN_INIT=skip or create it manually." >&2
    exit 1
  fi

  for np in ${core_list}; do
    log_file="output/benchmark/logs/cart_${case_name}_np${np}.log"
    echo "-- cart/${case_name}: convert np=${np}"
    start_ts="$(date +%s)"
    "${mpirun_cmd}" -np "${np}" ./amrvac -i "${bench_par}" > "${log_file}" 2>&1
    end_ts="$(date +%s)"
    elapsed="$((end_ts - start_ts))"
    printf '%s,%s,%s,%s,%s,%s\n' \
      "cart" "${case_name}" "${np}" "${elapsed}" "${bench_par}" "${log_file}" \
      >> output/benchmark/radsyn_benchmark_times.csv
  done
done
