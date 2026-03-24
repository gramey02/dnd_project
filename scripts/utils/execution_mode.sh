#!/bin/bash
#$ -N execution_mode
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o logs/out/execution_mode.out
#$ -e logs/err/execution_mode.err

set -euo pipefail

get_run_mode() {
    local run_mode="${RUN_MODE:-hpc}"
    run_mode="${run_mode,,}"

    if [[ "$run_mode" != "hpc" ]]; then
        echo "Error: RUN_MODE must be 'hpc', got '${RUN_MODE:-unset}'." >&2
        return 1
    fi

    printf '%s\n' "$run_mode"
}

run_indexed_jobs() {
    local count="$1"
    shift

    local target_script="$1"
    shift

    if (( count <= 0 )); then
        return 0
    fi

    local run_mode
    run_mode="$(get_run_mode)"

    if ! command -v qsub >/dev/null 2>&1; then
        echo "Error: RUN_MODE is 'hpc' but 'qsub' is not available in this shell." >&2
        return 1
    fi

    local output_root="${1:-$PWD}"
    local script_name
    script_name="$(basename "$target_script" .sh)"
    local log_dir="$output_root/logs/$script_name"
    mkdir -p "$log_dir"

    local -a qsub_args
    qsub_args=(-cwd -V -sync y -t "1-$count" -o "$log_dir" -e "$log_dir")

    if [[ -n "${HPC_QSUB_ARGS:-}" ]]; then
        local -a extra_qsub_args
        # shellcheck disable=SC2206
        extra_qsub_args=(${HPC_QSUB_ARGS})
        qsub_args+=("${extra_qsub_args[@]}")
    fi

    echo "Submitting ${count} HPC array tasks for ${script_name}..."
    qsub "${qsub_args[@]}" "$target_script" "$@"
}
