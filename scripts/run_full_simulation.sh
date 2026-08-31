#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXPECTED_R="4.4.2"
CONFIG="config/main_final_v003.R"
STUDY="main_final_v003"
COMMAND="${1:-all}"
TASK_ID="${2:-}"
CPUS="${CPUS:-${SLURM_CPUS_PER_TASK:-2}}"
WORK_ROOT="${CMO_WORK_ROOT:-$ROOT/generated/main_final_v003_work}"
PROJECT="$WORK_ROOT/simulation"
OUTPUT="$ROOT/generated/full-simulation"

command -v Rscript >/dev/null 2>&1 || {
  echo "Rscript is required." >&2
  exit 69
}
command -v sha256sum >/dev/null 2>&1 || {
  echo "GNU sha256sum is required." >&2
  exit 69
}
command -v flock >/dev/null 2>&1 || {
  echo "flock is required for task locking." >&2
  exit 69
}
ACTUAL_R="$(Rscript --vanilla -e 'cat(paste(R.version$major, R.version$minor, sep="."))')"
test "$ACTUAL_R" = "$EXPECTED_R" || {
  echo "R $EXPECTED_R is required; found R $ACTUAL_R." >&2
  exit 69
}
case "$CPUS" in
  ''|*[!0-9]*) echo "CPUS must be an integer of at least 2." >&2; exit 64 ;;
esac
test "$CPUS" -ge 2 || {
  echo "CPUS must be an integer of at least 2." >&2
  exit 64
}

prepare() {
  local prepare_lock="$WORK_ROOT/.prepare.lock"
  local prepare_fd
  mkdir -p "$WORK_ROOT"
  exec {prepare_fd}>"$prepare_lock"
  flock "$prepare_fd"
  if test ! -d "$PROJECT"; then
    cp -R "$ROOT/simulation" "$PROJECT"
  else
    cmp -s "$ROOT/simulation/SOURCE_SHA256.txt" \
      "$PROJECT/SOURCE_SHA256.txt" || {
      echo "The work copy differs from the current locked source." >&2
      exit 65
    }
  fi
  (
    cd "$PROJECT"
    sha256sum -c --status SOURCE_SHA256.txt
  ) || {
    echo "The simulation work copy failed its source manifest." >&2
    exit 65
  }

  local design="$PROJECT/design/$STUDY"
  local required=(
    units.csv tasks.csv queue.rds cells.rds truth.csv seed_registry.csv
    seeds.rds config_snapshot.rds fingerprint.rds design_fingerprint.rds
  )
  local present=0
  for file in "${required[@]}"; do
    test -f "$design/$file" && present=$((present + 1))
  done
  if test "$present" -gt 0 && test "$present" -lt "${#required[@]}"; then
    if test ! -d "$PROJECT/results/raw/$STUDY" && \
       test ! -d "$PROJECT/status/$STUDY"
    then
      rm -rf "$design"
    else
      echo "A partial design exists beside run state; inspect it manually." >&2
      exit 65
    fi
  fi
  CMO_ROOT="$PROJECT" Rscript --vanilla \
    "$PROJECT/R/prepare_design.R" --config="$CONFIG"
  flock -u "$prepare_fd"
}

run_task() {
  local task="$1"
  case "$task" in
    ''|*[!0-9]*) echo "Task must be an integer from 1 through 48." >&2; exit 64 ;;
  esac
  test "$task" -ge 1 && test "$task" -le 48 || {
    echo "Task must be an integer from 1 through 48." >&2
    exit 64
  }
  prepare
  local status_directory="$PROJECT/status/$STUDY"
  local lock_file="$status_directory/$(printf 'task_%05d.lock' "$task")"
  local lock_fd
  mkdir -p "$status_directory"
  exec {lock_fd}>"$lock_file"
  flock -n "$lock_fd" || {
    echo "Task $task is already running." >&2
    return 73
  }
  while test ! -f "$PROJECT/status/$STUDY/$(printf 'task_%05d.done' "$task")"
  do
    set +e
    CMO_ROOT="$PROJECT" \
      SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-$CPUS}" \
      Rscript --vanilla "$PROJECT/R/run_task.R" \
      --config="$CONFIG" --task="$task"
    status=$?
    set -e
    if test "$status" -ne 0 && test "$status" -ne 75; then
      exit "$status"
    fi
    if test "$status" -eq 75 && test -n "${SLURM_JOB_ID:-}"; then
      flock -u "$lock_fd"
      return 75
    fi
  done
  flock -u "$lock_fd"
}

finalize() {
  test -d "$PROJECT" || {
    echo "Run the prepare and task stages first." >&2
    exit 66
  }
  CMO_ROOT="$PROJECT" Rscript --vanilla \
    "$PROJECT/R/validate_results.R" --config="$CONFIG"
  CMO_ROOT="$PROJECT" Rscript --vanilla \
    "$PROJECT/R/summarize_results.R" --config="$CONFIG"
  CMO_ROOT="$PROJECT" Rscript --vanilla \
    "$PROJECT/R/build_article_outputs.R" --config="$CONFIG"

  local summary="$PROJECT/results/summary/$STUDY"
  local article="$summary/article"
  Rscript --vanilla "$ROOT/scripts/build_final_tables.R" \
    "$summary" "$article" "$PROJECT/config/main_final_v003.R"

  local stage="$WORK_ROOT/stage"
  rm -rf "$stage"
  mkdir -p "$stage"
  for file in \
    table_01_design.csv \
    table_02_benchmark.csv \
    table_03_root_inference.csv \
    table_04_special_power.csv \
    table_05_tail_strata.csv \
    table_07_ablation.csv \
    figure_01_tangency_contraction.png \
    figure_02_two_root_power.png
  do
    cp "$article/$file" "$stage/$file"
  done
  rm -rf "$OUTPUT"
  mv "$stage" "$OUTPUT"
  python3 "$ROOT/scripts/verify.py" --compare-simulation "$OUTPUT"
  echo "FULL SIMULATION REPRODUCTION: PASSED"
}

case "$COMMAND" in
  prepare)
    test -z "$TASK_ID" || { echo "prepare takes no task number." >&2; exit 64; }
    prepare
    ;;
  task)
    run_task "$TASK_ID"
    ;;
  finalize)
    test -z "$TASK_ID" || { echo "finalize takes no task number." >&2; exit 64; }
    finalize
    ;;
  all)
    test -z "$TASK_ID" || { echo "all takes no task number." >&2; exit 64; }
    prepare
    for task in $(seq 1 48); do
      run_task "$task"
    done
    finalize
    ;;
  *)
    echo "Usage: scripts/run_full_simulation.sh {prepare|task N|finalize|all}" >&2
    exit 64
    ;;
esac
