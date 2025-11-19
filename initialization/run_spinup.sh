#!/bin/bash
set -e

# Define log file
LOGFILE="run_log/run_log_$(date +%Y%m%d_%H%M%S).log"

# Define file
FILENAME="${1}"

# Start global timer
start_all=$(date +%s)

log_step() {
  local message="$1"
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $message" | tee -a "$LOGFILE"
}

run_step() {
  local cmd="$1"
  local description="$2" 

  log_step "Running: $description..."
  start=$(date +%s)

  # Run the command, redirect stdout and stderr into the log (and to console)
  if eval "$cmd" 2>&1 | tee -a "$LOGFILE"; then
    end=$(date +%s)
    log_step "$description completed in $((end - start)) seconds."
  else
    end=$(date +%s)
    log_step "ERROR: $description failed after $((end - start)) seconds."
    log_step "Script aborted."
    exit 1
  fi
}

log_step "Script started."

# Each step
run_step "Rscript spinup_treeinit.R '$FILENAME'" "Spinup merge script"

end_all=$(date +%s)
log_step "All scripts completed successfully in $((end_all - start_all)) seconds."
