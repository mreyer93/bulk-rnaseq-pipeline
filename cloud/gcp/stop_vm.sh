#!/usr/bin/env bash
# Stops the VM (keeps the data disk and boot disk, just stops paying for compute).
# The data disk keeps costing its per-GB storage fee whether the VM is running or not -
# see cloud/gcp/README.md for full teardown if you want to stop that too.
#
# Usage: ./cloud/gcp/stop_vm.sh

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./config.sh

gcloud compute instances stop "$VM_NAME" --project="$PROJECT_ID" --zone="$ZONE"
echo "Stopped. Restart any time with ./cloud/gcp/start_vm.sh - your data disk is untouched."
