#!/usr/bin/env bash
# One-time setup: create the persistent data disk that will hold references/CADD
# databases/pipeline outputs. This disk survives independently of the VM - you can
# stop/delete the VM to save on compute cost without losing this data, and re-attach it
# to a new VM later.
#
# Usage: ./cloud/gcp/create_disk.sh

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./config.sh

if [[ "$PROJECT_ID" == "your-gcp-project-id" ]]; then
    echo "Edit cloud/gcp/config.sh and set PROJECT_ID first." >&2
    exit 1
fi

if gcloud compute disks describe "$DISK_NAME" --project="$PROJECT_ID" --zone="$ZONE" &> /dev/null; then
    echo "Disk $DISK_NAME already exists in $ZONE - nothing to do."
    exit 0
fi

echo "Creating ${DISK_SIZE_GB}G $DISK_TYPE disk '$DISK_NAME' in $ZONE..."
gcloud compute disks create "$DISK_NAME" \
    --project="$PROJECT_ID" \
    --zone="$ZONE" \
    --size="${DISK_SIZE_GB}GB" \
    --type="$DISK_TYPE"

echo "Done. Run ./cloud/gcp/start_vm.sh next to create a VM with this disk attached."
