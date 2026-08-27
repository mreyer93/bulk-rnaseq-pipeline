#!/usr/bin/env bash
# Creates the VM (first run) or starts it again (if it already exists but is stopped),
# with the data disk from create_disk.sh attached. Uses a spot VM by default - much
# cheaper, but GCP can reclaim it at any time; Snakemake resumes cleanly from wherever
# it left off if that happens, just re-run the same snakemake command after restarting.
#
# Usage: ./cloud/gcp/start_vm.sh

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./config.sh

if [[ "$PROJECT_ID" == "your-gcp-project-id" ]]; then
    echo "Edit cloud/gcp/config.sh and set PROJECT_ID first." >&2
    exit 1
fi

if ! gcloud compute disks describe "$DISK_NAME" --project="$PROJECT_ID" --zone="$ZONE" &> /dev/null; then
    echo "Data disk '$DISK_NAME' doesn't exist yet - run ./cloud/gcp/create_disk.sh first." >&2
    exit 1
fi

if gcloud compute instances describe "$VM_NAME" --project="$PROJECT_ID" --zone="$ZONE" &> /dev/null; then
    STATUS=$(gcloud compute instances describe "$VM_NAME" --project="$PROJECT_ID" --zone="$ZONE" --format='value(status)')
    if [[ "$STATUS" == "RUNNING" ]]; then
        echo "VM $VM_NAME is already running."
    else
        echo "Starting existing VM $VM_NAME (was $STATUS)..."
        gcloud compute instances start "$VM_NAME" --project="$PROJECT_ID" --zone="$ZONE"
    fi
else
    echo "Creating VM $VM_NAME ($MACHINE_TYPE, spot) in $ZONE with $DISK_NAME attached..."
    gcloud compute instances create "$VM_NAME" \
        --project="$PROJECT_ID" \
        --zone="$ZONE" \
        --machine-type="$MACHINE_TYPE" \
        --image-family="$IMAGE_FAMILY" \
        --image-project="$IMAGE_PROJECT" \
        --boot-disk-size="${BOOT_DISK_SIZE_GB}GB" \
        --disk="name=$DISK_NAME,mode=rw,boot=no,auto-delete=no,device-name=$DISK_NAME" \
        --provisioning-model=SPOT \
        --instance-termination-action=STOP \
        --metadata="disk-name=$DISK_NAME" \
        --metadata-from-file=startup-script=./startup-script.sh
fi

echo ""
echo "Once it's up (may take ~1min for the startup script to finish formatting/mounting the disk):"
echo "  gcloud compute ssh $VM_NAME --project=$PROJECT_ID --zone=$ZONE"
echo ""
echo "First time logging in, see cloud/gcp/README.md for the one-time environment setup"
echo "(clone the repo, create the conda env, download references onto /mnt/data)."
