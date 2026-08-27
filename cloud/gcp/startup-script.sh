#!/usr/bin/env bash
# Runs as root on every VM boot (passed via --metadata-from-file startup-script=...).
# Handles infrastructure only (disk mount, base packages) - conda/mamba, cloning the
# repo, and creating the pipeline environment are manual first-login steps (see
# cloud/gcp/README.md), the same as the local setup in manual/requirements.md. This is
# deliberate: this script never touches your GitHub credentials or conda envs, so
# there's nothing sensitive baked into VM metadata.

set -euo pipefail
exec > >(tee -a /var/log/startup-script.log) 2>&1
echo "=== startup-script running at $(date) ==="

apt-get update -y
apt-get install -y git curl build-essential htop tmux

# Format and mount the attached data disk, idempotently. start_vm.sh attaches it with
# device-name=<DISK_NAME> and passes disk-name via instance metadata, so this always
# finds the right disk regardless of attachment order or how many other disks exist,
# even if DISK_NAME was customized in config.sh.
MOUNT_POINT=/mnt/data
DISK_NAME="$(curl -s -H 'Metadata-Flavor: Google' 'http://metadata.google.internal/computeMetadata/v1/instance/attributes/disk-name')"
DEVICE="/dev/disk/by-id/google-${DISK_NAME}"

if [[ -e "$DEVICE" ]]; then
    if ! blkid "$DEVICE" &> /dev/null; then
        echo "Formatting new data disk at $DEVICE..."
        mkfs.ext4 -F "$DEVICE"
    fi
    mkdir -p "$MOUNT_POINT"
    if ! mountpoint -q "$MOUNT_POINT"; then
        mount "$DEVICE" "$MOUNT_POINT"
    fi
    if ! grep -q "$MOUNT_POINT" /etc/fstab; then
        echo "$DEVICE $MOUNT_POINT ext4 defaults 0 2" >> /etc/fstab
    fi
    chmod 1777 "$MOUNT_POINT"
    echo "Data disk mounted at $MOUNT_POINT"
else
    echo "WARNING: no data disk found at $DEVICE - check the VM's attached-disks config" >&2
fi

echo "=== startup-script finished at $(date) ==="
