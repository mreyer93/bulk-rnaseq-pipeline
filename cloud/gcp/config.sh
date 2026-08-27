#!/usr/bin/env bash
# Shared settings for the cloud/gcp/*.sh scripts. Edit these, then run create_disk.sh
# and start_vm.sh. The other scripts source this file automatically.

PROJECT_ID="${PROJECT_ID:-your-gcp-project-id}"
ZONE="${ZONE:-us-central1-a}"

VM_NAME="${VM_NAME:-rnaseq-vm}"
DISK_NAME="${DISK_NAME:-rnaseq-data}"

# 16 vCPU / 64 GB. The binding constraint is STAR: a human genome index needs ~30-40 GB
# RAM to build and load, which is exactly why the local config uses Salmon instead.
# n2-highmem-8 (8 vCPU / 64 GB) is a cheaper alternative if you are RAM- not CPU-bound.
MACHINE_TYPE="${MACHINE_TYPE:-n2-standard-16}"

# Data disk: references + STAR index (~35 GB for human) + FASTQs + results.
# RNA-seq FASTQs are the bulk of this - size for your largest expected batch.
DISK_SIZE_GB="${DISK_SIZE_GB:-500}"
DISK_TYPE="${DISK_TYPE:-pd-balanced}"

BOOT_DISK_SIZE_GB="${BOOT_DISK_SIZE_GB:-50}"
IMAGE_FAMILY="${IMAGE_FAMILY:-debian-12}"
IMAGE_PROJECT="${IMAGE_PROJECT:-debian-cloud}"
