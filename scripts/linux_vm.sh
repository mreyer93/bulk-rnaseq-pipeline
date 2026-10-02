#!/usr/bin/env bash
# Run this pipeline inside a Linux virtual machine on an Apple Silicon Mac.
#
# The STAR paths (quantifier: star_salmon, and the TE analysis) need Linux: bioconda's
# Apple Silicon STAR exits successfully having read zero reads, because its input
# buffering relies on libstdc++ behaviour that macOS's libc++ does not have (STAR issues
# #2663 and #2142, no upstream fix). Everything else in the pipeline runs natively.
#
# This wraps Lima (https://lima-vm.io), which needs no administrator rights and uses
# Apple's Virtualization framework, so Linux runs at close to native speed. Your home
# directory is shared into the VM read-write at the same path, so the same config files,
# FASTQs and output directories work unchanged. Conda environments for the VM live on the
# VM's own disk, apart from the macOS ones.
#
# Usage:
#   scripts/linux_vm.sh setup                    # once: Lima, the VM, conda inside it
#   scripts/linux_vm.sh run [snakemake options]  # run snakemake in the VM from this directory
#   scripts/linux_vm.sh stop                     # stop the VM and give its memory back
#
# Settings, as environment variables: LIMA_VM (default bio), VM_CPUS (8),
# VM_MEMORY_GIB (40; a human or mouse STAR index needs about 32), VM_DISK_GIB (200).

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

VM="${LIMA_VM:-bio}"
LIMA_VERSION="v2.2.0"
export PATH="$HOME/.local/bin:$PATH"

[[ "$(uname -s)" == "Darwin" ]] || { echo "This is for macOS; on Linux run snakemake directly." >&2; exit 1; }
case "$PWD/" in "$HOME"/*) ;; *) echo "The repository must be under $HOME to be visible in the VM." >&2; exit 1 ;; esac

in_vm() { limactl shell "$VM" -- bash -lc "$1"; }

setup() {
    if ! command -v limactl >/dev/null 2>&1; then
        echo "==> Installing Lima $LIMA_VERSION into ~/.local"
        tmp=$(mktemp -d)
        base="https://github.com/lima-vm/lima/releases/download/$LIMA_VERSION"
        tarball="lima-${LIMA_VERSION#v}-Darwin-arm64.tar.gz"
        curl -fsSL -o "$tmp/$tarball" "$base/$tarball"
        curl -fsSL -o "$tmp/SHA256SUMS" "$base/SHA256SUMS"
        (cd "$tmp" && shasum -a 256 -c SHA256SUMS --ignore-missing)
        mkdir -p "$HOME/.local"
        tar -C "$HOME/.local" -xzmf "$tmp/$tarball"
        rm -rf "$tmp"
    fi
    if ! limactl list --quiet 2>/dev/null | grep -qx "$VM"; then
        echo "==> Creating VM $VM"
        limactl create --tty=false --name="$VM" --vm-type=vz --mount-type=virtiofs \
            --cpus="${VM_CPUS:-8}" --memory="${VM_MEMORY_GIB:-40}" --disk="${VM_DISK_GIB:-200}" \
            --containerd=none --mount-only "$HOME:w" template:ubuntu-24.04
    fi
    limactl start --tty=false "$VM" >/dev/null
    echo "==> Conda inside the VM"
    in_vm 'test -x ~/miniforge3/bin/conda || {
        curl -fsSL -o /tmp/miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh
        bash /tmp/miniforge.sh -b -p ~/miniforge3 >/dev/null; }'
    in_vm "cd '$PWD' && ~/miniforge3/bin/conda env update -q -n bulk-rnaseq -f envs/environment.yml >/dev/null"
    echo "VM $VM is ready. Run: scripts/linux_vm.sh run --configfile your_config.yaml --cores ${VM_CPUS:-8}"
}

run() {
    limactl start --tty=false "$VM" >/dev/null
    args=$(printf '%q ' "$@")
    in_vm "cd '$PWD' && export PATH=~/miniforge3/envs/bulk-rnaseq/bin:~/miniforge3/condabin:\$PATH && \
        snakemake -s workflow/Snakefile --use-conda --conda-frontend conda \
        --conda-prefix ~/snakemake-conda $args"
}

case "${1:-}" in
    setup) setup ;;
    run)   shift; run "$@" ;;
    stop)  limactl stop "$VM" ;;
    *)     sed -n '15,18p' "$0"; exit 1 ;;
esac
