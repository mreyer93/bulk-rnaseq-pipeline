# Running on GCP

Use this when you want the `star_salmon` path (STAR genomic alignment, nf-core/rnaseq's
default), or when a batch is too large for a laptop. STAR needs ~30-40 GB RAM for a human
index, which is the whole reason the local config uses Salmon instead.

A spot VM plus a persistent data disk that outlives it: start it, run, copy results off,
stop it. Same shape as the variant-calling pipeline's cloud setup.

## Prerequisites

- GCP project with billing enabled and the Compute Engine API on
- [`gcloud` CLI](https://cloud.google.com/sdk/docs/install), authenticated (`gcloud init`)
- Edit `cloud/gcp/config.sh` and set `PROJECT_ID`

## One-time setup

```bash
./cloud/gcp/create_disk.sh    # persistent disk for references, FASTQs, results
./cloud/gcp/start_vm.sh       # spot VM, disk mounted at /mnt/data
gcloud compute ssh rnaseq-vm --project=<project> --zone=us-central1-a
```

On the VM (done interactively so your GitHub credentials never enter VM metadata):

```bash
git clone https://github.com/<you>/bulk-rnaseq-pipeline.git
cd bulk-rnaseq-pipeline

wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p /mnt/data/miniforge3
source /mnt/data/miniforge3/etc/profile.d/conda.sh
mamba env create -f envs/environment.yml -n bulk-rnaseq
mamba env create -f envs/r.yml -n bulk-rnaseq-r
conda activate bulk-rnaseq

# references onto the data disk (Ensembl human shown; use whatever matches your project)
mkdir -p /mnt/data/reference && cd /mnt/data/reference
REL=112
curl -O https://ftp.ensembl.org/pub/release-$REL/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
curl -O https://ftp.ensembl.org/pub/release-$REL/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
curl -O https://ftp.ensembl.org/pub/release-$REL/gtf/homo_sapiens/Homo_sapiens.GRCh38.$REL.gtf.gz
gunzip *.gz
```

Confirm the machine really is big enough before building the STAR index:

```bash
free -g          # want >= 40 GB total for a human index
nproc
```

## Routine workflow

```bash
./cloud/gcp/start_vm.sh
gcloud compute ssh rnaseq-vm --project=<project> --zone=us-central1-a

# on the VM
source /mnt/data/miniforge3/etc/profile.d/conda.sh && conda activate bulk-rnaseq
cd bulk-rnaseq-pipeline
python3 scripts/samplesheet.py /mnt/data/run1/samples.csv condition condition,treated,control
snakemake --configfile /mnt/data/run1/config.yaml --use-conda --cores 16

# copy results off before stopping
gsutil -m cp -r /mnt/data/run1/results gs://<your-bucket>/run1/

# back on your laptop
./cloud/gcp/stop_vm.sh
```

Build the STAR index **once**, then set `reference.star_index` in the config so later runs
skip it entirely — indexing a human genome takes about an hour.

### Preemption

Spot VMs can be reclaimed at any time. Restart with `start_vm.sh` and re-run the identical
`snakemake` command: completed outputs are kept and it resumes where it stopped. Add
`--rerun-incomplete` if a rule was interrupted mid-write.

## Cost

`us-central1`, verify against the [pricing calculator](https://cloud.google.com/products/calculator):

- **Compute**: `n2-standard-16` spot ≈ $0.23/hr, billed only while running.
- **Storage**: 500 GB `pd-balanced` ≈ $50/month, billed **whether or not the VM runs**.

For occasional use, storage dominates. Either tear down fully between projects (below),
or keep references in a GCS bucket (~4-5× cheaper per GB) and copy them onto a fresh disk
per run.

## Teardown

```bash
gcloud compute instances delete rnaseq-vm --project=<project> --zone=us-central1-a
gcloud compute disks delete rnaseq-data --project=<project> --zone=us-central1-a
```
