# NCLEdetector containers

Reproducible Docker and Singularity/Apptainer images that let you run the
NCLEdetector minimal-workflow command-line tools directly, without installing the
conda environment yourself.

The image bundles:

- the exactly-pinned `ncledetector` conda environment (see `environment.lock.yml`),
- the NCLEdetector package and its console-scripts,
- the native tools it shells out to: `pulchra`, `stride` (Linux x86-64 ELF),
  the `jwalk` Python tool, and `perl` for `calc_Q.pl` / `calc_K.pl`.

> Platform: **linux/amd64** (the pulchra/stride binaries are x86-64). HPC
> Singularity is amd64 → native. Apple-silicon Docker runs it via emulation.

## Console-scripts available in the image

`run_nativeNCLE`, `run_OP_on_simulation_traj`,
`run_nonnative_entanglement_clustering`, `run_MSM`, `run_compare_sim2exp`,
`run_population_modeling`, `run_montecarlo`, `run_Foldingpathway`,
`convert_cor_psf_to_pdb`.

## Pull the prebuilt image

```bash
# Docker
docker pull ghcr.io/obrien-lab-psu/ncledetector:latest

# Singularity / Apptainer (HPC)
apptainer pull ncledetector.sif docker://ghcr.io/obrien-lab-psu/ncledetector:latest
```

## Data + configs

You provide your own config JSON per run. Configs contain **absolute data
paths**, so the simplest pattern is to **bind-mount your datastore 1:1** — then
your existing configs work unchanged:

```bash
# Docker: mount the datastore at the same path it uses on the host
docker run --rm \
  -v /scratch/me/NCLEdetector_Datastore:/scratch/me/NCLEdetector_Datastore \
  ghcr.io/obrien-lab-psu/ncledetector:latest \
  run_OP_on_simulation_traj \
    --config /scratch/me/NCLEdetector_Datastore/configs/OP.json \
    --Traj 0 \
    --DCD /scratch/me/NCLEdetector_Datastore/traj/t0.dcd

# Apptainer/Singularity: same idea with --bind
apptainer exec --bind /scratch/me/NCLEdetector_Datastore \
  ncledetector.sif \
  run_MSM --config /scratch/me/NCLEdetector_Datastore/configs/MSM.json
```

(Optional portable alternative: bind your data at `/data` and write configs with
`/data/...` paths.)

## Build locally

From the **repo root** (the build context must be the repo root so the package
and bundled binaries are included):

```bash
docker build -f container/Dockerfile -t ghcr.io/obrien-lab-psu/ncledetector:local .

# Singularity from the local Docker image, or from the def file:
apptainer build ncledetector.sif container/apptainer.def
```

## Verify the image

```bash
docker run --rm ghcr.io/obrien-lab-psu/ncledetector:latest \
  bash /opt/NCLEdetector/container/smoke_test.sh
# or
apptainer exec ncledetector.sif bash /opt/NCLEdetector/container/smoke_test.sh
```

## Reproducibility

`environment.lock.yml` is an exact export of the working `ncledetector` environment
(conda + pip, fully version-pinned). Regenerate it with:

```bash
conda env export -n ncledetector \
  | sed '/^\s*-\s*ncledetector==/d; /^prefix:/d' > container/environment.lock.yml
```

## CI

`.github/workflows/container.yml` builds and pushes the Docker image to GHCR on
version tags (`v*`), then builds a `.sif` from that same image and pushes it to
GHCR as an OCI artifact — so both containers stay in lockstep.
