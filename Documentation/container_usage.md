# Using EntDetect via Container (Apptainer/Singularity)

This guide shows how to run EntDetect without creating a local Conda environment.
Instead, pull the published container image from GHCR and run commands inside it.

---

## Why use the container

- No local dependency solving.
- Reproducible software stack across users and systems.
- Works well on HPC systems that support Apptainer or Singularity.

---

## Prerequisites

1. Apptainer or Singularity is installed.
2. You can access GHCR (public image pull, or authenticated pull if required by your site policy).

Check your runtime:

```bash
apptainer version
# or
singularity version
```

---

## Pull the image

From your project directory:

```bash
cd /path/to/EntDetect
apptainer pull -F entdetect-latest.sif docker://ghcr.io/obrien-lab-psu/entdetect:latest
```

If you want a specific released tag:

```bash
apptainer pull -F entdetect-vX.Y.Z.sif docker://ghcr.io/obrien-lab-psu/entdetect:vX.Y.Z
```

If you see "manifest unknown", the tag has not been published yet. Use latest or wait for the tag-triggered container workflow to finish.

---

## Optional GHCR login

If your site requires authenticated pulls:

```bash
apptainer registry login --username YOUR_GITHUB_USERNAME docker://ghcr.io
```

Enter your GitHub PAT (with package read scope) when prompted.

---

## Validate the container

Run the built-in smoke test:

```bash
apptainer exec entdetect-latest.sif bash /opt/EntDetect/container/smoke_test.sh
```

You should see: SMOKE TEST PASSED

---

## Run EntDetect scripts inside the container

Use bind mounts so your input/output paths are visible in the container.

Example using an HPC datastore:

```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore
apptainer exec \
  --bind "$DATASTORE:$DATASTORE" \
  --bind "$PWD:$PWD" \
  --pwd "$PWD" \
  entdetect-latest.sif \
  python scripts/run_nativeNCLE.py --help
```

Run a workflow command from the repo root:

```bash
apptainer exec \
  --bind "$DATASTORE:$DATASTORE" \
  --bind "$PWD:$PWD" \
  --pwd "$PWD" \
  entdetect-latest.sif \
  python scripts/run_OP_on_simulation_traj.py --config scripts/configs/workflow2_MSM_config.json
```

---

## Common issues

1. manifest unknown
   - The requested GHCR tag does not exist yet.
   - Pull latest or confirm the GitHub Actions build-container workflow completed successfully.

2. 403 Forbidden on pull
   - Log in to GHCR with apptainer registry login.
   - Verify token scope and organization package permissions.

3. File not found inside container
   - Add the correct --bind arguments.
   - Ensure your command paths are valid from the container working directory.

---

## When to still use Conda

Use local Conda if you are actively developing package code and frequently reinstalling dependencies. For stable analysis and reproducible runs, prefer the container.
