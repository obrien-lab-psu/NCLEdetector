# Using NCLEdetector via Containers (Docker and Apptainer/Singularity)

This guide shows how to run NCLEdetector without creating a local Conda environment.
Instead, pull the published container image from GHCR and run commands inside it.

---

## Why use the container

- No local dependency solving.
- Reproducible software stack across users and systems.
- Works well on HPC systems that support Apptainer or Singularity.

---

## Prerequisites

1. One container runtime is installed:
   - Docker, or
   - Apptainer/Singularity.
2. You can access GHCR (public image pull, or authenticated pull if required by your site policy).

Check your runtime:

```bash
docker --version
# or
apptainer version
# or
singularity version
```

---

## Pull the image (Docker)

From your project directory:

```bash
cd /path/to/NCLEdetector
docker pull ghcr.io/obrien-lab-psu/ncledetector:latest
```

If you want a specific released tag:

```bash
docker pull ghcr.io/obrien-lab-psu/ncledetector:vX.Y.Z
```

---

## Pull the image (Apptainer/Singularity)

From your project directory:

```bash
cd /path/to/NCLEdetector
apptainer pull -F ncledetector-latest.sif docker://ghcr.io/obrien-lab-psu/ncledetector:latest
```

If you want a specific released tag:

```bash
apptainer pull -F ncledetector-vX.Y.Z.sif docker://ghcr.io/obrien-lab-psu/ncledetector:vX.Y.Z
```

If you see "manifest unknown", the tag has not been published yet. Use latest or wait for the tag-triggered container workflow to finish.

---

## Optional GHCR login

For Docker, only needed if your organization/site policy requires authentication:

```bash
echo "YOUR_GITHUB_PAT" | docker login ghcr.io -u YOUR_GITHUB_USERNAME --password-stdin
```

For Apptainer/Singularity:

```bash
apptainer registry login --username YOUR_GITHUB_USERNAME docker://ghcr.io
```

Enter your GitHub PAT (with package read scope) when prompted.

---

## Validate the container

Docker:

```bash
docker run --rm ghcr.io/obrien-lab-psu/ncledetector:latest \
   bash /opt/NCLEdetector/container/smoke_test.sh
```

Apptainer/Singularity:

```bash
apptainer exec ncledetector-latest.sif bash /opt/NCLEdetector/container/smoke_test.sh
```

You should see: SMOKE TEST PASSED

---

## Run NCLEdetector scripts inside the container

Use bind mounts so your input/output paths are visible in the container.

Example using an HPC datastore with Docker:

```bash
DATASTORE=/scratch/ims86/NCLEdetector_Datastore
docker run --rm \
   -v "$DATASTORE:$DATASTORE" \
   -v "$PWD:$PWD" \
   -w "$PWD" \
   ghcr.io/obrien-lab-psu/ncledetector:latest \
   python scripts/run_nativeNCLE.py --help
```

Example using an HPC datastore with Apptainer/Singularity:

```bash
DATASTORE=/scratch/ims86/NCLEdetector_Datastore
apptainer exec \
  --bind "$DATASTORE:$DATASTORE" \
  --bind "$PWD:$PWD" \
  --pwd "$PWD" \
  ncledetector-latest.sif \
  python scripts/run_nativeNCLE.py --help
```

Run a workflow command from the repo root (Apptainer/Singularity):

```bash
apptainer exec \
  --bind "$DATASTORE:$DATASTORE" \
  --bind "$PWD:$PWD" \
  --pwd "$PWD" \
  ncledetector-latest.sif \
  python scripts/run_OP_on_simulation_traj.py --config scripts/configs/workflow2_MSM_config.json
```

---

## Common issues

1. manifest unknown
   - The requested GHCR tag does not exist yet.
   - Pull latest or confirm the GitHub Actions build-container workflow completed successfully.

2. 403 Forbidden on pull
   - Log in to GHCR with docker login or apptainer registry login.
   - Verify token scope and organization package permissions.

3. permission denied while trying to connect to the docker API
   - Your user may not have permission to access /var/run/docker.sock.
   - Use a system where Docker daemon access is enabled for your account.

4. File not found inside container
   - Add the correct --bind arguments.
   - Ensure your command paths are valid from the container working directory.

---

## When to still use Conda

Use local Conda if you are actively developing package code and frequently reinstalling dependencies. For stable analysis and reproducible runs, prefer the container.
