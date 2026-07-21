# NCLEdetector container image
#
# Reproducible image for running the NCLEdetector minimal-workflow command-line
# tools (run_nativeNCLE, run_OP_on_simulation_traj, run_MSM, ...) directly.
#
# Base is glibc-based (Ubuntu jammy) because the bundled pulchra/stride binaries
# are dynamically-linked Linux x86-64 ELF executables. Target platform is
# linux/amd64. (ARM hosts run it via emulation; HPC Singularity is amd64.)
#
# Build (from repo root):
#   docker build -f container/Dockerfile -t ghcr.io/obrien-lab-psu/ncledetector:local .
FROM condaforge/miniforge3:latest

# ---------------------------------------------------------------------------
# System dependencies
#   perl        -> calc_Q.pl / calc_K.pl (Q and K order parameters)
#   libgomp1 /
#   libstdc++6  -> runtime libs for the bundled pulchra / stride binaries
#   procps      -> convenience (ps) for debugging
# ---------------------------------------------------------------------------
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    perl procps ca-certificates libgomp1 libstdc++6 build-essential \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Create the exactly-pinned conda environment from the lock file
# ---------------------------------------------------------------------------
COPY container/environment.lock.yml /tmp/environment.lock.yml
RUN conda env create -n ncledetector -f /tmp/environment.lock.yml \
    && conda clean --all --yes \
    && rm -f /tmp/environment.lock.yml

# Put the env on PATH so console-scripts are directly callable in both
# `docker run` and `singularity exec` (Singularity honors the image ENV PATH).
ENV PATH=/opt/conda/envs/ncledetector/bin:$PATH \
    MAMBA_DEFAULT_ENV=ncledetector

# ---------------------------------------------------------------------------
# Install NCLEdetector itself from the checkout (deps already satisfied by lock)
# ---------------------------------------------------------------------------
COPY . /opt/NCLEdetector
WORKDIR /opt/NCLEdetector
RUN pip install --no-deps . \
    && chmod +x NCLEdetector/resources/pulchra NCLEdetector/resources/stride NCLEdetector/Jwalk/jwalk

# Default to a user-data working directory (bind-mount host data here or 1:1).
WORKDIR /work
ENTRYPOINT []
CMD ["bash"]
