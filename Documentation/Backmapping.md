# Back-map coarse-grained trajectories to all-atom structures

> **This step is computationally expensive and has already been completed** for the 1ZMR example system. The all-atom DCD trajectories in `$AA_TRAJ_DIR` are the result. This section describes the procedure for reference and for running on your own system.

Back-mapping reconstructs full all-atom models from Cα-only (coarse-grained) trajectory frames using the native all-atom structure as a template.

### 3a. Save individual Cα PDB frames from a trajectory

The frame extraction step is system-specific and typically done with MDAnalysis or VMD. Example with MDAnalysis:

```python
import MDAnalysis as mda

DATASTORE   = "/scratch/ims86/NCLEdetector_Datastore"
CG_TRAJ_DIR = f"{DATASTORE}/user_input/cg_trajectories"
OUTDIR      = f"{DATASTORE}/outputs/workflow3"

import os
os.makedirs(f"{OUTDIR}/cg_frames", exist_ok=True)

psf = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean_ca.psf"
dcd = f"{CG_TRAJ_DIR}/420_prod.dcd"

u = mda.Universe(psf, dcd)
protein = u.select_atoms("protein")

for i, ts in enumerate(u.trajectory[-335:]):   # last 335 frames
    out_pdb = f"{OUTDIR}/cg_frames/frame_{i:04d}.pdb"
    protein.write(out_pdb)
```

### 3b. Back-map each Cα frame to all-atom

```python
from NCLEdetector.change_resolution import BackMapping

DATASTORE = "/scratch/ims86/NCLEdetector_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow3"

# cg_pdb : a single CG (Cα-only) PDB frame extracted above
cg_pdb  = f"{OUTDIR}/cg_frames/frame_0000.pdb"
aa_pdb  = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb"
ID      = "1ZMR"

backMapper = BackMapping(outdir=f"{OUTDIR}/BackMapping")
backMapper.backmap(cg_pdb=cg_pdb, aa_pdb=aa_pdb, ID=ID)
```

### What back-mapping produces

Depending on configuration, outputs include:

- reconstructed all-atom structures (`.pdb`);
- intermediate PD2 / Pulchra outputs;
- OpenMM energy-minimization logs and energy-minimized final structures.

> **After back-mapping:** Inspect a representative subset of reconstructed structures in VMD before proceeding. Collate the validated per-frame PDBs into an all-atom DCD for downstream use.
