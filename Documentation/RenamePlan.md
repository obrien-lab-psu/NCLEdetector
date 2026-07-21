# Rename Plan: NCLEdetector -> NCLEdetector

## Goal
Rename the repository and package safely while keeping a rollback path and minimizing user disruption.

## Recommended Strategy
Use a two-phase approach:

1. Phase A: Git repository rename only (low risk)
2. Phase B: Code/package rename sweep (higher risk, done in a dedicated branch)

---

## Phase A: Commit, Push, and Rename the Git Repository

1. Confirm no active production jobs depend on this working tree.
2. Commit all current local changes.
3. Push to the current remote.
4. Rename the remote repository (GitHub UI or gh CLI).
5. Update local origin URL to the new remote path.
6. Verify fetch/push still works.
7. Rename local folder from NCLEdetector to NCLEdetector.
8. Re-open the renamed folder in VS Code.

Why first:
- This is fast and usually safe.
- You get a clean rollback point before broad text/path changes.

---

## Phase B: Rename Package/App References in a Dedicated Branch

Create a branch specifically for the rename sweep.

Scope:
- Directory rename: NCLEdetector -> NCLEdetector
- Import rename: NCLEdetector -> NCLEdetector
- Lowercase references: ncledetector -> ncledetector where appropriate
- Docs and config updates

Important caution:
- Do not run an unrestricted global sed over all files.
- Exclude binary/generated/artifact directories.
- Validate each replacement class (imports, package metadata, CLI names, docs links).

Suggested execution order:

1. Move package directory using git mv.
2. Update Python imports and internal references.
3. Update packaging metadata and project naming.
4. Update docs and scripts.
5. Run static checks and smoke tests.
6. Fix any broken entry points or module paths.
7. Push branch and open PR.

---

## High-Impact Files to Review During Rename

- setup.py
- pyproject.toml
- README.md
- NCLEdetector/__init__.py
- scripts/run_compare_sim2exp.py
- scripts/configs/workflow3_consistency_config.json
- container/Dockerfile

Also inspect:
- Documentation/
- assets/slurm/scripts/
- any shell wrappers, workflow configs, and installed entry points

---

## Compatibility Recommendation (Optional but Strongly Suggested)

For one release cycle, keep a compatibility shim package named NCLEdetector that re-exports from NCLEdetector.

Benefits:
- Prevents immediate breakage for existing user scripts/imports.
- Allows a smoother migration window.

---

## Validation Checklist After Rename

1. Editable install works.
2. CLI entry points resolve.
3. Core workflow script imports run.
4. Existing workflow configs still parse.
5. Docs/tutorial import examples are updated.
6. Container build path references still work.

---

## Rollback Plan

If issues appear:

1. Revert the rename branch (or reset to pre-rename commit).
2. Keep remote repo rename if desired, but defer package rename.
3. Re-attempt with narrower replacement scope and better test coverage.

---

## Notes

- Repository rename and package rename are independent tasks.
- Repository rename is easy; package rename is where most breakage risk lives.
- Treat the package rename as a normal feature migration with tests and review.
