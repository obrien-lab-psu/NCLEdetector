#!/usr/bin/env bash
# Smoke test for the EntDetect container image.
#
# Verifies the environment and bundled native tools are functional. Intended to
# run *inside* the container, e.g.:
#   docker run --rm ghcr.io/obrien-lab-psu/entdetect:latest bash /opt/EntDetect/container/smoke_test.sh
#   apptainer exec entdetect.sif bash /opt/EntDetect/container/smoke_test.sh
set -euo pipefail

echo "== EntDetect container smoke test =="

echo "-- console-scripts on PATH --"
for tool in run_nativeNCLE run_OP_on_simulation_traj \
            run_nonnative_entanglement_clustering run_MSM \
            run_compare_sim2exp run_population_modeling \
            run_montecarlo run_Foldingpathway convert_cor_psf_to_pdb; do
    command -v "$tool" >/dev/null || { echo "MISSING: $tool"; exit 1; }
    "$tool" --help >/dev/null 2>&1 || true
    echo "  ok: $tool"
done

echo "-- perl (calc_Q.pl / calc_K.pl) --"
perl --version >/dev/null
echo "  ok: perl"

echo "-- bundled native binaries --"
PKG=$(python -c "import EntDetect, os; print(os.path.dirname(EntDetect.__file__))")
"$PKG/resources/pulchra" -h >/dev/null 2>&1 || true   # prints usage / exits nonzero, but must be executable
"$PKG/resources/stride"  -h >/dev/null 2>&1 || true
test -x "$PKG/resources/pulchra" && echo "  ok: pulchra executable"
test -x "$PKG/resources/stride"  && echo "  ok: stride executable"

echo "-- python stack imports --"
python - <<'PY'
import EntDetect, topoly, mdtraj, MDAnalysis, freesasa, pyemma, numba, statsmodels
print("  ok: EntDetect, topoly, mdtraj, MDAnalysis, freesasa, pyemma, numba, statsmodels")
PY

echo "== SMOKE TEST PASSED =="
