#!/usr/bin/env python3
from EntDetect._logging import setup_logger

"""
Bundle commands from a file into SLURM jobs and submit them.

Reads a plain-text file where each line is one shell command, groups them into
bundles of --bundle lines, writes one SLURM script per bundle to ./tmp/, and
submits each with sbatch.

Usage example:
    python scripts/job_utilities/submit_cmds.py \\
        --cmds   /path/to/cmds.txt \\
        --tag    sliceGE \\
        --bundle 1

Arguments:
    --cmds    Path to file containing one shell command per line   [required]
    --tag     Short job-name tag appended to each SLURM job name  [required]
    --bundle  Number of commands to group per job (default: 1)
    --log_level  Logging verbosity: DEBUG, INFO, WARNING, ERROR (default: INFO)
"""

# ---------------------------------------------------------------------------
# Hard-coded SLURM template
# Line index 1  -> #SBATCH -J  (job name, filled in per bundle)
# Last line     -> the bundled commands (filled in per bundle)
# ---------------------------------------------------------------------------
_SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 72:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=40GB
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err


cd /storage/group/epo2/default/ims86/git_repos/EntDetect
source ~/.bashrc
conda activate entdetect

set -euo pipefail
{cmd}
"""


##########################################################################################################################################################
def main(argv=None):

    ###---------------------------------------------------------------------------------------------------------
    import sys
    import os
    import argparse
    import subprocess
    import time

    start_time = time.time()
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description="Bundle shell commands into SLURM jobs and submit them.")
    parser.add_argument("--cmds",    type=str, required=True,
                        help="Path to file containing one shell command per line")
    parser.add_argument("--tag",     type=str, required=True,
                        help="Short job-name tag appended to each SLURM job name")
    parser.add_argument("--bundle",  type=int, default=1,
                        help="Number of commands to group per job (default: 1)")
    parser.add_argument("--log_level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Logging verbosity (default: INFO)")
    parser.add_argument("--tmpdir", type=str, default="./tmp",
                        help="Directory to save temporary SLURM scripts (default: ./tmp)")
    args = parser.parse_args(argv)
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    import logging
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logger = setup_logger('submit_cmds', outdir='./', ID=args.tag, log_level=log_level)
    logger.info(f'args: {args}')
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # Input validation
    if not os.path.exists(args.cmds):
        logger.error(f'ERROR: Command file not found: {args.cmds}')
        sys.exit(1)

    os.makedirs('./tmp', exist_ok=True)
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # Read and bundle commands
    with open(args.cmds, 'r') as fh:
        all_cmds = [line.rstrip('\n') for line in fh if line.strip()]

    logger.info(f'Total commands: {len(all_cmds)} | bundle size: {args.bundle}')

    bundles = [all_cmds[i:i + args.bundle] for i in range(0, len(all_cmds), args.bundle)]
    logger.info(f'Total SLURM jobs to submit: {len(bundles)}')
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # Write and submit one SLURM script per bundle
    for bundle_i, bundle in enumerate(bundles):
        job_name  = f'c{bundle_i}{args.tag}'
        cmd_block = '\n'.join(bundle)

        slurm_content = _SLURM_TEMPLATE.format(job_name=job_name, cmd=cmd_block)

        slurm_path = f'tmp/{bundle_i}{args.tag}.slurm'
        with open(slurm_path, 'w') as fh:
            fh.write(slurm_content)
        logger.info(f'SAVED: {slurm_path} -> submitting')

        # result = subprocess.run(['sbatch', slurm_path], capture_output=True, text=True)
        # if result.returncode == 0:
        #     logger.info(f'  {result.stdout.strip()}')
        # else:
        #     logger.error(f'  sbatch failed: {result.stderr.strip()}')
    ###---------------------------------------------------------------------------------------------------------

    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0
##########################################################################################################################################################


if __name__ == "__main__":
    raise SystemExit(main())

