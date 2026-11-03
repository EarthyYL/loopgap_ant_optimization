#!/bin/bash
# Run 1's CMA-ES search as one Slurm job.  From the clone's root: sbatch opt/optimize.sh
# (docs/run1-longleaf.md has run 1's steps, from when this was sim.sh).  Resume after a kill: submit it again.
#SBATCH --job-name=loopgap-cmaes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=cmaes_%j.log
set -euo pipefail
# the clone belongs in the work directory (/hickory/...), not home: runs/opt is heavy I/O
cd "$SLURM_SUBMIT_DIR"

# Apptainer exists only on compute nodes, so pull here, once, by digest (a tag can move), and
# --disable-cache keeps a second copy out of the home quota.  Rebuilt the image?  New digest here.
[ -f loopgap.sif ] || apptainer pull --disable-cache loopgap.sif \
    docker://ghcr.io/earthyyl/loopgap@sha256:6a50f82fdf6ab7548a0fb7741b3cd90b7d33ffaf00736e5b984cdb2e49f77fbe

# --cleanenv keeps host modules out.  popsize is fixed: another popsize is another run, and
# resume only replays the same run.  Extra arguments go to opt/optimize.py; later ones win.
apptainer exec --cleanenv loopgap.sif /opt/openEMS/venv/bin/python -m opt.optimize \
    --popsize 8 --threads $(( SLURM_CPUS_PER_TASK / 8 )) "$@"
