#!/bin/bash
# The CMA-ES optimization as one Slurm job on one node.  Submit from the repo clone; see RUNNING.md.
# Resume after a kill: submit it again, unchanged.
#SBATCH --job-name=loopgap-cmaes
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# memory: one base-design solve peaks at ~200 MB; bigger discs mean more cells
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=cmaes_%j.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"   # the clone: runs/opt and data/cmaes.csv are written here

# Apptainer exists only on Longleaf's compute nodes, so the image is pulled here, once.
# --disable-cache: no second copy of it in the home quota.
[ -f loopgap.sif ] || apptainer pull --disable-cache loopgap.sif docker://ghcr.io/earthyyl/loopgap:9f5cdd4

# --cleanenv keeps host modules (anaconda) out of the container.
# popsize is fixed: another popsize is another run, and resume only replays the same run.
# Threads follow --cpus-per-task and don't change the run (with fewer than 8 CPUs pass --threads).
# Extra arguments go to optimize.py, and later ones win: the smoke test in RUNNING.md uses that.
apptainer exec --cleanenv loopgap.sif /opt/openEMS/venv/bin/python optimize.py \
    --popsize 8 --threads $(( SLURM_CPUS_PER_TASK / 8 )) "$@"
