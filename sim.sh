#!/bin/bash
# The CMA-ES optimization as one Slurm job on one node.  Submit from the repo clone; see README.md.
# Resume after a kill: submit it again, unchanged.
#SBATCH --job-name=loopgap-cmaes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# memory: one base-design solve peaks at ~200 MB; bigger discs mean more cells
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=cmaes_%j.log
set -euo pipefail
# The clone, which should be in the work directory (/hickory/...), not home: runs/opt (heavy I/O),
# data/cmaes.csv, the logs and loopgap.sif are all written here.  Apptainer mounts it by default.
cd "$SLURM_SUBMIT_DIR"

# Apptainer exists only on Longleaf's compute nodes, so the image is pulled here, once.
# --disable-cache: no second copy of it in the home quota.
# By digest, not by tag (:9f5cdd4): a tag can be moved to another image, a digest can't.
# Rebuilt the image?  Put its new digest here in the same commit.
[ -f loopgap.sif ] || apptainer pull --disable-cache loopgap.sif \
    docker://ghcr.io/earthyyl/loopgap@sha256:6a50f82fdf6ab7548a0fb7741b3cd90b7d33ffaf00736e5b984cdb2e49f77fbe

# --cleanenv keeps host modules (anaconda) out of the container.
# popsize is fixed: another popsize is another run, and resume only replays the same run.
# Threads follow --cpus-per-task and don't change the run (with fewer than 8 CPUs pass --threads).
# Extra arguments go to optimize.py, and later ones win: the smoke test in README.md uses that.
apptainer exec --cleanenv loopgap.sif /opt/openEMS/venv/bin/python optimize.py \
    --popsize 8 --threads $(( SLURM_CPUS_PER_TASK / 8 )) "$@"
