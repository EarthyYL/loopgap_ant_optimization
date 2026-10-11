# Running the loop-gap optimization on Longleaf

These instructions will set up a CMA-ES optimization run over five design (geometry) parameters of our PCB loop-gap antenna. 

Each individual design is one openEMS FDTD solve. The entire job will use one node, 32 cores, 8 solves at a time and at most ~800 solves (about 100 generations). 
I measured the solve time on Longleaf at 3-10 min per unit, depending on the design's size. The per-generation time is set by the slowest solve in the set (10 min), so expect the task to run just a bit under a day; the time limit is 2 days.

## 1. Get the code

Make sure to clone into the /users/ directory, this script has a lot of heavy I/O. Ideally, we would clone into ~/ and run targeted at /users/..., but this is much simpler.

```bash
cd /hickory/users/<o>/<n>/<onyen>   # your work directory
git clone --branch run1 https://github.com/EarthyYL/loopgap_ant_optimization.git
cd loopgap_ant_optimization
```

The first job downloads the container image (`loopgap.sif`) into the clone and later jobs
reuse it. Since it's containerized, you should not need to install anything else.

## 2. Smoke test (before the real run)

A short job: one generation of 2 solves. Seed 15's first two designs are small, so the solves
take about 5 min. The first job in a fresh clone also downloads the image (a few min).

```bash
# may need to add -A <account> if yours isn't the default
sbatch -c 8 --mem=8g -t 0:30:00 sim.sh --popsize 2 --threads 4 --maxfevals 1 --seed 15
```

The job log, `cmaes_<jobid>.log`, should end with a `every design: …` line. After the job completes, submit again (with lower resource asks). It should finish in seconds, as both solves will get reread from cached data:

```bash
sbatch -c 8 --mem=2g -t 0:05:00 sim.sh --popsize 2 --threads 4 --maxfevals 1 --seed 15
```

Delete `runs/`, `data/cmaes.csv` and the `cmaes_<jobid>.log` files if this all runs well, and proceed to step 3.

## 3. Start the real run

```bash
sbatch sim.sh                
```

## 4. Wait for it to finish

The run is hands-free: it stops by itself, after at most ~100 generations or earlier once the search converges, within a day or two. It's done when the job is gone from `squeue -u $USER` AND `cmaes_<jobid>.log` ends with an `every design: …` line (or `sacct -j <jobid> -o State` shows `COMPLETED`). If the second set of conditions is not fulfilled, then it is failed, not done. 

Check the log file after starting to verify that the container was setup correctly (`INFO:    Creating SIF file...` is the last line of the Apptainer setup). It will take a few minutes, then you should be good to let it run. I would appreciate it if you can check back in every once in a while to check for completion/failure.

If you want to peek at the work while it's running, `cmaes_<jobid>.log` is updated per-generation, showing progress in our objectives. Individual designs and results are kept in the `runs/opt/` folder, and `data/cmaes.csv` lists every design with its cost.

## 5. If the job dies (time limit, node failure, scancel)

Submit the same `sbatch sim.sh` again. The run has a fixed seed and every finished solve is kept in `runs/opt/`, so it will quickly replay back to where it stopped, without needing to do full solves.

Between resubmits, do **not** change the code, it will break.

If it fails at the same generation again, one design is crashing openEMS: don't resubmit, pack up the run as below and send it.

## 6. When it finishes (or if something goes wrong)

Pack up the whole run and send it to me:

```bash
cd /hickory/users/<o>/<n>/<onyen>
tar czf loopgap_results.tar.gz --exclude=loopgap.sif loopgap_ant_optimization
```

This includes the logs, `data/cmaes.csv` and every solve in `runs/opt/` (about 1.5 GB at
most).

## TLDR:
1. Grab the repo and container image
2. Run the smoke test
3. Run the simulation
    - Check in every once in a while to make sure it's alive.
     - If it died prematurely (no `every design: …` line at the end of the log, or `sacct` shows `CANCELLED`, `TIMEOUT`, `FAILED` or `NODE_FAIL`), resubmit it: see step 5. 
4. Regardless of success, once you are done, pack up the files and send them to me.
