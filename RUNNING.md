# Running the loop-gap optimization on Longleaf

One CMA-ES run over five geometry parameters of a PCB loop-gap antenna. Every design is
one openEMS FDTD solve: about 7 min on 4 cores for the base design, longer for large ones.
The job uses one node, 32 cores, 8 solves at a time and at most ~800 solves (about 100
generations; each takes as long as its slowest solve). The time limit is 2 days.
Disk: about 1.7 MB per solve in `runs/opt/` (~1.4 GB at most), plus the container image.

## Start

```bash
git clone --branch run1 https://github.com/EarthyYL/loopgap_ant_optimization.git
cd loopgap_ant_optimization
sbatch sim.sh                 # add -A <account> if yours isn't the default
```

The first job downloads the container image (`loopgap.sif`) into the clone; later jobs
reuse it. Nothing else needs installing: openEMS and Python are all inside the image.

## Watch

`cmaes_<jobid>.log` gets one line per generation:

```
gen   3   best     2.114   0/8 meet spec 2   sigma 0.231
```

- `best` above 0: no design meets the S11 spec yet. Below 0: some do, and `-best` is their B1.
- **If it stays at `0/8 meet spec 2` for 20 or more generations, the spec is probably out of
  reach. Cancel (`scancel <jobid>`) and tell Yao** rather than spending the allocation.
- Every design and its cost is in `data/cmaes.csv`.

## If the job dies (time limit, node failure, scancel)

Submit the same `sbatch sim.sh` again. The run has a fixed seed and every finished solve is
kept in `runs/opt/`, so it replays to where it stopped in minutes and carries on; at most one
generation of solves is lost.

Between resubmits, do **not** change the code (no `git pull`), `--popsize`, `--seed` or
`loopgap.sif`: any of those makes it a different run and the replay breaks.
Changing `--cpus-per-task` is fine.

If it fails at the same generation again, one design is crashing openEMS: send the log.

## When it finishes

The log ends with the best design and its B1, or "no design met spec 2". Send Yao
`data/cmaes.csv` and the `cmaes_*.log` files. Keep `runs/opt/` until then: the best designs
get re-checked from it.

## Smoke test (before the real run)

A short job on the interactive partition: one generation of 2 solves. With the fixed seed
these are two large designs (~880k cells), so allow up to an hour.

```bash
sbatch -p interact -c 8 --mem=8g -t 2:00:00 sim.sh --popsize 2 --threads 4 --maxfevals 1
```

(`--maxfevals 1` is one generation: pycma stops once the evaluations *exceed* it.)

The log should end with `stopped:` and a best design. Submit the same line again: it should
finish in seconds, with every solve read back from `runs/opt/`. That is the resume mechanism
working. Delete `runs/`, `data/cmaes.csv` and the logs before the real run.
