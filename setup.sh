#!/bin/bash

cd ~/loopgap_ant_optimization

module purge
module load PUT_THE_NAME_HERE(openEMS)

source # put the venv path

python optimize.py -o