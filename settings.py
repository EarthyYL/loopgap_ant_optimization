# -*- coding: utf-8 -*-
"""Where runs and results go, and how long a solve may run.  Values only."""

import os

repo_root = os.path.dirname(os.path.abspath(__file__))
Runs = os.path.join(repo_root, 'runs')
Data = os.path.join(repo_root, 'data')
Boards = os.path.join(repo_root, 'boards')   # KiCad projects, tracked in git

sim_time = 30e-9   # raise it if EndCriteria isn't reached
