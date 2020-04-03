#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 01 09:55:10 2019

Launches the Jupyter notebook hub in the FluxEngine tutorials directory.

@author: Tom Holding
"""

import os;
from fluxengine.core.fe_setup_tools import get_fluxengine_root;

tutorialsDirectory = os.path.join(get_fluxengine_root(), "tutorials");
os.system("jupyter notebook --notebook-dir="+tutorialsDirectory);