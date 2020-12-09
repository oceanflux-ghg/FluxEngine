#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 01 09:55:10 2019

Launches the Jupyter notebook hub in the FluxEngine tutorials directory.

@author: Tom Holding
"""

import os;
import shutil;
from fluxengine.core.fe_setup_tools import get_fluxengine_root;

tutorialsFEDirectory = os.path.join(get_fluxengine_root(), "tutorials");
tutorialsCopyPath = os.path.join(os.getcwd(), "FluxEngineTutorials", "fe_jupyter_notebooks");
print("Copying tutorial files to:", tutorialsCopyPath);
try:
    shutil.copytree(tutorialsFEDirectory, tutorialsCopyPath);
except FileExistsError:
    print("*** Detected existing files at:", tutorialsCopyPath);
    print("*** Continuing using existing files at this path.");

print("Starting jupyter notebook...");
os.system("jupyter notebook --notebook-dir="+"\""+tutorialsCopyPath+"\"");