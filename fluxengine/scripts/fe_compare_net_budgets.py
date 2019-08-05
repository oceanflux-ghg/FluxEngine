#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:50:09 2018

@author: tomholding
"""

from sys import argv;

from fluxengine.tools.lib_compare_net_budgets import calc_net_budget_percentages;

if __name__ == "__main__":
    if len(argv != 3):
        print("Usage: compare_net_budgets newOutputPath referenceOutputPath");
    else:
        outputPath = argv[1];
        referencePath = argv[2];    
        calc_net_budget_percentages(outputPath, referencePath);
