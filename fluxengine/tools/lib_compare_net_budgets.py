#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:50:09 2018

@author: tomholding
"""

from numpy import loadtxt;

#read the global flux budgets output.
def read_global_core_budgets(outputPath):
    budgets = loadtxt(outputPath, delimiter="'", dtype=str);
    colNames = str(budgets[0]).split(",")[2:11];
    budgetVals = [float(v) for v in budgets[-1].split(",")[2:11]];
    
    return (colNames, budgetVals);

#Compare global budgets for two (a new output and a reference) runs.
#Returns a dictionary of column names and percentage difference from reference.
def calc_net_budget_percentages(outputPath, referencePath, verbose=False):
    curBudgets = loadtxt(outputPath, delimiter="'", dtype=str);
    refBudgets = loadtxt(referencePath, delimiter="'", dtype=str);
    
    colNames = str(refBudgets[0]).split(",")[2:11];
    curGlobal = str(curBudgets[-1]).split(",")[2:11];
    refGlobal = str(refBudgets[-1]).split(",")[2:11];
    
    percentages = {};
    if verbose:
        print("Percentage difference between reference and new output:");
    for i, colName in enumerate(colNames):
        percent = float(refGlobal[i])/float(curGlobal[i]) * 100.0;
        if verbose:
            print("\t"+colName+": ", percent, "(%)");
        percentages[colName] = percent;
    
    return percentages;
