#!/usr/bin/env python
""" This code performs non-linear least squares fitting of different
unimodal functions to experimental thermal response curves."""
__author__ = 'Your Name (your.email.id@imperial.ac.uk)'
__version__ = '0.0.1'

# Scroll down to the section called "MAIN CODE" first.
import sys
import numpy as np
from numpy import log, exp, pi
import scipy as sc 
from lmfit import Minimizer, minimize, Parameters, Parameter, report_fit, fit_report
import pandas as pd
import csv

global data
data = pd.read_csv("Data/output.csv")

#############################
# F  U  N  C  T  I  O  N  S #
#############################


def cubic_eq(Temps, B0, E, E_h, T_h):
	"""Full schoolfield model for calculating trait values at a given temperature"""
	
	#e <- math.exp(1)
	k = 8.617 * 10 ** (-5)
	
	model = B0 + (B1*T) + (B2*(T**2) + (B3*T**3)
	
	
	return model

def cubicf(params, Temps, Data_subset):
	"""Schoolfield model, to be called by schoolfield_model()"""
	
	
	B0 = params['B0'].value
	#~ E = params['E'].value
	#E_l = params[''].value
	#T_l = params[''].value
	#~ E_h = params['E_h'].value
	#~ T_h = params['T_h'].value
	#temps = params['temps'].value
	#~ TraitVals = data["OriginalTraitValue"].values
	#trait vale????????????????????????????????
	
	ModelPred = cubic_eq(Temps, B0, E, E_h, T_h)
	#returns the residuals, needed to calculate the r square
	return(ModelPred - TraitVals)
	


