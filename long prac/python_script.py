#!/usr/bin/env python

""" Perform NLLS fitting of Phytoplankton performance curves...  """

from lmfit import minimize, Parameters, Parameter, report_fit
import scipy as sc
from matplotlib import pylab as pl
import pandas as pd
import csv
import numpy as np

# You should figure out the imports you need!

#####################################################
# S C H O O L F I E L D   F  U  N  C  T  I  O  N  S #
#####################################################
	
def schoolf(params, temps, traits):
	"""Schoolfield model, to be called by schoolfield_model()"""
	
	B0 = params['B0'].value
	E = params['E'].value
	E_D = params['E_D'].value
	T_pk = params['T_pk'].value
	
	model = B0 * np.exp(-E * ((1/(K*temps)) - (1/(K*283.15)))) \ ...?
	
	model_log = np.array(map(np.log, model), dtype=np.float64)
	
	return model_log - np.log(traits)


def schoolfield_model(temps, B0_start, E_start, T_pk_start, E_D_start, traits):
	"""NLLS fitting to the Schoolfield model; this function will 
	contain the lmfit.minimize calls to the schoolf() function. This is 
	where you can constrain the parameters."""
	
	
	params = Parameters()
	params.add('B0', value = B0_start)
	params.add('E', value = E_start, min = ?, max = ?)
	params.add('E_D', value = E_D_start, min = ?, max = ?)
	params.add('T_pk', value = T_pk_start, min = ?, max = ?)
	
	results_sf = minimize(schoolf, params, args=(temps, traits))
	
	return results_sf


########################
# F O R M A T  D A T A #
########################		

data = pd.DataFrame.from_csv("Data/output.csv")
subset = 

	for u in UniqueID:

		# Get the subset of the data for this ID, in appropriate format:
		
		tmp_Datasubset = []
		
		for i in original_dataset:
		
			if i[0] == u:
				tmp_Datasubset.append(data[unique_id == u] [['logB0','E','E_h','E_I','T_h','T_I']])
			else:
				pass

####
id_pd = data['unique_id']
id_list = id_pd.tolist()

	for u in id_list:
		gg[] = data.loc[(data.unique_id == u)] [['logB0','E','E_h','T_h','T_l','E_l','unique_id','index_panda','one_over_kt']]



	#for u in id_list:
		
		
		#tmp=[]
		#if u == 1:
			#tmp.append(data[data.unique_id == '1'] [['logB0','E','E_h','E_I','T_h','T_I']])
			

		# Get the subset of the data for this ID, in appropriate format:
		
		tmp_Datasubset = []
		
		#for i in data:
		
			#if i[0] == u:
			if u == 1:
				#print i
				#tmp_Datasubset.append(i)
				tmp_Datasubset.append(data[unique_id == u] [['logB0','E','E_h','E_I','T_h','T_I']])
			else:
				pass
				
		
		
		#==================================================
		T_pk_start = gg['T_h'][0]
		B0_start = gg['logB0'][0]
		E_start = gg['E'][0]
		E_D_start = gg['E_h'][0]	
		#==================================================
		
							
		# Extract the starting values for the NLLS fitting:
		
		T_pk_start = float(tmp_Datasubset[1][10])
		B0_start = float(tmp_Datasubset[1][11])
		E_start = float(tmp_Datasubset[1][12])
		E_D_start = float(tmp_Datasubset[1][13])
		
		# Extract temperature and  trait data
		temps = []
		traits = []
		
		for i in tmp_Datasubset:
			temps.append(float(i[5]))
			traits.append(float(i[6]))
		
		# Convert list to array
		temps = np.array(temps)
		traits = np.array(traits)
		
		###########################################
		# F I T  S C H O O L F I E L D  M O D E L #
		###########################################
	   
		try:
			# Fit model
			Schoolf_fit = schoolfield_model(temps, B0_start, E_start, T_pk_start, E_D_start, traits)			

