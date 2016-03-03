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


def schoolf_eq(Temps, B0, E, E_h, T_h):
	"""Full schoolfield model for calculating trait values at a given temperature"""
	
	#e <- math.exp(1)
	k = 8.617 * 10 ** (-5)
	
	model = B0 * np.exp(-E * ((1/(k*Temps)) - (1/(k*283.15)))) / 1+(np.exp(E_h * ((1/T_h * k)-(1/ k * Temps))))
	
	
	return model

#for i in data["unique_id"]:
#	if i == 1:
#		tmp_Datasubset = []
#		tmp_Datasubset.append(data.query('unique_id == i'))
#		print i


#a = data[(data.unique_id == 1)][["E"]]



def schoolf(params, Temps, Data_subset):
	"""Schoolfield model, to be called by schoolfield_model()"""
	
	
	B0 = params['B0'].value
	E = params['E'].value
	#E_l = params[''].value
	#T_l = params[''].value
	E_h = params['E_h'].value
	T_h = params['T_h'].value
	#temps = params['temps'].value
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	TraitVals = data["OriginalTraitValue"].values
	
	ModelPred = schoolf_eq(Temps, B0, E, E_h, T_h)
	#returns the residuals, needed to calculate the r square
	return(ModelPred - TraitVals)

def schoolfield_model(gg):
	"""NLLS fitting to the Schoolfield model; this function will 
	contain the lmfit.minimize calls to the schoolf() function. This is 
	where you can constrain the parameters."""

#	 Prepare the parameters and their bounds:

	#================================
	#data = pd.read_csv("Data/output.csv")
 	T_l_start = gg['T_l'].iloc[0]
	B0_start = gg['B0'].iloc[0]
	E_start = gg['E'].iloc[0]
	E_h_start = gg['E_h'].iloc[0]
	T_h_start = gg['T_h'].iloc[0]
	#temps = gg['Temps'][0]
	#=================================
 
	#	 The datase containing temperatures and trait values.
	#========================
	TraitVals = data["OriginalTraitValue"].values
	Temps = data["Temps"].values
	#====================================
	#TraitVals = data$TraitVals
	#Temps = ??

#	 Define parameters
#just a matter of how .... try and error thing?  TL, going to 280 290 k rreally
	params = Parameters()
	params.add('B0', value = B0_start, vary = True, min = -10, max = 1000)
	params.add('E', value=E_start, vary= True, min=0.0000000000000001, max=10)
	#params.add('E_l', value=E_l_start, vary = True, min=0.0000000000000001, max=10)
	params.add('E_h', value = E_h_start, vary = True, min=0.0000000000000001, max=10)
	params.add('T_h' ,value = T_h_start, vary = True, min=250, max=350)
	#params.add('T_l' ,value = T_l_start, vary = True, min=0.0000000000000001, max=10)
	
#	 Minimising the Model
	out = minimize(schoolf, params, args=(Temps, TraitVals),method="leastsq")
	par_out_school = out.params
#	 Calculates the r squared value to determine how well the model fits the data.
	r_squared_school = 1-out.residual.var()/sc.var(TraitVals)#variance of the resicudal, if logging the data, the out.residual will be logged       need to get the residuals and the sum of the residuals     the residueal 
	
	
	nvarys_school= out.nvarys


	ndata_school = out.ndata
	
	return(par_out_school, r_squared_school, nvarys_school, ndata_school,out.chisqr)
	

def AICrss(n, k, rss):
	"""Calculate the Akaike Information Criterion value, using:

	- n:   number of observations
	- k:   number of parameters
	- rss: residual sum of squares
	"""
	return n * log((2 * pi) / n) + n + 2 + n * log(rss) + 2 * k

def BICrss(n, k, rss):
	"""Calculate the Bayesian Information Criterion value, using:
	
	- n:   number of observations
	- k:   number of parameters
	- rss: residual sum of squares
	"""
	return n + n * log(2 * pi) + n * log(rss / n) + (log(n)) * (k + 1)

#~ ############################
#~ # M  A  I  N    C  O  D  E #
#~ ############################


def main():
	"""Performs fitting to the Gaussian-Gompertz, Schoolfield and Cubic model,
	and returns the best fits as a csv file to ./Results/results.csv"""
	#Produce an error is there is no dataset provided.
	
	
	
	
	#data = sc.genfromtxt(argv,dtype = None,delimiter = ',',deletechars='"')
	
	
	
	
	#input file "./Data/ThermResp_startvals.csv"
	# Define the Boltzmann constant (units of eV * K^-1).
	global k
	k = 8.617 * 10 ** (-5)
	data = pd.read_csv("Data/output.csv")
	id_pd = data['unique_id']
	global id_lsit
	id_list = id_pd.tolist()
	index_panda = data['index_panda']
	index_panda_list = index_panda.tolist()
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	global TraitVals 
	TraitVals = data["OriginalTraitValue"].values
	temp_kel = data['Temps']
	global temp_kel_list
	temp_kel_list = temp_kel.tolist()
	
	# Open csv to write parameter estimates to
	results = open("Results/results.csv", 'w')
	results_csv = csv.writer(results, delimiter=";")
	results_csv.writerow(['unique_id','B0','E','E_h','T_h','r_squared_school','nvarys_school','ndata_school','out.chisqr', 'Temps' ])
	
	qwe = id_list[-1]
	qwe = qwe+ 1
	
	for u in range(1, qwe):
		gg = data.loc[(data.unique_id == u)] [['FinalID','OriginalTraitValue','Temps','logB0','E','E_h','T_h','T_l','E_l','unique_id','index_panda','one_over_kt','B0']]
		dllm = data.loc [(data.unique_id == 1)] [[ 'Temps' ]]
		#gg = data[['FinalID','OriginalTraitValue','Temps','logB0','E','E_h','T_h','T_l','E_l','unique_id','index_panda','one_over_kt','B0']].query('data.unique_id == u')
		try:
			a = schoolfield_model(gg)
			results_csv.writerow( [id_list[u-1], a[0]['B0'].value, a[0]['E'].value ,a[0]['E_h'].value,a[0]['T_h'].value,a[1],a[2],a[3],a[4],temp_kel_list[u-1]])
			print 'dllm'
			
			
			
			
			
			
			
			
			
			#sum of square of the residual into the bic and aic
			#AIC and BIC
		except :
			print'gg'
			pass
			
	
	
	
	results.close()







#~ def main():
	#~ """Performs fitting to the Gaussian-Gompertz, Schoolfield and Cubic model,
	#~ and returns the best fits as a csv file to ./Results/results.csv"""
	#~ #Produce an error is there is no dataset provided.
	#~ 
	#~ 
	#~ 
	#~ 
	#~ #data = sc.genfromtxt(argv,dtype = None,delimiter = ',',deletechars='"')
	#~ 
	#~ 
	#~ 
	#~ 
	#~ #input file "./Data/ThermResp_startvals.csv"
	#~ # Define the Boltzmann constant (units of eV * K^-1).
	#~ global k
	#~ k = 8.617 * 10 ** (-5)
	#~ data = pd.read_csv("Data/output.csv")
	#~ id_pd = data['unique_id']
	#~ global id_lsit
	#~ id_list = id_pd.tolist()
	#~ index_panda = data['index_panda']
	#~ index_panda_list = index_panda.tolist()
	#~ #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	#~ global TraitVals 
	#~ TraitVals = data["OriginalTraitValue"].values

	#~ 
	#~ # Open csv to write parameter estimates to
	#~ results = open("Results/results.csv", 'w')
	#~ results_csv = csv.writer(results, delimiter=";")
	#~ results_csv.writerow(['unique_id','B0','E','E_h','T_h','r_squared_school','nvarys_school','ndata_school','out.chisqr' ])
	#~ 
	#~ qwe = id_list[-1]
	#~ qwe = qwe+ 1
	#~ 
	#~ for u in range(1, qwe):
		#~ gg = data.loc[(data.unique_id == u)] [['FinalID','OriginalTraitValue','Temps','logB0','E','E_h','T_h','T_l','E_l','unique_id','index_panda','one_over_kt','B0']]
		#~ #gg = data[['FinalID','OriginalTraitValue','Temps','logB0','E','E_h','T_h','T_l','E_l','unique_id','index_panda','one_over_kt','B0']].query('data.unique_id == u')
		#~ try:
			#~ a = schoolfield_model(gg)
			#~ results_csv.writerow( [id_list[u-1], a[0]['B0'].value, a[0]['E'].value ,a[0]['E_h'].value,a[0]['T_h'].value,a[1],a[2],a[3],a[4]])
			#~ print 'dllm'
			#~ 
			#~ 
			#~ 
			#~ 
			#~ 
			#~ 
			#~ 
			#~ 
			#~ 
			#~ #sum of square of the residual into the bic and aic
			#~ #AIC and BIC
		#~ except :
			#~ print'gg'
			#~ pass
			#~ 
	#~ 
	#~ 
	#~ 
	#~ results.close()









	
	#~ 
	#~ v = open("Results/results.csv", 'w')
	#~ r = csv.reader(v)
	#~ row0 = r.next()
	#~ row0.append('temp_kel')
	#~ for item in r:
		#~ item.append



	#==========================================================================================================================================		
	#~ df = pd.DataFrame()			
	#~ df = pd.DataFrame(index = index_panda_list ,columns=('unique_id','B0','E','E_h','T_h','r_squared_school','nvarys_school','ndata_school','out.chisqr'))
	#~ 
	#~ df=pd.DataFrame([1,a[0]['B0'].value,a[0]['E'].value,a[0]['E_h'].value,a[0]['T_h'].value,a[1],a[2],a[3],a[4]],columns=('unique_id','B0','E','E_h','T_h','r_squared_school','nvarys_school','ndata_school','out.chisqr'))
	#~ 
	
	
	
	#~ df_B0 = pd.DataFrame([a[0]['B0'].value], columns = 'B0')
	#~ df.append (df_B0,ignore_index=True)
	#~ df_E = pd.DataFrame([a[0]['E'].value], columns = 'E')
	#~ df.append (df_E,ignore_index=True)
	#~ df_E_h = pd.DataFrame([a[0]['E_h'].value], columns = 'E_h')
	#~ df.append (df_E_h,ignore_index=True)
	#~ df_T_h = pd.DataFrame([a[0]['T_h'].value], columns = 'T_h')
	#~ df.append (df_T_h,ignore_index=True)
	#~ df_r_suqared_school = pd.DataFrame([a[1]['r_suqared_school']], columns = 'r_suqared_school')
	#~ df.append (df_r_suqared_school,ignore_index=True)
	#~ df_nvarys_school = pd.DataFrame([a[2]['nvarys_school']], columns = 'nvarys_school')
	#~ df.append (df_nvarys_school,ignore_index=True)
	#~ df_ndata_school = pd.DataFrame([a[3]['ndata_school']], columns = 'ndata_school')
	#~ df.append (df_ndata_school,ignore_index=True)
	#~ df_out.chisqr = pd.DataFrame([a[4]['out.chisqr']], columns = 'out.chisqr')
	#~ df.append (df_out.chisqr,ignore_index=True) 
	#~ 
	
	#~ df.to_csv('../Results/result.csv')

	

#run school f model 




	#Open the csv file to write the output to.
	#??
	#results = open("../Results/results.csv", 'w')
	#results_csv = csv.writer(results, delimiter=",")


		
		
	# Here you will run the lmfitting over all unique data series. you 
	# will have to use try and except as the fitting won't always work with each model.

#~ if __name__ == "__main__":
	#~ #The input file name will be the minimum input, but you can add more inputs if you want 
	#~ main(sys.argv[0])	
