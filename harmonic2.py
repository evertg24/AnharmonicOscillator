# major reoverhaual

import numpy as np
from matplotlib.pyplot import *
rcParams["font.family"] = "Times New Roman"
style.use("ggplot")
from matplotlib.widgets import Slider
import threading
import pandas as pd 
import concurrent.futures



v_pot = lambda x, λ: (0.5 * x**2 + λ * x**4) # λ=0 for regular oscillator


class AnharmonicOscillator:
		
	def __init__(self, N, 
		nEigState=1, xlimit=6, nEigVectors=3, λ_final = 10., λ_increment = 0.5): 
		'''xlimit is the xrange; nEigVectors == length2'''
		self.N, self.nEigState, self.nEigVectors = N, nEigState, nEigVectors # 
		self.x, self.λ_0 = np.linspace(-1*xlimit, xlimit, self.N - 2), 0 #dies off at -6, 6; λ = 0 for regular harmonic
		self.h, self.λ_final, self.λ_increment = self.x[1] - self.x[0], λ_final, λ_increment # 

		self.T, self.V = self.generateEmpty_T_V_arrays() #initalize empty arrays
		self.diagnolizeKineticArray()
		self.T = -self.T / (2 * self.h**2)
		self.H = self.T + self.V
		self.sweepThroughLambdas() # this gives the calculated pointz
		#self.eigConntinuationFit()

		self.v = np.array(self.eigVectors)# v is the traditional notation for GS, v = collection of lowest eigvects at λ 
		# gram-shmitt 
		
		self.GS_dool() # GS on the vectors
		self.eigvectorFit, self.Hii = [], np.zeros(nEigVectors**2).reshape(nEigVectors, nEigVectors) # 
		self.eigenVectorContniuationFit() # do the fit
	def eigenVectorContniuationFit(self):
		self.λ = self.λ_0
		while self.λ < self.λ_final:		
		    i=0
		    j=0
		    self.diagnolizePotentialArray(self.λ)		    
		    self.H = self.T + self.V		    
		    i=0
		    j=0
		    while i<self.nEigVectors:
		        j=0
		        while j<self.nEigVectors:
		            dot1=np.dot(self.H, self.u[j])
		            dot2=np.dot(self.u[i], dot1)
		            self.Hii[i,j] = dot2
		            j=j+1
		        i=i+1
		    val, vec = np.linalg.eig(self.Hii)
		    self.eigvectorFit.append(np.amin(val))
		    self.λ += self.λ_increment
		    
			


	def sweepThroughLambdas(self):
		#================== double check this ===================================
		self.λ_list, self.eigVectors, self.eigValues = [], [], [] # the "s" means list; eigVects = list of eig vectors
		# sweeping
		self.λ = self.λ_0
		while self.λ <= self.λ_final:
			# 
			self.calculateEigValueAndVector() # 
			self.λ += self.λ_increment
		#========================================================================
	def generateEmpty_T_V_arrays(self):
		# empty arrays to be filled
	    T = np.zeros((self.N-2)**2).reshape(self.N-2, self.N-2) #N-2 by N-2 array of zeros
	    V = np.zeros((self.N-2)**2).reshape(self.N-2, self.N-2) #potential energy array
	    return T, V
	def diagnolizeKineticArray(self):
	    # void function diagnolizes T
	    for i in range(self.N-2): #create the finite differenc array
	        for j in range(self.N-2):
	            if i==j:
	                self.T[i, j]= -2
	            elif abs(i-j)==1:
	                self.T[i, j]=1
	            else:
	                self.T[i, j]=0
	def diagnolizePotentialArray(self, λ):
		# fixes for some value of λ
		for i in range(self.N-2):
			for j in range(self.N-2):
				if i==j:
					self.V[i, j] = v_pot(self.x[i], λ) # check index, x[i+1]
				else:
					self.V[i, j] = 0
	def diagnolizePotentialArrayForEC(self, λ):
		# does what name says for V
		for i in range(self.N-2):
			for j in range(self.N-2):
				if i==j:
					self.V[i, j] = v_pot(self.x[i], λ) # check index, x[i+1]
				else:
					self.V[i, j] = 0
	def GS_dool(self):
		'''doolings gram-schmitt'''
		self.v = self.v.T.reshape((21, 98)) # fixed the shape
		self.u = self.v # u is basis
		i = 1
		print("u = \n", self.u.T.reshape((21, 98)))
		while i < self.nEigVectors:
			j = 0
			while j < i:
				self.u[i]=self.u[i]-(np.dot(self.u[j],self.v[i])/np.dot(self.u[j], self.u[j]))*self.u[j]
				j += 1
			i += 1
		i=0
		while i<self.nEigVectors:
			# making them unit vectors
			self.u[i] = self.u[i]/np.sqrt(np.dot(self.u[i], self.u[i]))
			i += 1
	def calculateEigValueAndVector(self):
		''' does calculations and records them'''
		self.diagnolizePotentialArray(self.λ)
		self.H = self.T + self.V # hamiltonian

		# linalg.eig does not return in order
		val, vec = np.linalg.eig(self.H) # vec[z[0]]
		z = np.argsort(val)[:self.nEigState] # <==================================double check this<============================

		self.λ_list.append(self.λ)
		self.eigValues.append(val[z][0])
		self.eigVectors.append(vec[:, z]) #<<<<<<<<<<<<<<<< LOOK HERE, THIS WHERE YOU LEFT OFFF ON !!!!!!!
		







