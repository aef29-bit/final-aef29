import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
from scipy import integrate
import matplotlib.pyplot as plt


'''
Question 2:
Estimating good estimates for the parameters in an SEIR using least squares regression:

'''
N = 11000000 # current population of haiti
Rdata = np.matrix([0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,2,2,2,3,3, 3]) # Listing real life deaths
Idata = np.matrix([1,2,2,2,6,7, 8,8,8,8,8,15,15,15,16,16,17,18,19,24,25,30, 31,31,33,40, 41]) # real life infections
S = np.zeros(len(Idata))# initalize blank datasets
I = np.zeros(len(Idata))
R = np.zeros(len(Idata))
S[0] = (N-Idata.item(0))/N
I[0] = Idata.item(0)/N
R[0] = 1- S.item(0) -I.item(0)
w = np.zeros(len(Idata))
x = np.zeros(len(Idata))
y = np.zeros(len(Idata))
z = np.zeros(len(Idata))
w[0], x[0], y[0], z[0] = 0,0,0,0 # intial conditions for the partial derivatives

beta = 0.5 # two initial guesses for what beta and gamma might be
gamma = 1/14 # two initial guesses for what beta and gamma might be
vCurrent = [[beta],[gamma]] #starting guess for the parameters
# Now solve the 6th order system for SI, wxyz
alpha = 0.5
gradient = np.matrix([[0],[0]])
for m in range(100):
    for i in range(len(Idata)-1):
        t = i
        S[i+1] = - beta*S.item(i)*I.item(i)
        I[i+1] = beta*S.item(i)*I.item(i) - gamma*I.item(i)
        R[i+1] = gamma*I.item(i)
        w[i+1] = -S.item(i)*I.item(i) - beta*w.item(i)*I.item(i) -beta*S.item(i)*x.item(i)
        x[i+1] = S.item(i)*I.item(i) + beta*w.item(i)*I.item(i) + beta*S.item(i)*x.item(i) - gamma*x.item(i)
        y[i+1] = -beta*y.item(i)*I.item(i) - beta*S.item(i)*z.item(i)                                     
        z[i+1]= beta*y.item(i)
        gradient[0] = gradient.item(0) + 2*abs(Idata.item(i)-I.item(i))*x.item(i)
        gradient[1] = gradient.item(1)+ 2*abs(Idata.item(i)-I.item(i))*z.item(i)
        print(gradient)
    beta += alpha*gradient.item(0)
    gamma += alpha*gradient.item(1)
    vCurrent[0] = beta
    vCurrent[1] = gamma

  

