import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
from scipy import integrate
import matplotlib.pyplot as plt
import random

'''
Question 2:
Estimating good estimates for the parameters in an SEIR using least squares regression:

'''
P = 11000000 # current population of haiti
Rdata = np.matrix([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,2,2,2,3,3,3]) # Listing real life deaths
Rdata = Rdata/P
Idata = np.matrix([1,2,2,2,6,7,8,8,8,8,8,15,15,15,16,16,17,18,19,24,25,30,31,31,33,40,41]) # real life infections
Idata = Idata/P
(dummy, dataLength) = Rdata.shape
S = np.zeros(dataLength)# initalize blank datasets
I = np.zeros(dataLength)
R = np.zeros(dataLength)
S[0] = (P-1)/P
I[0] = 1/P
R[0] = 1- S.item(0) -I.item(0)

beta = 0.6 # two initial guesses for what beta and gamma might be
gamma = 1/12 # two initial guesses for what beta and gamma might be
tau = 0 # because first case is when time=0
fm = 0
alpha = .5 # cooling schedule, arbitratily chosen
T = 100 # intial temperature,arbitratily chosen
for j in range(10000):
    for i in range(dataLength-1):
        dS = -beta*S[i]*I[i]
        dI = beta*S[i]*I[i] - gamma*I[i]
        dR = gamma*I[i]
        I[i+1] = I[i] + dI
        fm += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the fitted model
    deltaBeta = .05 # step size
    deltaGamma = .001 # step size
    betaOrGamma = random.random()
    fhat = 0
    if betaOrGamma< 0.25:
        betaNew = beta + deltaBeta
        for i in range(dataLength - 1):
            dS = -betaNew*S[i]*I[i]
            dI = betaNew*S[i]*I[i] - gamma*I[i]
            dR = gamma*I[i]
            I[i+1] = I[i] + dI
            fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
            if fhat < fm:
                A = 1
            else:
                A = math.exp((fm-fhat)/T) # calculation of A is done wrong
            acceptanceProbability = random.random() # randomly chose acceptance probability
            if acceptanceProbability< A:
                fm = fhat
                beta = betaNew
            else:
              pass
    elif 0.5>betaOrGamma> 0.25:
        betaNew = beta - deltaBeta
        for i in range(dataLength - 1):
            dS = -betaNew*S[i]*I[i]
            dI = betaNew*S[i]*I[i] - gamma*I[i]
            dR = gamma*I[i]
            I[i+1] = I[i] + dI
            fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
            if fhat < fm:
                A = 1
            else:
                A = math.exp((fm-fhat)/T) # calculation of A is done wrong
            acceptanceProbability = random.random() # randomly chose acceptance probability
            if acceptanceProbability< A:
                fm = fhat
                beta = betaNew
            else:
              pass
    elif 0.75>betaOrGamma> 0.5:
        gammaNew = gamma + deltaGamma
        for i in range(dataLength - 1):
            dS = -beta*S[i]*I[i]
            dI = beta*S[i]*I[i] - gammaNew*I[i]
            dR = gammaNew*I[i]
            I[i+1] = I[i] + dI
            fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
            if fhat < fm:
                A = 1
            else:
                A = math.exp((fm-fhat)/T) # calculation of A is done wrong
            acceptanceProbability = random.random() # randomly chose acceptance probability
            if acceptanceProbability< A:
                fm = fhat
                gamma = gammaNew
            else:
              pass
    elif 1>betaOrGamma> 0.75:
        gammaNew = gamma - deltaGamma
        for i in range(dataLength - 1):
            dS = -beta*S[i]*I[i]
            dI = beta*S[i]*I[i] - gammaNew*I[i]
            dR = gammaNew*I[i]
            I[i+1] = I[i] + dI
            fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
            if fhat < fm:
                A = 1
            else:
                A = math.exp((fm-fhat)/T) # calculation of A is done wrong
            acceptanceProbability = random.random() # randomly chose acceptance probability
            if acceptanceProbability< A:
                fm = fhat
                gamma = gammaNew
            else:
              pass
print(beta,gamma)
t =  np.linspace(0, 27, 27)
S = np.zeros(27)# initalize blank datasets
I = np.zeros(27)
R = np.zeros(27)
S[0] = (P-1)/P
I[0] = 1/P
R[0] = 1- S.item(0) -I.item(0)
for i in range(26):
    dS = - beta*S[i]*I[i]
    dI = beta*S[i]*I[i] - gamma*I[i]
    dR = gamma*I[i]
    S[i+1] = S[i] + dS
    I[i+1] = I[i] + dI
    R[i+1] = R[i] + dR
    Dead = 3/41 * R
# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.figure()  # open the figure
fig,ax = plt.subplots()
realTime = np.linspace(0,dataLength, dataLength)
Idata = Idata.transpose()
ax.plot(t, I, 'm', alpha=0.5, lw=2, label='Model Data')
ax.plot(realTime,Idata, 'o')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Population')
legend = ax.legend(loc = 'upper left')
plt.savefig('SimAnnealing.png', bbox_inches ='tight')

plt.figure()  # open the figure
fig,ax = plt.subplots()
ax.plot(realTime,Idata, 'o')
plt.savefig('realCurve.png', bbox_inches ='tight')
beta = .5
gamma  = 1/30
t =  np.linspace(0, 160, 160)
S = np.zeros(160)# initalize blank datasets
I = np.zeros(160)
R = np.zeros(160)
S[0] = (P-1)/P
I[0] = 1/P
R[0] = 1- S.item(0) -I.item(0)
for i in range(159):
    dS = - beta*S[i]*I[i]
    dI = beta*S[i]*I[i] - gamma*I[i]
    dR = gamma*I[i]
    S[i+1] = S[i] + dS
    I[i+1] = I[i] + dI
    R[i+1] = R[i] + dR
    Dead = 3/41 * R
# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.figure()  # open the figure
fig,ax = plt.subplots()
realTime = np.linspace(0,dataLength, dataLength)
Idata = Idata.transpose()
deathRate = .034
Dead = deathRate*I
ax.plot(t, S, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I, 'm', alpha=0.5, lw=2, label='Infected')
ax.plot(t, Dead, 'r', alpha=0.5, lw=2, label='Dead')
ax.plot(t, R, 'g', alpha=0.5, lw=2, label='Recovered with immunity or Dead')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Population')
legend = ax.legend(loc = 'upper left')
plt.savefig('SimAnnealing160days.png', bbox_inches ='tight')
Ilist = I.tolist()
maxIday = Ilist.index(max(I))
maxDeath = max(Dead) *P
maxCases = max(I)*P
