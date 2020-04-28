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


#Total population, N.
P = 11000000
# Initial number of infected and recovered individuals, I0 and R0.
I0 = 1/P #beggining date march 19th
R0 =  0/P
# Everyone else, S0, is susceptible to infection initially.
S0 = (P - I0)/P
beta, gamma = 0.6, 1/12
# A grid of time points (in days)
t = np.linspace(0, 27, 27)

# The SIR model differential equations.
def deriv(y, t, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I
    dIdt = beta * S * I  - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

### Initial conditions vector
y0 = S0, I0, R0
### Integrate the SIR equations over the time grid, t.
### do it over different parameter values???
Rdata = np.matrix([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,2,2,2,3,3,3]) # Listing real life deaths
Rdata = Rdata/P # scale down to be proportion
Idata = np.matrix([1,2,2,2,6,7,8,8,8,8,8,15,15,15,16,16,17,18,19,24,25,30,31,31,33,40,41]) # real life infections
Idata = Idata/P #scale down to be proportion
(dummy, dataLength) = Rdata.shape
S = np.zeros(dataLength)# initalize blank datasets
I = np.zeros(dataLength)
R = np.zeros(dataLength)
S[0] = S0 # all except one are susceptile
I[0] = I0 # based on one case
R[0] = R0 # all others fit in Recovered category
# two initial guesses for what beta and gamma might be

beta, gamma = 0.6, 1/12 # two initial guesses for what beta and gamma might be

tau = 0 # because first case is when time=0
fm = 0 # initialize the lsr variable
T = [100]*1001
alpha = 0.7
for j in range(0, 1000):
    T[j+1] = T[j] *(alpha**(j+1))
plt.figure()  # open the figure
fig,ax = plt.subplots()
iteration = np.linspace(0,1001, 1001) # iteration vector
ax.plot(iteration[0:15], T[0:15]) # plot 
ax.set_xlabel('Iteration') # set label, legend, title
ax.set_ylabel('Temperature')
plt.title('Temperature v Iteration')
plt.savefig('Temp Graph.png', bbox_inches ='tight') # save the graph


consec = 0 # variable to keep track of consecutive failutes
bestBeta = [] # empty vector to add on the best betas
bestGamma = [] # empty vector to add on the best gammas
bestFit = [] # empty vector to add on the lowest LSE


    


for j in range(1000):
    if T[j] > 0:
        if consec < 700: # going thru 10000 trials
            ret = odeint(deriv, y0, t, args=(beta, gamma)) # evaluate sir using the inital, or newest best guess
            S, I, R = ret.T # store S,I, R
            for i in range(dataLength -1): #calculate LSE for all of the data points, excluding the first one
                fm += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the fitted model
            deltaBeta = .05 # step size
            deltaGamma = .001 # step size
            betaOrGamma = random.random() # number btwn 0 or 1 to see which step you will do
            fhat = 0 # initalizing lsr for new beta or gamma value
            if betaOrGamma< 0.25:
                betaNew = beta + deltaBeta # add to beta
                ret = odeint(deriv, y0, t, args=(betaNew, gamma))# evaluate sir using newest beta value
                S, I, R = ret.T # store S,I,R
                for i in range(dataLength -1): # calculate LSE for all of the data points, excluding the first one
                    fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
                    if fhat < fm:
                        A = 1 # ACCEPT
                        consec = 0
                    else:
                        A = math.exp((fm-fhat)/T[j]) # maybe accept
                    acceptanceProbability = random.random() # randomly chose acceptance probability
                    if acceptanceProbability< A:
                        fm = fhat # then we have a new lsr
                        beta = betaNew # and a new beta Value
                        bestBeta.append(beta)
                        bestGamma.append(gamma)
                        bestFit.append(fm)
                    else:
                        consec +=1
            elif 0.5>betaOrGamma> 0.25 and beta > 0:  #added to remove all the bad negative estimates
                betaNew = beta - deltaBeta # then subtract from beta
                ret = odeint(deriv, y0, t, args=(betaNew, gamma))# add commas for more parameters
                S, I, R = ret.T
                for i in range(dataLength - 1):
                    fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
                    if fhat < fm:
                        A = 1 # then ACCEPT
                    else:
                        A = math.exp((fm-fhat)/T[j]) # THEN maybe accept
                    acceptanceProbability = random.random() # randomly chose acceptance probability
                    if acceptanceProbability< A:
                        fm = fhat # adopt new lsr
                        beta = betaNew # and a new beta value
                        consec = 0
                        bestBeta.append(beta)
                        bestGamma.append(gamma)
                        bestFit.append(fm)
                    else:
                        consec +=1
            elif 0.75>betaOrGamma> 0.5:
                gammaNew = gamma + deltaGamma  # then add to gamma
                ret = odeint(deriv, y0, t, args=(beta, gammaNew))# add commas for more parameters
                S, I, R = ret.T
                for i in range(dataLength - 1):
                    fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
                    if fhat < fm:
                        A = 1 # accept
                        consec = 0
                    else:
                        A = math.exp((fm-fhat)/T[j]) # calculation of A is done wrong
                    acceptanceProbability = random.random() # randomly chose acceptance probability
                    if acceptanceProbability< A:
                        fm = fhat
                        gamma = gammaNew
                        bestBeta.append(beta)
                        bestGamma.append(gamma)
                        bestFit.append(fm)
                    else:
                        consec +=1
            elif 1>betaOrGamma> 0.75 and gamma > 0: #added to remove all the bad negative estimates
                gammaNew = gamma - deltaGamma # then subtract from gamma
                ret = odeint(deriv, y0, t, args=(beta, gammaNew))# add commas for more parameters
                S, I, R = ret.T
                for i in range(dataLength - 1):
                    fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
                    if fhat < fm:
                        A = 1
                        consec = 0
                    else:
                        A = math.exp((fm-fhat)/T[j]) # maybe accept, calculation of A is done wrong
                    acceptanceProbability = random.random() # randomly chose acceptance probability
                    if acceptanceProbability< A:
                        fm = fhat # adopt new lsr and new gamma value
                        gamma = gammaNew
                        bestBeta.append(beta)
                        bestGamma.append(gamma)
                        bestFit.append(fm)
                    else:
                        consec +=1


# Getting the best of the best, which doesnt depend on the end point

index = bestFit.index(min(bestFit)) # find day where number of infections is highest
bestBESTbeta = bestBeta[index] # finding best of bestBeta
bestBESTgamma = bestGamma[index] # finding best of best gamma
print(bestBESTbeta, bestBESTgamma, min(bestFit)) # reporting the best best parameters
# and their corresponding error



#Evaluate SIR using best parameters, and dead using the estimate of death rate
ret = odeint(deriv, y0, t, args=(bestBESTbeta, bestBESTgamma))
S, I, R = ret.T # solve for S, I, R, using best best parameters
Dead = 3/41 * R # using average death rate from the observed data



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
plt.title('First 27 days')

plt.savefig('SimAnnealing.png', bbox_inches ='tight')


## PLOT the data for the observed real time data to get an understanding of the trends before we are
## certain that this is a good fit
plt.figure()  # open the figure
fig,ax = plt.subplots()
ax.plot(realTime,Idata, 'o')
plt.savefig('realCurve.png', bbox_inches ='tight')
Time =  np.linspace(0, 160, 160)
ret = odeint(deriv, y0, Time, args=(bestBESTbeta, bestBESTgamma))# add commas for more parameters
S, I, R = ret.T
Dead = 3/41 * R # using average death rate from the observed data

# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.figure()  # open the figure
fig,ax = plt.subplots()
realTime = np.linspace(0,dataLength, dataLength) # time vector
Idata = Idata.transpose()
ax.plot(Time, S, 'b', alpha=0.5, lw=2, label='Susceptible') # plot all the different lines
ax.plot(Time, I, 'm', alpha=0.5, lw=2, label='Infected')
ax.plot(Time, Dead, 'r', alpha=0.5, lw=2, label='Dead')
ax.plot(Time, R, 'g', alpha=0.5, lw=2, label='Recovered with immunity or Dead')
ax.set_xlabel('Time (days)') # set label, legend, title
ax.set_ylabel('Population')
legend = ax.legend(loc = 'upper left')
plt.title('Estimated SIR Model Using Best Guesses')
plt.savefig('SimAnnealing160days.png', bbox_inches ='tight') # save the graph

Ilist = I.tolist() # reformatting to eliminate bugs
maxIday = Ilist.index(max(I)) # find day where number of infections is highest
maxDeath = max(Dead)*P
maxCases = max(I)*P

# printing relevant data parameters for the user
print('On day %d, the number of infections is the highest. ' %maxIday)
print('The maximum number of deaths is %d. ' %maxDeath)
print('The maximum number of cases is %d. ' %maxCases)



