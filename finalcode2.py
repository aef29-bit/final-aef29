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

##'''
P = 11000000 # current population of haiti
Rdata = np.matrix([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,2,2,2,3,3,3]) # Listing real life deaths
Rdata = Rdata/P # scale down to be proportion
Idata = np.matrix([1,2,2,2,6,7,8,8,8,8,8,15,15,15,16,16,17,18,19,24,25,30,31,31,33,40,41]) # real life infections
Idata = Idata/P #scale down to be proportion
(dummy, dataLength) = Rdata.shape
S = np.zeros(dataLength)# initalize blank datasets
I = np.zeros(dataLength)
R = np.zeros(dataLength)
S[0] = (P-1)/P # all except one are susceptible
I[0] = 1/P # based on one case
R[0] = 1- S.item(0) -I.item(0) # all others fit in Recovered category
# keep track of all the errors, started 10,000 trials
# beta, gamma  = 0.6, 1/12 fn = 8.6651411177292e-11
# beta, gamma  = 0.9, 1/10 fn = 9.25118160473985e-11, BAD guesses in negatives, higher error
# # beta, gamma  = 0.5, 1/10 fn = 8.923713165658255e-11 , BAD guesses in negatives, higher error
# # beta, gamma  = 2, 1/11 fn = 9.139455163446431e-11, BAD guesses in negatives, higher error
# # beta, gamma  = 1, 1/30 fn = 1.2248778563738362e-10 , BAD guesses in negatives, higher error


# with the starting guess, 0.6 and 1/12
# beta, gamma  =  1.0 0.10133333333333333, lsr = 8.799635676839043e-11
# if u make that your starting guess, error goes up
beta = 0.6 # two initial guesses for what beta and gamma might be
gamma = 1/12 # two initial guesses for what beta and gamma might be
tau = 0 # because first case is when time=0
fm = 0 # initialize the lsr variable
T = 100 # intial temperature,arbitratily chosen
consec = 0
if consec < 7000:
    for j in range(10000): # going thru 10000 trials
        for i in range(dataLength-1): # go thru all terms to get lsr
            dS = -beta*S[i]*I[i] # evaluate derivatives
            dI = beta*S[i]*I[i] - gamma*I[i]
            dR = gamma*I[i]
            I[i+1] = I[i] + dI
            fm += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the fitted model
        deltaBeta = .05 # step size
        deltaGamma = .001 # step size
        betaOrGamma = random.random() # number btwn 0 or 1 to see which step you will do
        fhat = 0 # initalizing lsr for new beta or gamma value
        if betaOrGamma< 0.25:
            betaNew = beta + deltaBeta # add to beta
            for i in range(dataLength - 1):
                dS = -betaNew*S[i]*I[i] # reevaluate derivatives using new beta value
                dI = betaNew*S[i]*I[i] - gamma*I[i]
                dR = gamma*I[i]
                I[i+1] = I[i] + dI
                fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
                if fhat < fm:
                    A = 1 # ACCEPT
                    consec = 0
                else:
                    A = math.exp((fm-fhat)/T) # maybe accept
                acceptanceProbability = random.random() # randomly chose acceptance probability
                if acceptanceProbability< A:
                    fm = fhat # then we have a new lsr
                    beta = betaNew # and a new beta Value
                    consec = 0
                else:
                    consec +=1
        elif 0.5>betaOrGamma> 0.25 and beta > 0:  #added to remove all the bad negative estimates
            betaNew = beta - deltaBeta # then subtract from beta
            for i in range(dataLength - 1):
                dS = -betaNew*S[i]*I[i] # reevaluate derivatives using new beta value
                dI = betaNew*S[i]*I[i] - gamma*I[i]
                dR = gamma*I[i]
                I[i+1] = I[i] + dI
                fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
                if fhat < fm:
                    A = 1 # then ACCEPT
                    consec = 0
                else:
                    A = math.exp((fm-fhat)/T) # THEN maybe accept
                acceptanceProbability = random.random() # randomly chose acceptance probability
                if acceptanceProbability< A:
                    fm = fhat # adopt new lsr
                    beta = betaNew # and a new beta value
                    consec = 0
                else:
                    consec +=1
        elif 0.75>betaOrGamma> 0.5:
            gammaNew = gamma + deltaGamma  # then add to gamma
            for i in range(dataLength - 1):
                dS = -beta*S[i]*I[i] # reeval derivs using new gamma value
                dI = beta*S[i]*I[i] - gammaNew*I[i]
                dR = gammaNew*I[i]
                I[i+1] = I[i] + dI
                fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
                if fhat < fm:
                    A = 1 # accept
                    consec = 0
                else:
                    A = math.exp((fm-fhat)/T) # calculation of A is done wrong
                acceptanceProbability = random.random() # randomly chose acceptance probability
                if acceptanceProbability< A:
                    fm = fhat
                    gamma = gammaNew
                    consec = 0
                else:
                    consec +=1
        elif 1>betaOrGamma> 0.75 and gamma > 0: #added to remove all the bad negative estimates
            gammaNew = gamma - deltaGamma # then subtract from gamma
            for i in range(dataLength - 1):
                dS = -beta*S[i]*I[i]
                dI = beta*S[i]*I[i] - gammaNew*I[i]
                dR = gammaNew*I[i]
                I[i+1] = I[i] + dI
                fhat += abs((I.item(i+1)-Idata.item(i+1))**2) # square of the error of the neighbors model
                if fhat < fm:
                    A = 1
                    consec = 0
                else:
                    A = math.exp((fm-fhat)/T) # maybe accept, calculation of A is done wrong
                acceptanceProbability = random.random() # randomly chose acceptance probability
                if acceptanceProbability< A:
                    fm = fhat # adopt new lsr and new gamma value
                    gamma = gammaNew
                    consec = 0
                else:
                    consec +=1
print(beta,gamma)
print(fm)
t =  np.linspace(0, 27, 27) # make time vector
S = np.zeros(27)# initalize blank datasets
I = np.zeros(27)
R = np.zeros(27)
S[0] = (P-1)/P # first item is certain by assumptions of initial state
I[0] = 1/P
R[0] = 1- S.item(0) -I.item(0)
for i in range(26): # make SIR vectors using new beta and gamma, for all except first term
    dS = - beta*S[i]*I[i]
    dI = beta*S[i]*I[i] - gamma*I[i]
    dR = gamma*I[i]
    S[i+1] = S[i] + dS
    I[i+1] = I[i] + dI
    R[i+1] = R[i] + dR
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
plt.savefig('SimAnnealing.png', bbox_inches ='tight')


## PLOT the data for the observed real time data to get an understanding of the trends before we are
## certain that this is a good fit
plt.figure()  # open the figure
fig,ax = plt.subplots()
ax.plot(realTime,Idata, 'o')
plt.savefig('realCurve.png', bbox_inches ='tight')
t =  np.linspace(0, 160, 160)
S = np.zeros(160)# initalize blank datasets
I = np.zeros(160)
R = np.zeros(160)
S[0] = (P-1)/P # assumptions about the initial state
I[0] = 1/P
R[0] = 1- S.item(0) -I.item(0)
for i in range(159): # calculate vector over 160 days
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
realTime = np.linspace(0,dataLength, dataLength) # time vector
Idata = Idata.transpose()
deathRate = 3/41 # observed in the data
Dead = deathRate*I
ax.plot(t, S, 'b', alpha=0.5, lw=2, label='Susceptible') # plot all the different lines
ax.plot(t, I, 'm', alpha=0.5, lw=2, label='Infected')
ax.plot(t, Dead, 'r', alpha=0.5, lw=2, label='Dead')
ax.plot(t, R, 'g', alpha=0.5, lw=2, label='Recovered with immunity or Dead')
ax.set_xlabel('Time (days)') # set label, legend, title
ax.set_ylabel('Population')
legend = ax.legend(loc = 'upper left')
plt.savefig('SimAnnealing160days.png', bbox_inches ='tight') # save the graph
Ilist = I.tolist() # reformatting to eliminate bugs
maxIday = Ilist.index(max(I)) # find day where number of infections is highest
maxDeath = max(Dead)*P
maxCases = max(I)*P
print('On day %d, the number of infections is the highest. ' %maxIday)
print('The maximum number of deaths is %d. ' %maxDeath)
print('The maximum number of cases is %d. ' %maxCases)



