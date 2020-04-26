import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
from scipy import integrate
import matplotlib.pyplot as plt
'''
Question 1: Seing How Changing Different Parameters Affects the Day Gap, The height of the I curve
and the amount of days untill the peak in case
'''

'''
Part a:
Using good estimates for the parameters provided by a variety of different sources:

'''

# Total population, N.
N = 11000000
# Initial number of infected and recovered individuals, I0 and R0.
I0 = 2 #beggining date march 19th
E0 = 0
R0 =  0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - E0 - R0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days). #CHANGE WORDS HERE
beta, sigma, gamma = 0.78735, 1/7, 0.154 # sigma = 1/ average duration
# A grid of time points (in days)
t = np.linspace(0, 160, 160)

# The SIR model differential equations.
def deriv(y, t, N, beta, sigma, gamma):
    S, E, I, R = y
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - sigma*E
    dIdt = sigma*E- gamma * I
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt

# Initial conditions vector
y0 = S0, E0, I0, R0
# Integrate the SIR equations over the time grid, t.
# do it over different parameter values???
ret = odeint(deriv, y0, t, args=(N, beta, sigma, gamma))# add commas for more parameters
S, E, I, R = ret.T
Dead1 = .034*R
DeadinEnd1 = Dead1[-1]
AliveRecovered1 = (1-.034)*R

RecoveredinEnd1 = AliveRecovered1[-1]
print(' The model estimates that %d are dead by the end of the period of 160 days.' %DeadinEnd1)
print(' The model estimates that %d are recovered by the end of the period of 160 days.' %RecoveredinEnd1)

# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.figure()  # open the figure
fig,ax = plt.subplots() 
ax.plot(t, S/1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, E/1000, 'y', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, I/1000, 'm', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/1000, 'g', alpha=0.5, lw=2, label='Recovered with immunity or Dead')
ax.plot(t, Dead1/1000, linestyle='--', alpha=0.5, lw=2, label='Dead')

ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
legend = ax.legend(loc = 'upper left')
plt.title('SEIR 1st Graph')
plt.savefig('SEIR1.png', bbox_inches ='tight')
# the day gap between peak exposure and peak number of infected is relevant
Elist = E.tolist()
Ilist = I.tolist()
maxEday1 = Elist.index(max(E))
maxIday1 = Ilist.index(max(I))
dayGap1 = maxIday1 - maxEday1
maxDeath1 = max(Dead1)
maxCases1 = max(I)

'''
Part b:
Lowering and increasing beta based on
https://www.businessinsider.com/coronavirus-contagious-r-naught-average-patient-spread-2020-3

'''
beta, sigma, gamma = 0.5, 1/7, 0.154 # beta is estimated number of contacts per day
# Initial conditions vector
y0 = S0, E0, I0, R0
# Integrate the SIR equations over the time grid, t.
# do it over different parameter values???
ret = odeint(deriv, y0, t, args=(N, beta, sigma, gamma))# add commas for more parameters
S, E, I, R = ret.T
Dead2 = .034*R
DeadinEnd2 = Dead2[-1]
AliveRecovered2 = (1-.034)*R

RecoveredinEnd2 = AliveRecovered2[-1]
print(' The model estimates that %d are dead by the end of the period of 160 days.' %DeadinEnd2)
print(' The model estimates that %d are recovered by the end of the period of 160 days.' %RecoveredinEnd2)

# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.figure()  # open the figure
fig,ax = plt.subplots() 
ax.plot(t, S/1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, E/1000, 'y', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, I/1000, 'm', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/1000, 'g', alpha=0.5, lw=2, label='Recovered with immunity or Dead')
ax.plot(t, Dead2/1000, linestyle='--', alpha=0.5, lw=2, label='Dead')

ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
legend = ax.legend(loc = 'upper left')
plt.title('SEIR 2nd Graph')

plt.savefig('SEIR2.png', bbox_inches ='tight')
# the day gap between peak exposure and peak number of infected is relevant
Elist = E.tolist()
Ilist = I.tolist()
maxEday2 = Elist.index(max(E))
maxIday2 = Ilist.index(max(I))
dayGap2 = maxIday2 - maxEday2
maxDeath2 = max(Dead2)
maxCases2 = max(I)
beta, sigma, gamma = 2, 1/7, 0.154
# Initial conditions vector
y0 = S0, E0, I0, R0
# Integrate the SIR equations over the time grid, t.
# do it over different parameter values???
ret = odeint(deriv, y0, t, args=(N, beta, sigma, gamma))# add commas for more parameters
S, E, I, R = ret.T
Dead3 = .034*R
DeadinEnd3 = Dead3[-1]
AliveRecovered3 = (1-.034)*R

RecoveredinEnd3 = AliveRecovered3[-1]
print(' The model estimates that %d are dead by the end of the period of 160 days.' %DeadinEnd3)
print(' The model estimates that %d are recovered by the end of the period of 160 days.' %RecoveredinEnd3)

# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.figure()  # open the figure
fig,ax = plt.subplots() 
ax.plot(t, S/1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, E/1000, 'y', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, I/1000, 'm', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/1000, 'g', alpha=0.5, lw=2, label='Recovered with immunity or Dead')
ax.plot(t, Dead3/1000, linestyle='--', alpha=0.5, lw=2, label='Dead')

ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
legend = ax.legend(loc = 'upper left')
plt.title('SEIR 3rd Graph')

plt.savefig('SEIR3.png', bbox_inches ='tight')
# the day gap between peak exposure and peak number of infected is relevant
Elist = E.tolist()
Ilist = I.tolist()
maxEday3 = Elist.index(max(E))
maxIday3 = Ilist.index(max(I))
dayGap3 = maxIday3 - maxEday3
maxDeath3 = max(Dead3)
maxCases3 = max(I)



'''
Part c:
increasing sigma based on
12 days: Average duration of survivors
https://www.telegraph.co.uk/news/2020/03/12/coronavirus-kills-average-185-days/
18.5 days: average time untill death
https://www.telegraph.co.uk/news/2020/03/12/coronavirus-kills-average-185-days/
'''
beta, sigma, gamma = 0.78735, 1/12, 0.154 # sigma = 1/ average duration
# Initial conditions vector
y0 = S0, E0, I0, R0
# Integrate the SIR equations over the time grid, t.
# do it over different parameter values???
ret = odeint(deriv, y0, t, args=(N, beta, sigma, gamma))# add commas for more parameters
S, E, I, R = ret.T
Dead4 = .034*R
DeadinEnd4 = Dead4[-1]
AliveRecovered4 = (1-.034)*R

RecoveredinEnd4 = AliveRecovered4[-1]
print(' The model estimates that %d are dead by the end of the period of 160 days.' %DeadinEnd4)
print(' The model estimates that %d are recovered by the end of the period of 160 days.' %RecoveredinEnd4)

# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.figure()  # open the figure
fig,ax = plt.subplots() 
ax.plot(t, S/1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, E/1000, 'y', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, I/1000, 'm', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/1000, 'g', alpha=0.5, lw=2, label='Recovered with immunity or Dead')
ax.plot(t, Dead4/1000, linestyle='--', alpha=0.5, lw=2, label='Dead')

ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
legend = ax.legend(loc = 'upper left')
plt.title('SEIR 4th Graph')

plt.savefig('SEIR4.png', bbox_inches ='tight')
# the day gap between peak exposure and peak number of infected is relevant
Elist = E.tolist()
Ilist = I.tolist()
maxEday4 = Elist.index(max(E))
maxIday4 = Ilist.index(max(I))
dayGap4 = maxIday4 - maxEday4
maxDeath4 = max(Dead4)
maxCases4 = max(I)
beta, sigma, gamma = 0.78735, 1/18.5, 0.154 # sigma = 1/ average duration

# Initial conditions vector
y0 = S0, E0, I0, R0
# Integrate the SIR equations over the time grid, t.
# do it over different parameter values???
ret = odeint(deriv, y0, t, args=(N, beta, sigma, gamma))# add commas for more parameters
S, E, I, R = ret.T
Dead5 = .034*R
DeadinEnd5 = Dead5[-1]
AliveRecovered5 = (1-.034)*R

RecoveredinEnd5 = AliveRecovered5[-1]
print(' The model estimates that %d are dead by the end of the period of 160 days.' %DeadinEnd5)
print(' The model estimates that %d are recovered by the end of the period of 160 days.' %RecoveredinEnd5)

# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.figure()  # open the figure
fig,ax = plt.subplots() 
ax.plot(t, S/1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, E/1000, 'y', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, I/1000, 'm', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/1000, 'g', alpha=0.5, lw=2, label='Recovered with immunity or Dead')
ax.plot(t, Dead5/1000, linestyle='--', alpha=0.5, lw=2, label='Dead')

ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
legend = ax.legend(loc = 'upper left')
plt.title('SEIR 5th Graph')
plt.savefig('SEIR5.png', bbox_inches ='tight')
# the day gap between peak exposure and peak number of infected is relevant
Elist = E.tolist()
Ilist = I.tolist()
maxEday5 = Elist.index(max(E))
maxIday5 = Ilist.index(max(I))
dayGap5 = maxIday5 - maxEday5
maxDeath5 = max(Dead5)
maxCases5 = max(I)
'''
Part d:
 γ is the average rate of recovery or death in infected populations
Lowering and increasing gamma based on
http://jtd.amegroups.com/article/view/36385/html
'''
beta, sigma, gamma = 0.78735, 1/12,0.0721
# Initial conditions vector
y0 = S0, E0, I0, R0
# Integrate the SIR equations over the time grid, t.
# do it over different parameter values???
ret = odeint(deriv, y0, t, args=(N, beta, sigma, gamma))# add commas for more parameters
S, E, I, R = ret.T
Dead6 = .034*R
DeadinEnd6 = Dead6[-1]
AliveRecovered6 = (1-.034)*R

RecoveredinEnd6 = AliveRecovered6[-1]
print(' The model estimates that %d are dead by the end of the period of 160 days.' %DeadinEnd6)
print(' The model estimates that %d are recovered by the end of the period of 160 days.' %RecoveredinEnd6)

# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.figure()  # open the figure
fig,ax = plt.subplots() 
ax.plot(t, S/1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, E/1000, 'y', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, I/1000, 'm', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/1000, 'g', alpha=0.5, lw=2, label='Recovered with immunity or Dead')
ax.plot(t, Dead6/1000, linestyle='--', alpha=0.5, lw=2, label='Dead')

ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
legend = ax.legend(loc = 'upper left')
plt.title('SEIR 6th Graph')
plt.savefig('SEIR6.png', bbox_inches ='tight')
# the day gap between peak exposure and peak number of infected is relevant
Elist = E.tolist()
Ilist = I.tolist()
maxEday6 = Elist.index(max(E))
maxIday6 = Ilist.index(max(I))
dayGap6 = maxIday6 - maxEday6
maxDeath6 = max(Dead6)
maxCases6 = max(I)



beta, sigma, gamma = 0.78735, 1/12, 0.238

# Initial conditions vector
y0 = S0, E0, I0, R0
# Integrate the SIR equations over the time grid, t.
# do it over different parameter values???
ret = odeint(deriv, y0, t, args=(N, beta, sigma, gamma))# add commas for more parameters
S, E, I, R = ret.T
Dead7 = .034*R
DeadinEnd7 = Dead5[-1]
AliveRecovered7 = (1-.034)*R

RecoveredinEnd7 = AliveRecovered7[-1]
print(' The model estimates that %d are dead by the end of the period of 160 days.' %DeadinEnd7)
print(' The model estimates that %d are recovered by the end of the period of 160 days.' %RecoveredinEnd7)

# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.figure()  # open the figure
fig,ax = plt.subplots() 
ax.plot(t, S/1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, E/1000, 'y', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, I/1000, 'm', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/1000, 'g', alpha=0.5, lw=2, label='Recovered with immunity or Dead')
ax.plot(t, Dead7/1000, linestyle='--', alpha=0.5, lw=2, label='Dead')

ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
legend = ax.legend(loc = 'upper left')
plt.title('SEIR 7th Graph')
plt.savefig('SEIR7.png', bbox_inches ='tight')
# the day gap between peak exposure and peak number of infected is relevant
Elist = E.tolist()
Ilist = I.tolist()
maxEday7 = Elist.index(max(E))
maxIday7 = Ilist.index(max(I))
dayGap7 = maxIday7 - maxEday7
maxDeath7 = max(Dead7)
maxCases7 = max(I)



