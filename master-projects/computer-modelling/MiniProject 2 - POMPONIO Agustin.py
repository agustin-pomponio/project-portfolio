# ==== MiniProject 2 ======================================
# ==== Computer Modelling =================================
# ==== Politechnika Krakowska =============================
# ==== Author: Agustin Pomponio ===========================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ===========================================================================
# Case 1 : isothermal conditions
# ===========================================================================

# The kinetics of degradation of compound X is a first-order kinetics, with 
# constant k as a function of temperature. 

# k has two limit values of temperature that correspond to k=0. By setting k 
# equal to zero we can find those two values, which are Tmin=12°C and Tmax=30°C. 
# For temperatures outside this range, k adopts negative or null values
# (with no physical sense)

# The functionality of k with temperature follows a quadratic function with a 
# maximum. In the case of an isothermal reaction, the fastest rate of 
# biodegradation of X will occur at the temperature in which k is maximum 

T_range = np.arange(13,30,0.5) # creates an array with a step of 0.5°C between 13 and 30°C

a = 21
b = 100
c = 0.81
k = -((T_range-a)**2)/b + c

plt.figure()
plt.plot(T_range,k)
plt.title("Fig. 1 - Functionality of k with temperature")
plt.xlabel("Temperature, °C")
plt.ylabel("Reaction constante (k), 1/day")
plt.grid()

# Finds the temperature that corresponds to the maximum value of k
T_opt = T_range[np.argmax(k)]

def model_constant(y,t):
    X = y[0]
    k = -((T-a)**2)/b + c
    dXdt = -k*X     # variation of toxic substance concentration
    return [dXdt]

t_range = np.linspace(0,15,300)     # timespan (in days)
X0 = 1                              # Initial concentration (in mg/L)

T = T_opt
ic = [X0]
results1 = odeint(model_constant, ic, t_range)
X1 = results1[:,0]

k_max = -((T_opt-a)**2)/b + c

# Looks for the moment in which the concentration of X is 0.5% of the initial
# concentration. At that moment, it is considered that the toxine has been
# biodegraded
i = 0
limit_biodeg = 0.005
# goes through the array until a concentration lower than the limit is found
while X1[i] > limit_biodeg:        
    i += 1
    
t_biodeg = t_range[i]       # time to achieve biodegradation

print("Biodegradation in isothermal case occurs most rapidly at {:.1f}°C".format(T_opt))
print("At that temperature, the maximum biodegradation constant is observed: {:.2f} 1/day".format(k_max))
print('(See Figure 1)')
print()

# ============================================================================
# Case 2 : daily temperature fluctuations 
# ============================================================================

# Temperature fluctuates daily between 13 and 30°C. As it was shown, variations
# in temperature have an impact in the value of k. 
# Behaviour of temperature throughout the different times of the day is done 
# using a trigonometrical function, so as to picture the periodical ups-and-downs
# that we can find during the day.
# The moment of maximum temperature are found in the afternoon, close to midday,
# while the coldest ones are at around 3-5 a.m., before or during sunrise.

# Modelling of temperature was done following an example provided by KhanAcademy
# Link: https://youtu.be/RX0DY9eRp8g?si=XdjqvBxDBVcgQ6FX

def temperature(t):
    T = A * np.sin(2*np.pi/B * (t-C)) + D
    return T

# Parameters of the sinoidal function:
 # Amplitude(A): absolut value of the difference between the average and the limits
 # Period(B): time it takes for the sinusoidal function to complete a "cycle"
 # Phase shift(C): 
 # Vertical shift(D): the average value between the limits (13 and 30°C) 

A = 8.5
B = 1      # It takes one day to complete a cycle of temperature
C = 9/24   # The average value is reached at 9h (1 hour = 1/24 days)
# With this phase shift the minimum is reached at 3h and the minimum at 15h
D = 21.5

def model_fluctuations(y,t):
    X = y[0]
    T = temperature(t)
    k = -((T-a)**2)/b + c
    
    dXdt = -k*X

    return [dXdt]

results2 = odeint(model_fluctuations, ic, t_range)
X2 = results2[:,0]

i=0
while X2[i] > limit_biodeg:
    i += 1
    
t_biodeg_fluc = t_range[i]

print("When temperature is constant and optimal, degradation of X takes {:.1f} days".format(t_biodeg))
print("When temperature fluctuates daily between 13 and 30°C, X degrades in {:.1f} days".format(t_biodeg_fluc))
print('(See Figures 3 and 4)')

plt.figure()
plt.plot(t_range, temperature(t_range),'b-')
plt.xlim(0,10)
plt.title("Fig.2 - Temperature modelling")
plt.xlabel("time, days")
plt.ylabel("Temperature, °C")
major_ticks = np.linspace(0,10,11)
plt.xticks(major_ticks)
plt.grid()

plt.figure()
plt.plot(t_range, X1, 'g-', label='Isothermal (Topt={:.1f}°C)'.format(T_opt))
plt.title("Fig.3 - Concentration of toxine")
plt.xlabel("time, days")
plt.ylabel("Concentraton, mg/L")
plt.plot(t_range, X2, 'b-', label="Temperature fluctuations")
plt.plot(t_range, np.linspace(limit_biodeg,limit_biodeg, 300), 'r--', label='Biodegradation limit')
plt.legend(loc='upper right')
plt.grid()
plt.ylim(0,1)
plt.xlim(0,15)

plt.figure()
plt.plot(t_range, X1, 'g-', label='Isothermal (Topt={:.1f}°C)'.format(T_opt))
plt.title("Fig.4 - Concentration of toxine (zoomed)")
plt.xlabel("time, days")
plt.ylabel("Concentraton, mg/L")
plt.plot(t_range, X2, 'b-', label="Temperature fluctuations")
plt.plot(t_range, np.linspace(limit_biodeg,limit_biodeg, 300), 'r--', label='Biodegradation limit')
plt.legend(loc='upper right')
plt.grid()
plt.ylim(0,0.1)
plt.xlim(0,15)

# From the plots it can be seen that the isothermal case at the optimum T will
# surely be faster than when temperature fluctuates daily. While the isothermal
# biodegradation curve is smooth, the curve in the fluctuation case experiences
# multiple "oscilations" as a result of the temperature variation. When the T
# is in the surroundings of the optimum (21°C), higher biodegradation rates
# are observed. On the other hand, when temperature moves away of the optimum, 
# therefore approaching the limits (Tmin=13°C and Tmax=30°C), either cooling
# or heating, a lower degradation rate is observed. All these effects make 
# the biodegradation time around 2x the time needed in the ideal isothermal
# and optimal case.