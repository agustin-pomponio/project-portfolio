# ==== MiniProject 1 ======================================
# ==== Computer Modelling =================================
# ==== Politechnika Krakowska =============================
# ==== Author: Agustin Pomponio ===========================

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Mass of each commercial mixture [grams] - Commercial mixtures are consituted
# by the compound itself (ethanol, methanol, lactic acid) + water
mass_EtOH = 500
mass_MeOH = 500
mass_LA = 1000

    # EtOH = ethanol
    # MeOH = methanol
    # LA = lactic acid
    # H2O = water
    # ML = methyl lactate
    # EL = ethyl lactate

# Composition of each commercial mixture [%wt]
comp_EtOH = 0.96
comp_MeOH = 1
comp_LA = 0.85
# The remaining amounts (1-comp[i]) correspond to water

# Density of commercial solutions [g/mL] - Already include water content
den_EtOH = 0.789
den_MeOH = 0.792
den_LA = 1.209

# Molecular weight of each component [g/mol]
MW_EtOH = 46.068
MW_MeOH = 32.04
MW_LA = 90.08
MW_H2O = 18.015

# Mass of each component [g]
m_EtOH = mass_EtOH * comp_EtOH
m_MeOH = mass_MeOH * comp_MeOH
m_LA = mass_LA * comp_LA
m_H2O = mass_EtOH * (1-comp_EtOH) + mass_LA * (1-comp_LA)

# Individual and total volume [mL] - Water is included in the volume of each
# commercial mixture
vol_EtOH = mass_EtOH / den_EtOH
vol_MeOH = mass_MeOH / den_MeOH
vol_LA = mass_LA / den_LA
volume = vol_EtOH + vol_MeOH + vol_LA

# Initial molar concentration of each component [mol/L] (mass/(MW*volume))
# Mulitplying by 1000 is needed to convert mL to L
C_EtOH = m_EtOH * 1000 / (MW_EtOH * volume)
C_MeOH = m_MeOH * 1000 / (MW_MeOH * volume)
C_LA = m_LA * 1000 / (MW_LA * volume)
C_H2O = m_H2O * 1000 / (MW_H2O * volume)

print("Initial concentration of the reactants (in mol/L):")
print("[EtOH]o = {:.3f}".format(C_EtOH))
print("[MeOH]o = {:.3f}".format(C_MeOH))
print("[LA]o = {:.3f}".format(C_LA))
print("[H2O]o = {:.3f}".format(C_H2O))
print(" ")

# Model for resolution of exercise 1:    
def model(X):
    x1, x2 = X                      # Reaction extent of each reaction
    
    ML = x1                         # Concentration variation of each component
    EL = x2                         # during the reactions
    MeOH = C_MeOH - x1
    LA = C_LA - x1 - x2
    EtOH = C_EtOH - x2
    H2O = C_H2O + x1 + x2
    
    eq1 = ( ML * H2O ) / ( LA * MeOH ) - K1
    eq2 = ( EL * H2O ) / ( LA * EtOH ) - K2
    
    return [eq1, eq2]

K1 = 11
K2 = 25

guess = [0.01, 0.01]

x1, x2 = fsolve(model, guess)

ML = x1                         
EL = x2                         
MeOH = C_MeOH - x1
LA = C_LA - x1 - x2
EtOH = C_EtOH - x2
H2O = C_H2O + x1 + x2

print("Reaction extent for each reaction at equilibrium is (in mol/L):")
print("x1 = {:.4f}".format(x1))
print("x2 = {:.4f}".format(x2))
print()
print("Concentration of each compound at equilibrium (in mol/L):")
print("[EtOH] = {:.4f}".format(EtOH))
print("[MeOH] = {:.4f}".format(MeOH))
print("[LA] = {:.4f}".format(LA))
print("[H2O] = {:.4f}".format(H2O))
print("[ML] = {:.4f}".format(ML))
print("[EL] = {:.4f}".format(EL))
print()

# Checks that the number of moles of reactants and products are the same
reactants = C_MeOH + C_EtOH + C_LA + C_H2O
equilibrium = ML + EL + MeOH + LA + EtOH + H2O

print('Result check: \n[reactants] = {:.6f} \n[equilibrium] = {:.6f}'.format(reactants,equilibrium))
print()


# ==== Exercise 2: equimolar products ======================

# This excercise will be solved starting from the basis of setting arbitrarily
# the mass of lactic acid 85% solution. For a given amount of lactic acid, both
# methanol and ethanol amounts will be varied in a range of mass between 100
# and 1000 grams. For all these scenarios, the aim is to find all those combinations
# which result in an equimolar mixture of ethyl and metyhl lactate at the 
# equilibrium condition. 

# Changing the mass of commercial solutions of alcohols and lactic acid also
# has an impact on the concentration of water. This must be taken into account
# when going through the arrays.

# The problem admits multiple solutions, almost infinit. The more trios of 
# methanol-ethanol-lactic acid are evaluated, the more amount of solutions can
# be found. Therefore, the solutions given by this code are just a limited 
# amount of all the possible ones, since the code is written using a linspace
# array for the alcohols, and an array of for predefined mass values for acid.

# The array of mass of alcohols will be defined for a range between 200 and 800
# grams each. For each array, 200 elements will be distributed equally distant.

# The mass of lactic acid was defined as 1000 grams in exercise 1. For this 
# 2nd exercise, +/- 250 grams will also be evaluated, to show the capacity of 
# the code to solve scenarios with variable amounts of the three commercial
# mixtures.

# If the arrays contained more amount of elements, more solutions could be 
# found, but also more computational time would be required to run the code 
# and arrive to all the possible solutions.

# Creates array of mass values for lactic acid
LA_mass_values = np.array([750, 1000, 1250])

# Defines a range of mass values for each alcohol
EtOH_mass_range = np.linspace(200,800,200)
MeOH_mass_range = np.linspace(200,800,200)

tol = 0.0005        # Defines tolerance to evaluate equimolarity condition

# Creates empty arrays to later append the results 
EtOH_equimolar = np.array([])
MeOH_equimolar = np.array([])
LA_equimolar = np.array([])
ML_equimolar = np.array([])
EL_equimolar = np.array([])
react_equimolar = np.array([])
prod_equimolar = np.array([])

error_arr = np.array([])

for mass_LA in LA_mass_values:
    for mass_MeOH in MeOH_mass_range:
        for mass_EtOH in EtOH_mass_range:
            
            volume = (mass_LA / den_LA) + (mass_MeOH / den_MeOH) + (mass_EtOH / den_EtOH)
            
            C_EtOH = mass_EtOH * comp_EtOH * 1000 / (MW_EtOH * volume)
            C_MeOH = mass_MeOH * comp_MeOH * 1000 / (MW_MeOH * volume)
            C_LA = mass_LA * comp_LA * 1000 / (MW_LA * volume)
            mass_H2O = mass_EtOH * (1-comp_EtOH) + mass_LA * (1-comp_LA)
            C_H2O = mass_H2O * 1000 / (MW_H2O * volume)
            
            guess = [0.01, 0.01]
        
            x1, x2 = fsolve(model, guess)
            ML = x1          # Concentration of esters in the equilibrium
            EL = x2
            
            MeOH = C_MeOH - x1
            LA = C_LA - x1 - x2
            EtOH = C_EtOH - x2
            H2O = C_H2O + x1 + x2
        
            error = abs(ML - EL)     # Calculates difference between concentrations
            if error < tol:     # Checks if the difference is smaller than the tolerance 
            
            # For each case in which the condition of equimolarity is fulfilled,
            # the mass of the commercial mixtures and the concentration of the
            # esters are appended to the solution arrays. 
            
                EtOH_equimolar = np.append(EtOH_equimolar, mass_EtOH)
                MeOH_equimolar = np.append(MeOH_equimolar, mass_MeOH)
                LA_equimolar = np.append(LA_equimolar, mass_LA)
                ML_equimolar = np.append(ML_equimolar, ML)
                EL_equimolar = np.append(EL_equimolar, EL)
            
                error_arr= np.append(error_arr, error)
                
                react_equimolar = np.append(react_equimolar, C_MeOH + C_EtOH + C_LA + C_H2O)
                prod_equimolar = np.append(prod_equimolar, ML + EL + MeOH + LA + EtOH + H2O)

# Creates a txt file with all the combinations found
with open('ex2_project1.txt', 'w') as file:
    file.write('Mass of LA 85%\t| Mass of EtOH 96%\t| Mass of MeOH 100%\t| [ML]eq\t| [EL]eq\t| Conc. diff.\t| [react]\t| [prod]\n')
    file.write('[g]\t\t| [g]\t\t\t| [g]\t\t\t| [mol/L]\t| [mol/L]\t| [mol/L]\t| [mol/L]\t| [mol/L]\n\n')
    for LA, EtOH, MeOH, ML, EL, error, react, prod in zip(LA_equimolar,EtOH_equimolar, MeOH_equimolar, ML_equimolar, EL_equimolar, error_arr, react_equimolar, prod_equimolar):
        file.write('{}\t\t| {:.2f}\t\t| {:.2f}\t\t| {:.3f}\t\t| {:.3f}\t\t| {:.4e}\t| {:.3f}\t| {:.3f}\n'.format(LA,EtOH,MeOH,ML,EL,error, react, prod))
        
# Plot generation
# Creates empty arrays to later append the solutions found for each LA mass value
plot_array_x1 = np.array([])
plot_array_x2 = np.array([])
plot_array_x3 = np.array([])

plot_array_y1 = np.array([])
plot_array_y2 = np.array([])
plot_array_y3 = np.array([])

plt.title("Ethanol and methanol amounts that yield [EL]=[ML]")
plt.xlabel("Mass of ethanol 96%wt")
plt.ylabel("Mass of methanol 100%wt")

# For loop that appends the solutions for each LA mass value
i=0
for LA_conc_equim in LA_equimolar:
    if LA_conc_equim == LA_mass_values[0]:
        plot_array_x1 = np.append(plot_array_x1, EtOH_equimolar[i])
        plot_array_y1 = np.append(plot_array_y1, MeOH_equimolar[i])
    elif LA_conc_equim == LA_mass_values[1]:
        plot_array_x2 = np.append(plot_array_x2, EtOH_equimolar[i])
        plot_array_y2 = np.append(plot_array_y2, MeOH_equimolar[i])
    else:
        plot_array_x3 = np.append(plot_array_x3, EtOH_equimolar[i])
        plot_array_y3 = np.append(plot_array_y3, MeOH_equimolar[i])
    i += 1

# Plots the obtained arrays
plt.plot(plot_array_x1, plot_array_y1, 'ro-', label="750g LA 85%wt")
plt.plot(plot_array_x2, plot_array_y2, 'bo-', label="1000g LA 85%wt")
plt.plot(plot_array_x3, plot_array_y3, 'go-', label="1250g LA 85%wt")

plt.legend(loc="upper left")

# Running the code may delay some time (around 15 seconds) due to the number
# of iterations that the nested for-loop has to make in exercise 2