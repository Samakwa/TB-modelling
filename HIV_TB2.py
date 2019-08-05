'''
A program to model treatment options of HIV/TB Coinfection and Estimate the Geopatial Allocation of Resources Among West African Countries,

First Part is with Guinea

Created by Sampson Akwafuo
Oct 2017

Edited by: Sultanah Alshammari
Edited by:
'''

# Import zeros function
from numpy import zeros
import numpy as np
# Import CSV usability
import csv
import pandas as pd

#zeros = np.zeros((1, 200), np.float64)

# Set model simulation time (   )
run_time = 200

# Define transmission parameters
BetaGen_rate = 0.2275  # General Transmission rate
BetaHosp_rate = 0.4269  # Transmission rate during Hospitalization
BetaFuneral_rate = 0.0445  # Transmission rate during funeral
Vac_rate = 0.005  # Vaccination rate. Varies; 0, 0.005 and 0.01
VacEffic = 0.5487  # Efficacy of Vaccination
GammaHosp = 0.309  # Time lag before Hospitalization 1/GammaHosp = 3.24 days
Theta1 = 0.197  # Fraction of Infection Popn going to Hospital
GammaInfected = 0.067  # Typical Infection Duration, 1/GammaInfected = 15days
GammaDeath = 0.0963  # Time from Infection to death  1/Gammad =10.38 days
Delta1 = 0.750  # Unhospitalized fatality rate
GammaDeathHosp = 0.1597  # time from Hospitalization to death, 1/gammaDH = 6.26
Delta2 = 0.6728  # Hospitalized fatality rate
GammaInfecHosp = 0.0630  # Time from Infection to recovery, 1/GammaIH =15.88
GammaFuneral = 0.3072  # Traditional Funeral Duration, 1/GammaF = 3.255
R0V_GU = 1.52  # Initial Reproduction number for Guinea

# Set the initial number of people in the model
N_GU_0 = 11745189  # Total Population
# N_SL_0      = 6092075   # To be used for cross-country modeling
# N_LI_0      = 4294077


# Set the initial proportion of the sub-groups
S_GU_0 = 0.84  # Initial proportion of Susceptible
I_GU_0 = 0.09  # Initial proportion of Infected
H_GU_0 = 0.04  # Initial proportion of Hospitlaized
F_GU_0 = 0.02  # Initial proportion of Funeral Wait
R_GU_0 = 0.01  # Initial proportion of Recovered
V_GU_0 = 0.00  # Initial proportion of Vaccinated

# Pre allocate vectors for model variables...
N_GU = zeros((1, run_time))
S_GU = zeros((1, run_time))
V_GU = zeros((1, run_time))
I_GU = zeros((1, run_time))
H_GU = zeros((1, run_time))
F_GU = zeros((1, run_time))
R_GU = zeros((1, run_time))
R0V_GU = zeros((1, run_time))

# Print which iteration we've run
print('Running model...')
# Open a CSV file to write to
openfile = open("HIV_TB_Modeling_Results.csv", "w")
excelwriter = csv.writer(openfile)
# Write the first row of the output file (a column for each population group)
# excelwriter.writerow(["Time", "Susceptible_GU", "Vaccinated_GU", "Infected_GU", "Hospitalized_GU", "FuneralWait_GU", "Recovered_GU"])

# This loop simulates the model
for j in range(run_time):

    # On the first iteration, do this part
    if j == 0:
        # Define initial conditions...
        # Total people in each population group
        N_GU[0, j] = N_GU_0
        # Susceptibles in each country
        S_GU[0, j] = N_GU_0
        # Susceptibles
        S_GU[0, j] = N_GU_0 * S_GU_0
        # Infected
        I_GU[0, j] = N_GU_0 * I_GU_0
        # Hospitalized
        H_GU[0, j] = N_GU_0 * H_GU_0
        # Funeral Wait List
        F_GU[0, j] = N_GU_0 * F_GU_0
        # Recovered
        R_GU[0, j] = N_GU_0 * R_GU_0
        # Vaccinated
        V_GU[0, j] = N_GU_0 * V_GU_0
        # Write a new row for the initial conditions (same format as title row)
        excelwriter.writerow([2017 + j, S_GU[0, j], V_GU[0, j], I_GU[0, j], H_GU[0, j], F_GU[0, j], R_GU[0, j]])
        # ---------- end-if ---------------#
    # On any other iteraction do this part
    else:
        # Use the previous time point, were everything has been calculated
        i = j - 1
        # Balance populations...
        # Calculate difference equations...
        # Calculate Susceptibles at time t+1
        S_GU[0, j] = S_GU[0, i] - (((BetaGen_rate * S_GU[0, i] * I_GU[0, i]) + (
                    BetaHosp_rate * S_GU[0, i] * H_GU[0, i]) + (BetaFuneral_rate * S_GU[0, i] * F_GU[0, i])) / N_GU[
                                       0, i])
        # Calculate Vaccinated Population at time t+1
        V_GU[0, j] = V_GU[0, i] + ((Vac_rate * S_GU[0, i]) - ((
                    ((BetaGen_rate * I_GU[0, i]) + (BetaHosp_rate * H_GU[0, i]) + (BetaFuneral_rate * F_GU[0, i])) * (
                        VacEffic * S_GU[0, i] * V_GU[0, i]))) / N_GU[0, i])
        # Calculate Infected at time t+1
        I_GU[0, j] = I_GU[0, i] + ((GammaHosp * Theta1 + (GammaInfected * (1 - Theta1) * (1 - Delta1)) + (
                    GammaDeath * (1 - Theta1) * Delta1)) * I_GU[0, i])
        # Calculate Number of Hospitalized Individuals at time t+1
        H_GU[0, j] = H_GU[0, i] + (
                    (GammaHosp * Theta1 * I_GU[0, i]) - ((GammaDeathHosp * Delta2) + (GammaInfecHosp * (1 - Delta2))) *
                    H_GU[0, i])
        # Calculate Number of Individuals on Funeral_Waiting List at time t+1
        F_GU[0, j] = F_GU[0, i] + (
                    (GammaDeath * (1 - Theta1) * Delta1 * I_GU[0, i]) + (GammaDeathHosp * Delta2 * H_GU[0, i]) - (
                        GammaFuneral * F_GU[0, i]))
        # Calculate Number of Recovered Individuals at time t+1
        R_GU[0, j] = R_GU[0, i] + (((GammaInfected * (1 - Theta1) * (1 - Delta1)) * I_GU[0, i]) + (
                    GammaDeathHosp * (1 - Delta2) * H_GU[0, i]) - (GammaFuneral * F_GU[0, i]))
        # Calculate Basic Reproduction Number, R0, When Vaccination is Introduced
        R0V_GU[0, j] = R0V_GU[0, i] + (Vac_rate * R0V_GU[0, i])

        # Print which iteration we've run
        print('Iteration %i completed' % j)

        # ---------- end-else ---------------#

    # Determine population size
    N_GU[0, j] = S_GU[0, j] + V_GU[0, j] + I_GU[0, j] + H_GU[0, j] + F_GU[0, j] + R_GU[0, j]

    # Write a new row on each iteration
    excelwriter.writerow([2020 + j, S_GU[0, j], V_GU[0, j], I_GU[0, j], H_GU[0, j], F_GU[0, j], R_GU[0, j]])
# ---------- end-loop ---------------#

# Print which iteration we've run
#print('Iteration %i completed' % j)

# Close the csv file
openfile.close()

