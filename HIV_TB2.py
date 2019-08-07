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
alpha = 714

Beta1 = 4.3  # General Transmission rate
Beta2 = 0.051  # HIV Transmission rate, 0.055, 0.08
enn = 1.02
NatDeathRate = 1/70
Tee2 = Tee3 = 2 #TreatRate for Individuals with Active TB 2 yr-1
ActiveDeathRate = 1/8 #yr
Delta =1.03
TransmitRateRecovTB = 0.9 #Beta11

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
N = 11745189  # Total Population
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
S_N = zeros((1, run_time))
TI = zeros((1, run_time))
TR = zeros((1, run_time))
HI= zeros((1, run_time))
A = zeros((1, run_time))
THI = zeros((1, run_time))
THR = zeros((1, run_time))
AT = zeros((1, run_time))
RTH = zeros((1, run_time))

lamdaTB = (Beta1 *(TI + THI + AT))/N
lamdaHIV = (Beta2 *(HI+ THI + AT + RTH + (enn *(A+AT))))/N

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
        S_N[0, j] = S_N[0, i] - ((alpha - lamdaTB * S_N[0, i]) -(lamdaHIV * S_N[0, i]) - (NatDeathRate*S_N[0, i]))

        # Calculate Active Infected TB at time t+1
        TI[0, j] = TI[0, i] - ((Tee2 + ActiveDeathRate + NatDeathRate + (Delta*lamdaHIV))*TI[0,i])

        # Calculate Recovered at time t+1
        TR[0, j] = TR[0, i] + (Tee2 *TI[0, i] )-((TransmitRateRecovTB*lamdaTB) +lamdaHIV+NatDeathRate)*TR[0,i]



        # Calculate HIV Infected  Individuals (With No AIDS)at time t+1
        HI[0, j] = HI[0, i] + (lamdaHIV * S_N[0, i]) - (( HIVtoAIDSRate1 + TBInfectRateForHI +NatDeathRate)* HI[0, i] +
                                                        (AIDSTreatRate*AT[0, i]) +(lamdaHIV*TR[0, i])


        # Calculate HIV Infected  Individuals (With  AIDS)t at time t+1
        A[0, j] = A[0, i] + ((HIVtoAIDSRate1*HI[0, i]) -(AIDSTreatRate*A[0, i])-(NatDeathRate+AIDSDeathRate)*A[0, i])

        # Calculate Coinfected Individuals (Pre-AIDS) at time t+1
        THI [0, j] = THI[0, i] + (((HIVInfectRateforActiveTB * TI[0, i])+ (TBInfectRateForHI*HI[0, i] + Alpha2*AT[0, i]
                                                                           -(Tee3 +  Pee2+ NatDeathRate + TBDeathRateforCoinfect)*THI[0, i]

        # Calculate TB_recovered Individuals (Co-infected, Pre-AIDS) at time t+1
        RTH[0, j] = RTH[0, i] + ((TBTreatRateForCoInfect *THI[0, i]) +(ReinfectRateCoinfect*lamdaTB + Pee3 + NatDeathRate)*RTH[0, i]



        # Calculate HIV-infected Individuals (Co-infected, AIDS) at time t+1
        AT[0, j] = AT[0, i] + ((Pee2 * THI[0, i]) + (Pee3*RTH[0, i])- ((Alpha2 + NatDeathRate +AidsTBDeathRate ) * AT[0, i])
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

