'''
A program to model treatment options of HIV/TB Coinfection and Estimate the Geospatial Allocation of Resources Among West African Countries,
Created by Sampson Akwafuo
Oct 2019

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
NatDeathRate =  1/70
Tee2 = Tee3 = 2 #TreatRate for Individuals with ActiveTB 2 yr-1 Tee3 -TBTreatRateForCoInfect
ActiveDeathRate = 1/8 #yr
Delta =1.03
TransmitRateRecovTB = 0.9 #Beta11
Pee1 = 0.1 #HIVtoAIDSRate1
Pee2 = 0.25 # 0.25yr-1 AIDSProgressionRateforCoinfect
Pee3 = 0.125 #0.125yr-1 AIDSProgressionRateforReinfect, after recovCoinfect
mew = 1.07 #TBInfectRateForHI
Alpha1 = Alpha2 = 0.33 #HIV Treatment rate for AIDS stage
dA= 0.3 # AIDSDeathRate 0.3 Death rate for AIDS stage dA
dT =  1/8# dT TBDeathRateforCoinfect
dTA = 0.33 #0.33yr-1 AidsTBDeathRate
Beta22 = 1.1 #ReinfectRateCoinfect




# Set the initial number of people in the model
N_0 = 101558800  # Total Population of 15-65 in Nigeria= 101558800
# N_SL_0      = 6092075   # To be used for cross-country modeling
# N_LI_0      = 4294077


# Set the initial proportion of the sub-groups
S_N_0 = 0.84  # Initial proportion of Susceptible
TI_0 = 0.09  # Initial proportion of TB-Infected
TR_0 = 0.04  # Initial proportion of TB_Recovered
HI_0 = 0.02  # Initial proportion of HIV_Infected, No AIDS
A_0 = 0.01  # Initial proportion of AIDS individuals
THI_0 = 0.00  # Initial proportion of Co-infected
THR_0 = 0.00  # Initial proportion of Co-infected, TB_Recovered
AT_0 = 0.00  # Initial proportion of Co-infected, Active TB and AIDS Symptoms


# Pre allocate vectors for model variables...
N_Tot = zeros((1, run_time))
S_N = zeros((1, run_time))
TI = zeros((1, run_time))
TR = zeros((1, run_time))
HI= zeros((1, run_time))
A = zeros((1, run_time))
THI = zeros((1, run_time))
THR = zeros((1, run_time))
AT = zeros((1, run_time))
#RTH = zeros((1, run_time))

lamdaTB = (Beta1 *(TI + THI + AT))/N_0
lamdaHIV = (Beta2 *(HI+ THI + AT + THR + (enn *(A+AT))))/N_0

# Print which iteration we've run
print('Running model...')
# Open a CSV file to write to
openfile = open("HIV_TB_Modeling_Results.csv", "w")
excelwriter = csv.writer(openfile)
# Write the first row of the output file (a column for each population group)
# excelwriter.writerow(["Time", "Susceptible_Popn", "Active TB", "TB_Recovered", "HIV_Only", "AIDS_Symptoms", "Coinfected_No_AIDS", "Coinfected_AIDSnTB"])

# This loop simulates the model
for j in range(run_time):

    # On the first iteration, do this part
    if j == 0:
        # Define initial conditions...
        # Total people in each population group
        #S_N[0, j] = S_N_0

        # Susceptibles
        S_N[0, j] = N_0 * S_N_0
        # TB_Infected
        TI[0, j] = N_0 * TI_0
        # Hospitalized
        TR[0, j] = N_0 * TR_0
        # Funeral Wait List
        HI[0, j] = N_0 * HI_0
        # Recovered
        A[0, j] = N_0 * A_0
        # Vaccinated
        THI[0, j] = N_0 * THI_0
        THR[0, j] = N_0 * THR_0
        AT[0, j] = N_0 * AT_0
        # Write a new row for the initial conditions (same format as title row)
        excelwriter.writerow([2019 + j, S_N[0, j], TI[0, j], TR[0, j], HI[0, j], A[0, j], THI[0, j], THR[0, j], AT[0, j]])
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
        HI[0, j] = HI[0, i] + (lamdaHIV * S_N[0, i]) - (( Pee1 + (mew*lamdaTB) +NatDeathRate)* HI[0, i] +
                                                        (Alpha1*AT[0, i]) +(lamdaHIV*TR[0, i]))


        # Calculate HIV Infected  Individuals (With  AIDS)t at time t+1
        A[0, j] = A[0, i] + ((Pee1*HI[0, i]) -(Alpha1*A[0, i])-(NatDeathRate + dA)*A[0, i])

        # Calculate Coinfected Individuals (Pre-AIDS) at time t+1
        THI [0, j] = THI[0, i] + (Delta * lamdaHIV * TI[0, i]) + (mew * lamdaTB*HI[0, i]) + (Alpha2*AT[0, i]) \
                     -((Tee3 +  Pee2+ NatDeathRate + dT)*THI[0, i])

        # Calculate TB_recovered Individuals (Co-infected, Pre-AIDS) at time t+1
        THR[0, j] = THR[0, i] + (( Tee3 *THI[0, i]) -((Beta22*lamdaTB + Pee3 + NatDeathRate)*THR[0, i]))


        # Calculate HIV-infected Individuals (Co-infected, AIDS) at time t+1
        AT[0, j] = AT[0, i] + ((Pee2 * THI[0, i]) + (Pee3*THR[0, i])- ((Alpha2 + NatDeathRate + dTA) * AT[0, i]))

        # Print which iteration we've run
        print('Iteration %i completed' % j)

        # ---------- end-else ---------------#

    # Determine population size
    N_Tot[0, j] = S_N[0, j]+ TI[0, j]+ TR[0, j]+ HI[0, j]+ A[0, j]+ THI[0, j]+ THR[0, j]+ AT[0, j]

    # Write a new row on each iteration
    excelwriter.writerow([2020 + j, S_N[0, j], TI[0, j], TR[0, j], HI[0, j], A[0, j], THI[0, j], THR[0, j], AT[0, j]])
# ---------- end-loop ---------------#

# Print which iteration we've run
#print('Iteration %i completed' % j)

# Close the csv file
openfile.close()

