# Thylakoid_Model
WhjWood 05/03/2019
A model of the Spinach thylakoid membrane for Monte Carlo simulations.
### In Construction ###
Main code in Thylakoid_model_2.py
Simulates the interactions of PSII, LHCII and PSI in the thylakoid membrane.
Parametrised using AFM,SIM and published data.

How to use Thylakoid_model_2.py

in the if __name__== '__main__' clause (from line 1412):

set up the constants:

GRANA_SIZE = 170 # width of grana, nm.
- there are premade model files for thylakoids with 170 nm and 190 nm. Therefore, these number should be used.
if you wish to use a different size thylakoid, you will have to create a new population (model) using the function create_initial_population(GRANA_RADIUS).

Number_of_iterations = 11000001 # number of Monte Carlo steps, Note that data is only collected after 10M iterations.
- this should be left as is (I may remove this from the global scope)

DATE = "05_07_19"  # a reference date in which the simulations are run.
EXPERIMENT = "Stacking_Only_170"   # a reference for which experiment is being run.
- DATE and EXPERIMENT help store the data. The data and results are stored in a folder [EXPERIMENT] + "_" + [DATE].
eg. "Stacking_Only_170_05_07_19"

Stacking_Interaction_Energy = 4 # stacking interaction strength, kT (default = 4).
LHCII_Binding_Interaction_Energy = 0 # intralayer LHCII interaction strength, kT (default = 2).
PSI_interaction_energy = 0 # PSI - LHCII interaction strength, kT (default = 0, SII = 2).

- These are the interaction strengths for the Monte Carlo simulation.



