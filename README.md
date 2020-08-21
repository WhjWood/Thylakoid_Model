# Thylakoid_Model
WhjWood 05/03/2019
A model of the Spinach thylakoid membrane for Monte Carlo simulations.

Main code in Thylakoid_model_2.py
Simulates the interactions of PSII, LHCII and PSI in the thylakoid membrane.
Parametrised using AFM,SIM and published data.

The simplest way to run simulations is to use the streamlit app

- A guide for how to do this can be found here:
<in progress>

 How to use Thylakoid_model_2.py in Python

in the if __name__== '__main__' clause (from line 1621):

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


Running the simulation

POPULATION1, POPULATION2 = Run_Simulation(GRANA_SIZE,DATE,EXPERIMENT,Number_of_iterations,Stacking_Interaction_Energy,LHCII_Binding_Interaction_Energy,PSI_interaction_energy)
- This runs the simulation using the constants defined above.
- 10^6 Monte Carlo steps are used for equilibration before data is collected. The energy versus number of iterations for monitoring purposes is stored in Equilibration_Energy.csv in the simulation folder.
- You will be kept updated with the progress of the simulation eg.
Equilibration 40% complete
 LHCII in grana =  0.7293447293447294 / 0.7293447293447294
 LHCII NN =  8.60058381242 / 8.60058381242
- All data will be stored in the experiment folder in a folder called Data
- takes 4-10 days on a standard PC


Analysing the data

Run_analysis(GRANA_SIZE,DATE,EXPERIMENT)
- This should be run after the simulation.
- creates 2 files in the simulation folder: Results.csv contains average nearest neighbour distances between PSII, LHCII and PSI, the percentage of LHCII in the grana, and the densities of the grana and stromal lamellae regions. 
Nearest_Neighbour_Results records individual nearest neighbour distances in order to create distribution histograms.

Run_graph_antenna_analysis(GRANA_SIZE,DATE,EXPERIMENT,PSII=True,PSI=True)
- This creates the chlorophyll networks which are used to measure the connectivity of PSII and the antenna sizes of PSII and PSI for topological analysis.
- outputs are saved in the simulation folder
- you can choose which photosystem to analyse eg. PSII = True runs the analysis for PSII
- 1-2 days running time



   Below is an example of an entire experiment plus analysis

GRANA_SIZE = 170 # width of grana, nm.
Number_of_iterations = 11000001 # number of Monte Carlo steps, Note that data is only collected after 10M iterations.
DATE = "08_07_19"  # a reference date in which the simulations are run.  
EXPERIMENT = "Test_170"   # a reference for which experiment is being run.
Stacking_Interaction_Energy = 4 # stacking interaction strength, kT (default = 4).
LHCII_Binding_Interaction_Energy = 2 # intralayer LHCII interaction strength, kT (default = 2).
PSI_interaction_energy = 0 # PSI - LHCII interaction strength, kT (default = 0, SII = 2).

POPULATION1, POPULATION2 = Run_Simulation(GRANA_SIZE,DATE,EXPERIMENT,Number_of_iterations,Stacking_Interaction_Energy,LHCII_Binding_Interaction_Energy,PSI_interaction_energy)

Run_analysis(GRANA_SIZE,DATE,EXPERIMENT)
Run_graph_antenna_analysis(GRANA_SIZE,DATE,EXPERIMENT,PSII=False,PSI=True)




