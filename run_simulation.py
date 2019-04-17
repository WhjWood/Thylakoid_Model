import Thylakoid_model_2 as Tm
import pickle
import time

GRANA_SIZE = 170 # width of grana, nm
DATE = "10_04_19"  # a reference date in which the simulations are run 
EXPERIMENT = "SI"   # a reference for which experiment is being run
Number_of_iterations = 1000 # number of Monte Carlo steps
Stacking_Interaction_Energy = 4 # stacking interaction strength, kT
LHCII_Binding_Interaction_Energy = 2 # intralayer LHCII interaction strength, kT

t0 = time.time()
POPULATION1, POPULATION2 = Tm.Run_Simulation(GRANA_SIZE,DATE,EXPERIMENT,Number_of_iterations,Stacking_Interaction_Energy,LHCII_Binding_Interaction_Energy)
print("Completed in ", time.time()-t0, "seconds")

D = Tm.Nearest_Neighbour_Distances(POPULATION1, "LHCII")