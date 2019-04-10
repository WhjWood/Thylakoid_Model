import Thylakoid_model_2 as Tm
import pickle
import time

t0 = time.time()
POPULATION1, POPULATION2 = Tm.Run_Simulation(170,"11_05_18","FREE_LHCII",1000)
print(time.time()-t0)