#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 10:12:52 2020

@author: williamwood

Streamlit App for the Thylakoid_Model

"""
import math, random, time
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as dist
import pandas as pd
import pickle
import networkx as nx
import streamlit as st


import Thylakoid_model_2 as TM

# generally constants
t0 = time.time()
Experiment_Setting = st.sidebar.selectbox("Experiment", ["State I","State II","Custom"])
if Experiment_Setting == "State I":
    EXPERIMENT = "SI"
    GRANA_SIZE = 190
    Stacking_Interaction_Energy = 4 # stacking interaction strength, kT (default = 4). 
    LHCII_Binding_Interaction_Energy = 2 # intralayer LHCII interaction strength, kT (default = 2).
    PSI_interaction_energy = 0 # PSI - LHCII interaction strength, kT (default = 0, SII = 3).
elif Experiment_Setting == "State II":
    EXPERIMENT = "SII"
    GRANA_SIZE = 170
    Stacking_Interaction_Energy = 0 # stacking interaction strength, kT (default = 4). 
    LHCII_Binding_Interaction_Energy = 0 # intralayer LHCII interaction strength, kT (default = 2).
    PSI_interaction_energy = 3 # PSI - LHCII interaction strength, kT (default = 0, SII = 3.
    
elif Experiment_Setting == "Custom":
    grana_radius = st.sidebar.radio("Grana Radius (nm)", [190,170])
    GRANA_SIZE = int(grana_radius)
    Stacking_Interaction_Energy = st.sidebar.text_input("Stacking interaction (KT)",0)
    try:
        Stacking_Interaction_Energy = float(Stacking_Interaction_Energy)
    except:
        st.write("Unable to convert stacking interaction energy to float")
        
    LHCII_Binding_Interaction_Energy = st.sidebar.text_input("LHCII Lateral interaction (KT)",0)
    try:
        LHCII_Binding_Interaction_Energy = float(LHCII_Binding_Interaction_Energy)
    except:
        st.write("Unable to convert lateral interaction energy to float")
    
    PSI_interaction_energy = st.sidebar.text_input("PSI-LHCII interaction (KT)",0)
    
    try:
        PSI_interaction_energy = float(PSI_interaction_energy)
    except:
        st.write("Unable to convert PSI-LHCII interaction energy to float")
    EXPERIMENT = "Custom_" + str(GRANA_SIZE )+ "_" + str(Stacking_Interaction_Energy)+"_"+str(LHCII_Binding_Interaction_Energy)+"_"+str(PSI_interaction_energy)


Run_Analysis = st.sidebar.checkbox("Run Post Simulation Analysis")
Run_Network_Analysis = st.sidebar.checkbox("Run Network Analysis")
 # width of grana, nm.
Number_of_iterations = 11000001 # number of Monte Carlo steps, Note that data is only collected after 10M iterations.

today = date.today()
DATE = str(date.today().year) + "_" + str(date.today().month)+"_" + str(date.today().day)# a reference date in which the simulations are run.



if st.button("Run Simulation"):
    with st.spinner('Running Simulation'):
        POPULATION1, POPULATION2 = TM.Run_Simulation(GRANA_SIZE,DATE,EXPERIMENT,Number_of_iterations,Stacking_Interaction_Energy,LHCII_Binding_Interaction_Energy,PSI_interaction_energy)

    if Run_Analysis:
        with st.spinner('Running Post-Simulation Analysis'):
            TM.Run_analysis(GRANA_SIZE,DATE,EXPERIMENT)
    
    if Run_Network_Analysis:
        with st.spinner('Running Networn Analysis'):
                TM.Run_graph_antenna_analysis(GRANA_SIZE,DATE,EXPERIMENT,PSII=True,PSI=True)