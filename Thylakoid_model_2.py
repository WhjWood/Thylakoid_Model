import math, random, time
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as dist
import pandas as pd
import pickle
import networkx as nx


##WHJWood Thylakoid model 2.0 05/03/19


global Complexes
global Complexes_Chl
Complexes_Chl = pickle.load(open("Complexes_Chl.p",'rb')) # for antenna analysis
Complexes = pickle.load(open("Complexes.p", 'rb'))


#Particle class definition

class Particle(object):
    
    # class variables 
    
    # offsets of PSII-bound LHCIIs
    C2S2M2_LHCIIa = np.matrix([[-10.0,-8.0]]).T
    C2S2M2_LHCIIb = np.matrix([[-10.0,1.0]]).T
    C2S2M2_LHCIIc = np.matrix([[11.0,8.0]]).T
    C2S2M2_LHCIId = np.matrix([[11.0,0.0]]).T 
    C2S2_LHCIIa = np.matrix([[-7.0,7.0]]).T
    C2S2_LHCIIb = np.matrix([[7.0,-7.0]]).T
    PSI_LHCII = np.matrix([[-8.0,-7.0]]).T
    
    
    def __init__(self,x,y,theta, particle_type):            
        self.location = np.matrix([[int(x)],[int(y)]])
        self.rotation = float(theta)
        self.Ptype = str(particle_type)
        self.Pmatrix = self.particle_matrix()

    
    def particle_matrix(self):
        """takes the list [translation, rotation] and produces the particle matrix for P""" 
        X1 = self.rotate(self.rotation,Complexes[self.Ptype],0)
        X2 = self.translate(X1,self.location)
        return self.remove_repeat_cols(X2).astype(int)
    
    def remove_repeat_cols(self,A,C=10000):
        """returns only unique columns of A"""
        return np.matrix(np.unique(A,axis=1))
#         if A.shape[0] == 2:
# 
#             MIN = np.min(A)
#             MIN_M = np.repeat(np.matrix([[MIN],[MIN]]),A.shape[1],axis=1)
#             unique_A = pd.unique(self.unique_2D_mapping(A-MIN_M))
#             X = np.matrix(unique_A%C)
#             Y = np.matrix(MIN+ (unique_A-X)/C)
#             X+=MIN
#     
#             new_A = np.concatenate((Y,X),axis=0)
#             return new_A
#         else:
#            print("Something went wrong. Matrix has incorrect dimensions")
    
    def move_particle(self,d_translation,d_rotation):
        self.location += d_translation
        self.rotation += d_rotation
    
    
    def bound_LHCII(self):
        """returns the coordinates of bound LHCIIs in matrix form"""
        if self.Ptype == "C2S2M2": 
            LHCa = self.rotate(self.rotation,Particle.C2S2M2_LHCIIa,0)+self.location
            LHCb = self.rotate(self.rotation,Particle.C2S2M2_LHCIIb,0)+self.location
            LHCc = self.rotate(self.rotation,Particle.C2S2M2_LHCIIc,0)+self.location
            LHCd = self.rotate(self.rotation,Particle.C2S2M2_LHCIId,0)+self.location
            return np.concatenate((LHCa,LHCb,LHCc,LHCd), axis=1)
        
        elif self.Ptype == "C2S2":
            LHCa = rotate(self.rotation,Particle.C2S2_LHCIIa,0)+self.location
            LHCb = rotate(self.rotation,Particle.C2S2_LHCIIb,0)+self.location
            return np.concatenate((LHCa,LHCb), axis=1)
        
        elif self.Ptype == "LHCII":
            return self.location
        
        else:
            pass
    
    def PSI_site(self):
        return  self.rotate(self.rotation,Particle.PSI_LHCII)+self.location
    
    
    def bound_Chorophylls(self):
        """Returns a matrix containing locations of the chlorophylls""" 
        X1 = self.rotate(self.rotation,Complexes_Chl[self.Ptype],0)
        X2 = self.translate(X1,self.location)
        return self.remove_repeat_cols(X2).astype(int)
        

    
    @staticmethod
    def rotate(theta,A, axis=0):
        # might remove axis later
        return np.round(np.matrix([[np.cos(theta),-np.sin(theta)],[np.sin(theta), np.cos(theta)]])*np.matrix(A))
    

    @staticmethod
    def translate(T, A):
        return A+np.repeat(T,A.shape[1],axis=1)

    @staticmethod
    def unique_2D_mapping(A,C=10000):
        """converts a 2D matrix to a 1d vector"""
        return np.array(A[0,:]+C*A[1,:])[0] # this maps x and y onto a unique z and returns array 

class SudoParticle(Particle):
    # Same as particle but takes a location matrix
    # used for particle steps
    def __init__(self,Location,theta, particle_type):            
        self.location = Location
        self.rotation = float(theta)
        self.Ptype = str(particle_type)
        self.Pmatrix = self.particle_matrix()


### Functions for overlap detection ###    


def unique_2D_mapping(A,C=10000):
        """converts a 2D matrix to a 1d vector"""
        return np.array(A[0,:]+C*A[1,:])[0] # this maps x and y onto a unique z and returns array 

def N_unique(A):
    """returns the number of unique columns in matrix A"""
    MIN = np.min(A)
    return pd.unique(unique_2D_mapping(A-MIN)).shape[0]

def remove_repeat_cols(A):
    return np.matrix(np.unique(A,axis=1))


def collision2(X,P):
    return np.unique(unique_2D_mapping(np.concatenate((X,P),axis=1))).shape[0] < X.shape[1]+P.shape[1] 

def collision3(X,P):
    return N_unique(np.concatenate((X,P),axis=1)) < X.shape[1]+P.shape[1] 


### Utility functions ###

def Save_Population(population, path):
    """Saves population (Particle class) data as csv in form x,y,theta,Ptype"""
    Population = pd.DataFrame(columns = ["x","y","theta","Ptype"])
    for particle in population:
        xy_coordinate = particle.location.T.tolist()[0]
        theta = particle.rotation
        PTYPE = particle.Ptype
        row = pd.DataFrame([[xy_coordinate[0],xy_coordinate[1],theta,PTYPE]],columns = ["x","y","theta","Ptype"])
        Population = Population.append(row,ignore_index=True)
    Population.to_csv(path)

def Load_Population(path):
    """Loads data from csv (x,y,theta,Ptype) and returns a list of Particle class instances"""
    POP_df = pd.read_csv(path) # dataframe containing particle data
    population = []
    for index, row in POP_df.iterrows():        
        x = int(row['x'])
        y = int(row['y'])
        theta = float(row['theta'])
        PTYPE = str(row['Ptype'])
        particle = Particle(x,y,theta,PTYPE)
        population.append(particle)
    return population

def Plot_Population(Population,DirName):
    #plotting
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_xlabel("Distance (nm)")
    plt.axis([-300,300,-300,300])
    plt.axis('equal')
    plt.axes().set_aspect('equal')
    
    Ptype_list = [p.Ptype for p in Population]
    
    if ("C2S2M2" in Ptype_list) or ("C2S2" in Ptype_list) or ("PSII_mono" in Ptype_list):
        MPSII_4_plot = np.concatenate(tuple(p.Pmatrix for p in Population if (p.Ptype=="C2S2M2") or (p.Ptype=="C2S2") or (p.Ptype=="PSII_mono")), axis=1)
        plt.scatter(np.array(MPSII_4_plot[0,:]),np.array(MPSII_4_plot[1,:]),s=0.1,color='g')
    if "LHCII" in Ptype_list:
        MLHCII_4_plot = np.concatenate(tuple(p.Pmatrix for p in Population if (p.Ptype=="LHCII")), axis=1)
        plt.scatter(np.array(MLHCII_4_plot[0,:]),np.array(MLHCII_4_plot[1,:]),s=0.1,color='c')
    if "B6F" in Ptype_list:
        MB6F = np.concatenate(tuple(p.Pmatrix for p in Population if (p.Ptype=="B6F")), axis=1)
        plt.scatter(np.array(MB6F[0,:]),np.array(MB6F[1,:]),s=0.1,color='m')
    if "PSI" in Ptype_list:
        MPSI = np.concatenate(tuple(p.Pmatrix for p in Population if (p.Ptype=="PSI")), axis=1)
        plt.scatter(np.array(MPSI[0,:]),np.array(MPSI[1,:]),s=0.1,color='r')
    if "ATP" in Ptype_list:
        MATP = np.concatenate(tuple(p.Pmatrix for p in Population if (p.Ptype=="ATP")), axis=1)
        plt.scatter(np.array(MATP[0,:]),np.array(MATP[1,:]),s=0.1,color='k')
    
    plt.savefig(DirName+".png",format='png',dpi=1000)
    plt.clf()
    plt.cla()
    plt.close('all')

def list_to_matrix(coords):
    return np.matrix(coords).T


def rotate(theta,A, axis):
        return np.round(np.matrix([[np.cos(theta),-np.sin(theta)],[np.sin(theta), np.cos(theta)]])*np.matrix(A))


def population_matrix(X):
    """X is a list of Particle objects"""
    return np.concatenate([x.Pmatrix for x in X],axis=1)

 
def coordinate_matrix(X):
    return np.concatenate([np.matrix(x) for x in X],axis=1)

 
def pdb_to_matrix(file_name):
    particle_file = open(file_name,'r')
    COORDS = []
    for line in particle_file:
        if line[0:4] == "ATOM" or line[0:6] == "HETATM" or line[0:6] == "ANISOU":
            # divide by 10 to get nm or leave as is for Angstroms
            COORDS.append((int(round(float(line[30:38])/10)),int(round(float(line[38:46])/10))))#,int(float(line[46:54])/10.0)))

    X = remove_repeat_cols(list_to_matrix(COORDS))
    return (X-X.mean(1)).astype(int)


def create_initial_population(GRANA_RADIUS): 
    """Creates a population and returns list of Particle instances"""
    
    print("Creating initial Population.\n")
    global grana_radius
    grana_radius = GRANA_RADIUS
    PSI_PSII_ratio = 1
    
    NPSII = round(1.14*(10**(-3))*np.pi*(grana_radius**2))
    NB6F = round(0.5*NPSII)
    B6F_GRANA = 0.5 #fraction of b6f allocated to grana
    NPSII_grana = round(NPSII-(NB6F*B6F_GRANA))
    NPSII = 1.2*NPSII_grana
    NPSI= round(PSI_PSII_ratio*2*NPSII)
    stroma_radius =  round(math.sqrt((NPSI/(math.pi*2.45*(10**(-3))))+(grana_radius**2)))
    
    #the following are made global for use in analysis
    global grana_area
    global stroma_area
    grana_area = math.pi*(grana_radius**2)
    stroma_area = math.pi*(stroma_radius**2)-grana_area
    
    NC2S2 = round(0.5*NPSII_grana) # % total PSII
    NC2S2M2 = round(0.5*NPSII_grana) # % total PSII
    LHCII_per_RC = 5
    LHCII_grana = 0.7
    NLHCIIg = round(LHCII_grana*(NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4)))
    NB6F_grana = round(B6F_GRANA*NB6F)
    
    
    NB6F_SL = round((1-B6F_GRANA)*NB6F)
    NATP = round(0.7*2*NPSII)
    NPSII_SL = round(0.2*NPSII_grana) #2x10% NPSII or 10% reaction centres
    NLHCII_SL =  round((1-LHCII_grana)*(NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4)))
    
    # Numbers as dictionary
    N_Particles_grana ={"C2S2M2":NC2S2M2,"C2S2":NC2S2,"LHCII":NLHCIIg,"B6F": NB6F_grana}
    N_Particles_SL ={"PSI":NPSI, "LHCII":NLHCII_SL, "B6F": NB6F_SL, "ATP":NATP, "PSII_mono":NPSII_SL}

    # grana boundaries
    
    BOUNDARIES1 = [(round(grana_radius*np.cos(T*(2*np.pi/10000.0))),round(grana_radius*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES2 = [(round((grana_radius+1)*np.cos(T*(2*np.pi/10000.0))),round((grana_radius+1)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    
    BOUNDARIES_grana  = list(set(BOUNDARIES1[:]) | set(BOUNDARIES2[:]) )                                            
    BOUNDARIES_grana = remove_repeat_cols(list_to_matrix(BOUNDARIES_grana))
    
    #Stromal lamellae boundaries
    BOUNDARIES3 = [(round(grana_radius*np.cos(T*(2*np.pi/10000.0))),round(grana_radius*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES4 = [(round((grana_radius-1)*np.cos(T*(2*np.pi/10000.0))),round((grana_radius-1)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES5 = [(round((stroma_radius)*np.cos(T*(2*np.pi/10000.0))),round((stroma_radius)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES6 = [(round((stroma_radius+1)*np.cos(T*(2*np.pi/10000.0))),round((stroma_radius+1)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    
    BOUNDARIES_SL  = list(set(BOUNDARIES3[:]) | set(BOUNDARIES4[:]) | set(BOUNDARIES5[:]) | set(BOUNDARIES6[:])) 
    BOUNDARIES_SL = remove_repeat_cols(list_to_matrix(BOUNDARIES_SL))
    
    #LHCII BOUNDARIES
    BOUNDARIES_LHCII  = list(set(BOUNDARIES5[:]) | set(BOUNDARIES6[:])) 
    BOUNDARIES_LHCII = remove_repeat_cols(list_to_matrix(BOUNDARIES_LHCII))
    global Complexes
    Complexes = pickle.load(open("Complexes.p", 'rb'))
    
    ## create directory for saving intial population
    Dir_Name = "Initial"
    
    import os
    if not os.path.exists(Dir_Name):
        os.makedirs(Dir_Name)
        print("Created directory: "+Dir_Name)
    else:
        print("The directory", Dir_Name,"already exists")
        Continue = "X"
        while Continue not in {"Y", "y", "N", "n"}:
            Continue = input("Do you want to continue with risk of data overwrite? (Y/N) : ")
        if Continue=="N" or Continue =="n":
            print("Aborting script")
            import sys
            sys.exit()
        else:
            print("Continuing with initial population build") # continue with risk of overwrite
    

    
    POPULATION = []# stores objects from every Ptype  
    
    ## Grana particles
    
    for Ptype in ["C2S2M2","C2S2","B6F","LHCII"]:
        print("Adding"+" "+Ptype+": "+str(N_Particles_grana[Ptype])+" particles into grana.")
        NParticles = 0
        while NParticles < N_Particles_grana[Ptype]:
            x =  np.random.randint(0,2*grana_radius)-grana_radius
            y =  np.random.randint(0,2*grana_radius)-grana_radius
            
            while np.sqrt(x**2 + y**2) >= grana_radius:
                x =  np.random.randint(0,2*grana_radius)-grana_radius
                y =  np.random.randint(0,2*grana_radius)-grana_radius
            
            New_Particle = Particle(x,y,theta=0, particle_type=Ptype)
            #coord = [np.matrix([[x],[y]]),0.0, np.random.random()*2*np.pi]                                                   
            
            if len(POPULATION)==0:#only boundaries taken into account
                COLL = collision3(BOUNDARIES_LHCII,New_Particle.Pmatrix)
                sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
                while (COLL == True) or (sqr_dist>= grana_radius**2):
                    coll2, bound, New_Particle_2  = particle_step2(np.matrix([[1000],[2000]]),BOUNDARIES_LHCII,New_Particle) # np.matrix([[1000],[2000]]) is a dummy population
                    sqr_dist_2 = New_Particle_2.location[0,0]**2 + New_Particle_2.location[1,0]**2
                    if sqr_dist_2 < grana_radius**2:
                        New_Particle = New_Particle_2
                        COLL = coll2
                    sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
                    
            
            else:# other particles involved 
                COLL = collision3(np.concatenate((population_matrix(POPULATION),BOUNDARIES_LHCII),axis =1),New_Particle.Pmatrix)
                sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
                while (COLL == True) or (sqr_dist>= grana_radius**2):
                    coll2, bound, New_Particle_2  = particle_step2(population_matrix(POPULATION),BOUNDARIES_LHCII,New_Particle)
                    sqr_dist_2 = New_Particle_2.location[0,0]**2 + New_Particle_2.location[1,0]**2
                    if sqr_dist_2 < grana_radius**2:
                        New_Particle = New_Particle_2
                        COLL = coll2
                    sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
            
            NParticles += 1
            POPULATION.append(New_Particle)
        
    ## Stromal lamellae particles
    
    for Ptype in ["PSI", "B6F", "LHCII",  "ATP", "PSII_mono"]:
        print("Adding"+" "+Ptype+": "+str(N_Particles_SL[Ptype])+" particles into stromal lamellae.")
        NParticles = 0
        while NParticles < N_Particles_SL[Ptype]:
            r = np.random.randint(grana_radius+1,stroma_radius-1)
            theta = np.random.random()*2*np.pi
            x = int(r*np.sin(theta))
            y = int(r*np.cos(theta))
        
            while (np.sqrt(x**2 + y**2) <= grana_radius) or (np.sqrt(x**2 + y**2) >= stroma_radius):
                r = np.random.randint(grana_radius+1,stroma_radius-1)
                theta = np.random.random()*2*np.pi
                x = int(r*np.sin(theta))
                y = int(r*np.cos(theta))
            
            New_Particle = Particle(x,y,theta=0, particle_type=Ptype)
            #coord = [np.matrix([[x],[y]]),0.0, np.random.random()*2*np.pi]                                                   
            sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
            COLL = collision3(np.concatenate((population_matrix(POPULATION),BOUNDARIES_LHCII),axis =1),New_Particle.Pmatrix)
            while COLL == True or (grana_radius**2 > sqr_dist) or (sqr_dist > stroma_radius**2):
                coll2, bound, New_Particle_2 = particle_step2(population_matrix(POPULATION),BOUNDARIES_LHCII,New_Particle)
                sqr_dist_2 = New_Particle_2.location[0,0]**2 + New_Particle_2.location[1,0]**2
                if (grana_radius**2 < sqr_dist_2 < stroma_radius**2):
                    New_Particle = New_Particle_2
                    COLL = coll2
                sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
        
            NParticles += 1
            POPULATION.append(New_Particle)
    
    Save_Population(POPULATION, Dir_Name+"/POPULATION_initial_"+str(grana_radius))
    Plot_Population(POPULATION,Dir_Name+"/POPULATION_initial_plot")
    print("Created initial Population.\n")        
    return POPULATION

def create_test_population(GRANA_RADIUS,Particle_Numbers): 
    """Creates a population and returns list of Particle instances.
    in this version, you can hard code the number of particles of each type by passing a list
    Particle_Numbers= [NC2S2M2,NC2S2,NLHCIIg, NB6F_grana,NPSI,NLHCII_SL,NB6F_SL, NATP, NPSII_SL]"""
    
    print("Creating initial Population.\n")
    global grana_radius
    grana_radius = GRANA_RADIUS
    
    PSI_PSII_ratio = 1
    
    NPSII = round(1.14*(10**(-3))*np.pi*(grana_radius**2))
    NB6F = round(0.5*NPSII)
    B6F_GRANA = 0.5 #fraction of b6f allocated to grana
    NPSII_grana = round(NPSII-(NB6F*B6F_GRANA))
    NPSII = 1.2*NPSII_grana
    NPSI= round(PSI_PSII_ratio*2*NPSII)
    stroma_radius =  round(math.sqrt((NPSI/(math.pi*2.45*(10**(-3))))+(grana_radius**2)))
    
    #the following are made global for use in analysis
    global grana_area
    global stroma_area
    grana_area = math.pi*(grana_radius**2)
    stroma_area = math.pi*(stroma_radius**2)-grana_area
    
    NC2S2 = round(0.5*NPSII_grana) # % total PSII
    NC2S2M2 = round(0.5*NPSII_grana) # % total PSII
    LHCII_per_RC = 5
    LHCII_grana = 0.7
    NLHCIIg = round(LHCII_grana*(NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4)))
    NB6F_grana = round(B6F_GRANA*NB6F)
    
    
    NB6F_SL = round((1-B6F_GRANA)*NB6F)
    NATP = round(0.7*2*NPSII)
    NPSII_SL = round(0.2*NPSII_grana) #2x10% NPSII or 10% reaction centres
    NLHCII_SL =  round((1-LHCII_grana)*(NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4)))
    
    
    [NC2S2M2,NC2S2,NLHCIIg, NB6F_grana,NPSI,NLHCII_SL,NB6F_SL, NATP, NPSII_SL] = Particle_Numbers
    
    # Numbers as dictionary
    N_Particles_grana ={"C2S2M2":NC2S2M2,"C2S2":NC2S2,"LHCII":NLHCIIg,"B6F": NB6F_grana}
    N_Particles_SL ={"PSI":NPSI, "LHCII":NLHCII_SL, "B6F": NB6F_SL, "ATP":NATP, "PSII_mono":NPSII_SL}

    # grana boundaries
    
    BOUNDARIES1 = [(round(grana_radius*np.cos(T*(2*np.pi/10000.0))),round(grana_radius*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES2 = [(round((grana_radius+1)*np.cos(T*(2*np.pi/10000.0))),round((grana_radius+1)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    
    BOUNDARIES_grana  = list(set(BOUNDARIES1[:]) | set(BOUNDARIES2[:]) )                                            
    BOUNDARIES_grana = remove_repeat_cols(list_to_matrix(BOUNDARIES_grana))
    
    #Stromal lamellae boundaries
    BOUNDARIES3 = [(round(grana_radius*np.cos(T*(2*np.pi/10000.0))),round(grana_radius*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES4 = [(round((grana_radius-1)*np.cos(T*(2*np.pi/10000.0))),round((grana_radius-1)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES5 = [(round((stroma_radius)*np.cos(T*(2*np.pi/10000.0))),round((stroma_radius)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES6 = [(round((stroma_radius+1)*np.cos(T*(2*np.pi/10000.0))),round((stroma_radius+1)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    
    BOUNDARIES_SL  = list(set(BOUNDARIES3[:]) | set(BOUNDARIES4[:]) | set(BOUNDARIES5[:]) | set(BOUNDARIES6[:])) 
    BOUNDARIES_SL = remove_repeat_cols(list_to_matrix(BOUNDARIES_SL))
    
    #LHCII BOUNDARIES
    BOUNDARIES_LHCII  = list(set(BOUNDARIES5[:]) | set(BOUNDARIES6[:])) 
    BOUNDARIES_LHCII = remove_repeat_cols(list_to_matrix(BOUNDARIES_LHCII))
    global Complexes
    Complexes = pickle.load(open("Complexes.p", 'rb'))
    
    ## create directory for saving intial population
    Dir_Name = "Initial_test"
    
    import os
    if not os.path.exists(Dir_Name):
        os.makedirs(Dir_Name)
        print("Created directory: "+Dir_Name)
    else:
        print("The directory", Dir_Name,"already exists")
        Continue = "X"
        while Continue not in {"Y", "y", "N", "n"}:
            Continue = input("Do you want to continue with risk of data overwrite? (Y/N) : ")
        if Continue=="N" or Continue =="n":
            print("Aborting script")
            import sys
            sys.exit()
        else:
            print("Continuing with initial population build") # continue with risk of overwrite
    

    
    POPULATION = []# stores objects from every Ptype  
    
    ## Grana particles
    
    for Ptype in ["C2S2M2","C2S2","B6F","LHCII"]:
        print("Adding"+" "+Ptype+": "+str(N_Particles_grana[Ptype])+" particles into grana.")
        NParticles = 0
        while NParticles < N_Particles_grana[Ptype]:
            x =  np.random.randint(0,2*grana_radius)-grana_radius
            y =  np.random.randint(0,2*grana_radius)-grana_radius
            
            while np.sqrt(x**2 + y**2) >= grana_radius:
                x =  np.random.randint(0,2*grana_radius)-grana_radius
                y =  np.random.randint(0,2*grana_radius)-grana_radius
            
            New_Particle = Particle(x,y,theta=0, particle_type=Ptype)
            #coord = [np.matrix([[x],[y]]),0.0, np.random.random()*2*np.pi]                                                   
            
            if len(POPULATION)==0:#only boundaries taken into account
                COLL = collision3(BOUNDARIES_LHCII,New_Particle.Pmatrix)
                sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
                while (COLL == True) or (sqr_dist>= grana_radius**2):
                    coll2, bound, New_Particle_2  = particle_step2(np.matrix([[1000],[2000]]),BOUNDARIES_LHCII,New_Particle) # np.matrix([[1000],[2000]]) is a dummy population
                    sqr_dist_2 = New_Particle_2.location[0,0]**2 + New_Particle_2.location[1,0]**2
                    if sqr_dist_2 < grana_radius**2:
                        New_Particle = New_Particle_2
                        COLL = coll2
                    sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
                    
            
            else:# other particles involved 
                COLL = collision3(np.concatenate((population_matrix(POPULATION),BOUNDARIES_LHCII),axis =1),New_Particle.Pmatrix)
                sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
                while (COLL == True) or (sqr_dist>= grana_radius**2):
                    coll2, bound, New_Particle_2  = particle_step2(population_matrix(POPULATION),BOUNDARIES_LHCII,New_Particle)
                    sqr_dist_2 = New_Particle_2.location[0,0]**2 + New_Particle_2.location[1,0]**2
                    if sqr_dist_2 < grana_radius**2:
                        New_Particle = New_Particle_2
                        COLL = coll2
                    sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
            
            NParticles += 1
            POPULATION.append(New_Particle)
        
    ## Stromal lamellae particles
    
    for Ptype in ["PSI", "B6F", "LHCII",  "ATP", "PSII_mono"]:
        print("Adding"+" "+Ptype+": "+str(N_Particles_SL[Ptype])+" particles into stromal lamellae.")
        NParticles = 0
        while NParticles < N_Particles_SL[Ptype]:
            r = np.random.randint(grana_radius+1,stroma_radius-1)
            theta = np.random.random()*2*np.pi
            x = int(r*np.sin(theta))
            y = int(r*np.cos(theta))
        
            while (np.sqrt(x**2 + y**2) <= grana_radius) or (np.sqrt(x**2 + y**2) >= stroma_radius):
                r = np.random.randint(grana_radius+1,stroma_radius-1)
                theta = np.random.random()*2*np.pi
                x = int(r*np.sin(theta))
                y = int(r*np.cos(theta))
            
            New_Particle = Particle(x,y,theta=0, particle_type=Ptype)
            #coord = [np.matrix([[x],[y]]),0.0, np.random.random()*2*np.pi]                                                   
            sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
            COLL = collision3(np.concatenate((population_matrix(POPULATION),BOUNDARIES_LHCII),axis =1),New_Particle.Pmatrix)
            while COLL == True or (grana_radius**2 > sqr_dist) or (sqr_dist > stroma_radius**2):
                coll2, bound, New_Particle_2 = particle_step2(population_matrix(POPULATION),BOUNDARIES_LHCII,New_Particle)
                sqr_dist_2 = New_Particle_2.location[0,0]**2 + New_Particle_2.location[1,0]**2
                if (grana_radius**2 < sqr_dist_2 < stroma_radius**2):
                    New_Particle = New_Particle_2
                    COLL = coll2
                sqr_dist = New_Particle.location[0,0]**2 + New_Particle.location[1,0]**2
        
            NParticles += 1
            POPULATION.append(New_Particle)
    
    Save_Population(POPULATION, Dir_Name+"/POPULATION_initial_test_"+str(grana_radius))
    Plot_Population(POPULATION,Dir_Name+"/POPULATION_initial_test_plot")
    print("Created initial Population.\n")        
    return POPULATION


def particle_step2(POPULAT,BOUNDARIES,PARTICLE):
    
    
    """in this version of the function, POPULATION is a matrix not a list"""
    dxdy = np.matrix(random.choice([[0,1],[0,-1],[1,0],[-1,0]])).T                 
    d0 = random.choice([np.pi/33.0,-1*np.pi/33.0])

    
    PARTICLE2 = SudoParticle(PARTICLE.location+dxdy,PARTICLE.rotation+d0,PARTICLE.Ptype)   
    #coord = [x0+np.matrix(random.choice([[0,1,0],[0,-1,0],[1,0,0],[-1,0,0],[0,0,1],[0,0,-1]])).T,(th0+random.choice([np.pi/12,-1*np.pi/12]))%2*np.pi,(th0+random.choice([np.pi/12,-1*np.pi/12]))%2*np.pi ] 
    COLL = collision3(np.concatenate((POPULAT,BOUNDARIES),axis =1),PARTICLE2.Pmatrix)
    BOUND = boundary_check(PARTICLE.Ptype,PARTICLE.location[0,0],dxdy[0,0],PARTICLE.location[1,0],dxdy[1,0])
    return COLL, BOUND, PARTICLE2  


def time_step(Population1,Population2,L,BOUNDARIES,layer,current_energy):
    """chooses a random particle from population.Population is a list of particles. L is length of population"""
    if layer == 1:
        Population = Population1
        other_population = Population2
    else:
        Population = Population2
        other_population = Population1
    x = np.random.randint(0,L) # index of particle
    
    coll, bound, new_particle = particle_step2(population_matrix(Population[0:x]+Population[x+1:]),BOUNDARIES,Population[x])
    coll = (bound or coll)
    # coll = collision?, bound = gone from one domain to another?
    if coll == False:
        new_energy = Hamiltonian(Population[0:x]+[new_particle]+Population[x+1:],other_population)

        dE =  current_energy - new_energy # change in energy

        
        
        # Metropolis algorithm
        if dE <= 0:
            Population[x] = new_particle
            energy = new_energy
        else:
            if random.random() <= math.e**(-dE):
                Population[x] = new_particle
                energy = new_energy
            else:
                energy = current_energy
    else:
        energy = current_energy
        
    return Population, energy


def Alike_particles(Population, Ptype):
    return [p for p in Population if p.Ptype == Ptype]


# def NN_vals(pop):
#     T = [np.array(x[0]) for x in pop]
#     for t1 in T:
#         vals = []
#         for t2 in T:
#             if list(t1)!=list(t2):
#                 vals.append(np.sqrt((t1[0][0]-t2[0][0])**2+(t1[1][0]-t2[1][0])**2))
#         print(min(vals))                  


def load_coordinates(filename, ptype):
    """loads from the files: line: x y theta"""
    coords = open(filename,'r').readlines()
    POP = []
    for coord in coords:
        f=coord.split()
        POP.append(Particle(float(f[0]),float(f[1]),float(f[2]),ptype))
    
    return POP 


def distances(M1,M2):
    """returns a 1D matrix of distances between M1 and M2"""
    D = dist.cdist(M1.T,M2.T)
    return np.absolute(np.ndarray.flatten(D))


def intra_layer(D):
    """calculates intralayer potentials from 1D matrix of distances"""
    return E_intra*(D[(6<=D)&(D<=8)].shape[0])/2.0 # division by 2 accounts for repeat measures 


def inter_layer(X1,X2):
    """X1 and X2 are the LHCII locations in layers 1 and 2"""
    D = distances(X1,X2)
    return np.sum(E_inter*np.square((D[D<LHCII_radius]-LHCII_radius)/LHCII_radius))/2.0 # division by 2 accounts for repeat measures 


def boundary_check(Ptype,x,dx,y,dy):
    """this function returns true if a particle tries to move from grana to SL or vice versa"""
    #NON_MOVING = {"C2S2M2","C2S2","B6Fg","B6Fs","PSI","ATP","PSIIs"} # the particle types which are restricted to their respective regions
    if Ptype != "LHCII": #in NON_MOVING:
        if ((x**2+y**2 <= grana_radius**2) and ((x+dx)**2+(y+dy)**2 >grana_radius**2)) or ((x**2+y**2 > grana_radius**2) and ((x+dx)**2+(y+dy)**2 <= grana_radius**2)):
            return True # the movement is ount of the region
        else:
            return False
    else:
        return False


def PSII_Hamiltonian(Population1,Population2):
    LHCIIs_Layer1_tuple = tuple(p.bound_LHCII() for p in Population1 if (p.Ptype=="C2S2M2") or (p.Ptype=="C2S2") or (p.Ptype=="LHCII"))
    LHCIIs_Layer2_tuple = tuple(p.bound_LHCII() for p in Population2 if (p.Ptype=="C2S2M2") or (p.Ptype=="C2S2") or (p.Ptype=="LHCII"))
    
    
    LHCIIs_Layer1 = np.concatenate(LHCIIs_Layer1_tuple, axis=1)
    LHCIIs_Layer2 = np.concatenate(LHCIIs_Layer2_tuple, axis=1)
    
    # Distances between LHCIIs in each layer
    D1 = distances(LHCIIs_Layer1,LHCIIs_Layer1)
    D2 = distances(LHCIIs_Layer2,LHCIIs_Layer2)
    
    # intra particle interaction energy in each layer
    Eintra1 = intra_layer(D1)
    Eintra2 = intra_layer(D2)
    
    # interlayer interaction energy (stacking interaction)
    
    # Stacking interactions only happen within grana radius
    LHCIIs_Layer1_tuple_stacking = tuple([p for p in LHCIIs_Layer1_tuple if np.multiply(p,p).sum() < grana_radius**2])
    LHCIIs_Layer2_tuple_stacking = tuple([p for p in LHCIIs_Layer2_tuple if np.multiply(p,p).sum() < grana_radius**2])
     
    LHCIIs_Layer1 = np.concatenate(LHCIIs_Layer1_tuple_stacking, axis=1)
    LHCIIs_Layer2 = np.concatenate(LHCIIs_Layer2_tuple_stacking, axis=1)
    
    Einter1 = inter_layer(LHCIIs_Layer1,LHCIIs_Layer2)
    
    return Eintra1 + Eintra2 + Einter1


## functions for PSI Hamiltonian ##

def intra_layer_PSI_LHCII(D):
    return E_intra_PSI*(D[D<=1].shape[0])


def PSI_Hamiltonian(Population1,Population2):
    LHCIIs_Layer1 = np.concatenate(tuple(p.bound_LHCII() for p in Population1 if (p.Ptype=="LHCII")), axis=1)
    LHCIIs_Layer2 = np.concatenate(tuple(p.bound_LHCII() for p in Population2 if (p.Ptype=="LHCII")), axis=1)
    
    PSIs_Layer1 =  np.concatenate(tuple(p.PSI_site() for p in Population1 if (p.Ptype=="PSI")), axis=1)
    PSIs_Layer2 =  np.concatenate(tuple(p.PSI_site() for p in Population2 if (p.Ptype=="PSI")), axis=1)
    
    # Distances between LHCIIs in each layer
    D1 = distances(LHCIIs_Layer1,PSIs_Layer1)
    D2 = distances(LHCIIs_Layer2,PSIs_Layer2)
    
    # intra particle interaction energy in each layer
    Eintra1 = intra_layer_PSI_LHCII(D1)
    Eintra2 = intra_layer_PSI_LHCII(D2)
    
    return Eintra1 + Eintra2

def Hamiltonian(Population1,Population2):

    H = 0 # total Energy
    # only calculate energy where interactions exist
    if (E_intra != 0) or (E_inter != 0):
        H += PSII_Hamiltonian(Population1,Population2)
    if E_intra_PSI != 0:
        H += PSI_Hamiltonian(Population1,Population2)
    return H 

def Run_Simulation(GRANA_RADIUS=170,DATE="11_05_18",EXPERIMENT="FREE_LHCII",TMAX=1100001,Stacking_Interaction_Energy = 4, LHCII_Binding_Interaction_Energy = 2,PSI_Interaction_Energy=0,POPULATION1=[],POPULATION2=[]):
    
    ##parameters which describe the model
    ##see Wood et al., 2019/2020 for more details
    
    
    global grana_radius
    grana_radius = GRANA_RADIUS
    PSI_PSII_ratio = 1
    
    NPSII = round(1.14*(10**(-3))*np.pi*(grana_radius**2))
    NB6F = round(0.5*NPSII)
    B6F_GRANA = 0.5 #fraction of b6f allocated to grana
    NPSII_grana = round(NPSII-(NB6F*B6F_GRANA))
    NPSII = 1.2*NPSII_grana
    NPSI= round(PSI_PSII_ratio*2*NPSII)
    stroma_radius =  round(math.sqrt((NPSI/(math.pi*2.45*(10**(-3))))+(grana_radius**2)))
    
    #the following are made global for use in analysis
    global grana_area
    global stroma_area
    grana_area = math.pi*(grana_radius**2)
    stroma_area = math.pi*(stroma_radius**2)-grana_area
    
    NC2S2 = round(0.5*NPSII_grana) # % total PSII
    NC2S2M2 = round(0.5*NPSII_grana) # % total PSII
    LHCII_per_RC = 5
    LHCII_grana = 0.7
    NLHCIIg = round(LHCII_grana*(NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4)))
    NB6F_grana = round(B6F_GRANA*NB6F)
    
    
    NB6F_SL = round((1-B6F_GRANA)*NB6F)
    NATP = round(0.7*2*NPSII)
    NPSII_SL = round(0.2*NPSII_grana) #2x10% NPSII or 10% reaction centres
    NLHCII_SL =  round((1-LHCII_grana)*(NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4)))
    
    
    #Probabilities of moves in 1us (for diffusion)
    #translations
    Pt_C2S2M2 = 0.59
    Pt_C2S2 = 0.63
    Pt_B6F = 0.71
    Pt_LHCII = 0.77
    Pt_PSIImon = 0.71
    Pt_PSI = 0.67
    Pt_ATP = 0.77
    
    #rotations
    Pr_C2S2M2 = 0.06
    Pr_C2S2 = 0.09
    Pr_B6F = 0.33
    Pr_LHCII = 0.43
    Pr_PSIImon = 0.26
    Pr_ATP = 0.67
    Pr_PSI = 0.17
    
 
    # energy parameters
    global LHCII_radius
    global interaction_radius
    global E_intra
    global E_inter
    global E_intra_PSI
    
    LHCII_radius = 3.5 # nm
    interaction_radius = 5.5 # nm
    E_intra =  LHCII_Binding_Interaction_Energy#KT , intralayer interactions
    E_inter = Stacking_Interaction_Energy #KT, interactions over the stromal gap
    E_intra_PSI = PSI_Interaction_Energy #KT. binding o LHCII to PSI
    # grana boundaries
    
    BOUNDARIES1 = [(round(grana_radius*np.cos(T*(2*np.pi/10000.0))),round(grana_radius*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES2 = [(round((grana_radius+1)*np.cos(T*(2*np.pi/10000.0))),round((grana_radius+1)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    
    BOUNDARIES_grana  = list(set(BOUNDARIES1[:]) | set(BOUNDARIES2[:]) )                                            
    BOUNDARIES_grana = remove_repeat_cols(list_to_matrix(BOUNDARIES_grana))
    
    #Stromal lamellae boundaries
    BOUNDARIES3 = [(round(grana_radius*np.cos(T*(2*np.pi/10000.0))),round(grana_radius*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES4 = [(round((grana_radius-1)*np.cos(T*(2*np.pi/10000.0))),round((grana_radius-1)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES5 = [(round((stroma_radius)*np.cos(T*(2*np.pi/10000.0))),round((stroma_radius)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    BOUNDARIES6 = [(round((stroma_radius+1)*np.cos(T*(2*np.pi/10000.0))),round((stroma_radius+1)*np.sin(T*(2*np.pi/10000.0)))) for T in range(10000)]
    
    BOUNDARIES_SL  = list(set(BOUNDARIES3[:]) | set(BOUNDARIES4[:]) | set(BOUNDARIES5[:]) | set(BOUNDARIES6[:])) 
    BOUNDARIES_SL = remove_repeat_cols(list_to_matrix(BOUNDARIES_SL))
    
    #LHCII BOUNDARIES
    BOUNDARIES_LHCII  = list(set(BOUNDARIES5[:]) | set(BOUNDARIES6[:])) 
    BOUNDARIES_LHCII = remove_repeat_cols(list_to_matrix(BOUNDARIES_LHCII))
    global Complexes
    Complexes = pickle.load(open("Complexes.p", 'rb'))

    
    # loading in populations
    #POPULATION1 = pickle.load(open("Population1_initial", 'rb'))
    #POPULATION2 = pickle.load(open("Population2_initial", 'rb'))
    
    #grana size = 170
    if POPULATION1==[]:
        if grana_radius == 170:
            POPULATION1 = Load_Population("Initial/POPULATION1_random_170")
        elif grana_radius == 190:#grana size = 190
            POPULATION1 = Load_Population("Initial/POPULATION1_random_190")
        else:
            print("Incorrect grana size")
            print("Aborting script")
            import sys
            sys.exit()
    else:
        print("Using loaded POPULATION1...")
    
    if POPULATION2==[]:
        if grana_radius == 170:
            POPULATION2 = Load_Population("Initial/POPULATION2_random_170")
        elif grana_radius == 190:#grana size = 190
            POPULATION2 = Load_Population("Initial/POPULATION2_random_190")
        else:
            print("Incorrect grana size")
            print("Aborting script")
            import sys
            sys.exit()
    else:
        print("Using loaded POPULATION2...")
    
    
    print("Model Initialised. Running Simulation\n")
    # Create the directory for the experiment
    Dir_Name = EXPERIMENT+"_"+DATE
    
    import os
    if not os.path.exists(Dir_Name):
        os.makedirs(Dir_Name)
        print("Created directory: "+Dir_Name)
    else:
        print("The directory", Dir_Name,"already exists")
        Continue = "X"
        while Continue not in {"Y", "y", "N", "n"}:
            Continue = input("Do you want to continue with risk of data overwrite? (Y/N) : ")
        if Continue=="N" or Continue =="n":
            print("Aborting script")
            import sys
            sys.exit()
        else:
            print("Continuing with simulation") # continue with risk of overwrite
    
    # A directory for the data itself (Pickle files of the populations)
    if not os.path.exists(Dir_Name+"/Data"):
        os.makedirs(Dir_Name+"/Data")
    else:
        pass
    
    
    Nparticles = len(POPULATION1)

    
    Tmax= TMAX
    ENERGY = [Hamiltonian(POPULATION1,POPULATION2)]
    T_energy = [0]
    
    X_array = np.random.random(Tmax)
    for t in np.arange(Tmax):
        if (t%1000000 == 0) and (t<10000000):
            print("Equilibration "+str(int(t*(10/1000000)))+"% complete")
            print("Energy = ", str(ENERGY[-1]))
            print(" LHCII in grana = ",LHCII_Localisation_Analysis(POPULATION1),"/",LHCII_Localisation_Analysis(POPULATION2))
            print(" LHCII NN = ",np.mean(Nearest_Neighbour_Distances(POPULATION1, ['LHCII'])),"/",np.mean(Nearest_Neighbour_Distances(POPULATION2, ['LHCII'])))
        elif t==10000000:
            print("Collecting data")
            
    
        
        # the following is a decision tree that leads to the perturbation
        X = X_array[t]
        if X<0.5: # Layer one
            
            layer = 1
            POPULATION1, new_energy = time_step(POPULATION1,POPULATION2, Nparticles, BOUNDARIES_LHCII, layer, ENERGY[-1])
            ENERGY.append(new_energy)
            T_energy.append(t)
        else: # Layer two
            layer = 2
            POPULATION2 , new_energy = time_step(POPULATION1, POPULATION2,Nparticles,BOUNDARIES_LHCII, layer, ENERGY[-1])
            ENERGY.append(new_energy)
            T_energy.append(t)
    
        
        # Checkpoint saves are there to start again if code is interupted
        # after 10M iterations data is collected every 1000 iterations
        if (t%1000000==0) and (10000000>t): # saves every millionth iteration
            Save_Population(POPULATION1, Dir_Name+"/Data/POPULATION1_Checkpoint_"+str(t))
            Save_Population(POPULATION2, Dir_Name+"/Data/POPULATION2_Checkpoint_"+str(t))
            
        elif (t>=10000000) and (t%1000==0):
            Save_Population(POPULATION1, Dir_Name+"/Data/POPULATION1_"+str(t))
            Save_Population(POPULATION2, Dir_Name+"/Data/POPULATION2_"+str(t))
            
            # Here I want to add a single file with all results in
            
                
    ENERGY_DF = pd.DataFrame({"Iters" : T_energy, "Energy":ENERGY})
    ENERGY_DF.to_csv(Dir_Name+"/Equilibration_Energy.csv")
    Plot_Population(POPULATION1,Dir_Name+"/POPULATION1_Post_Simulation_plot")
    Plot_Population(POPULATION2,Dir_Name+"/POPULATION2_Post_Simulation_plot")
    return POPULATION1, POPULATION2


def Nearest_Neighbour_Distances(POPULATION, ptypes):
    """returns centre distances of the nearest neighbours of particles in POPULATION.
    ptypes is a list or set of Ptypes to be considered E.G 'LHCII'. """
    Particles = [p.location.T.tolist()[0] for p in POPULATION if p.Ptype in ptypes]
    Nearest_Neighbours = []
    for P1 in Particles:
        Sqr_Distances = []
        for P2 in Particles:
            if P1 != P2:
                Sqr_Distances.append((P1[0]-P2[0])**2 + (P1[1]-P2[1])**2)
        Nearest_Neighbours.append(np.sqrt(np.min(Sqr_Distances)))
    return Nearest_Neighbours


def LHCII_Localisation_Analysis(POPULATION):
    """Returns the fraction of LHCII in the grana of POPULATION"""
    LHCII_locations = [p.location.T.tolist()[0] for p in POPULATION if p.Ptype == "LHCII"]
    
    NLHCII_grana = len([i for i in LHCII_locations if (i[0]**2 + i[1]**2 <= grana_radius**2)])
    NLHCII_sl = len([j for j in LHCII_locations if (j[0]**2 + j[1]**2 > grana_radius**2)])
        
    
    Fraction_LHCII_in_grana = float(NLHCII_grana)/(NLHCII_sl+NLHCII_grana)
    return Fraction_LHCII_in_grana

def Density_Analysis(POPULATION):
    """Returns the area fraction of protein in the grana and the stromal lamellae"""
    
    # lower case variables are the number of lattice sites taken by the particles
    nc2s2m2 = NC2S2M2*(Complexes["C2S2M2"].shape[1])
    nc2s2 = NC2S2M2*(Complexes["C2S2M2"].shape[1])
    nlhcii = (NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4))*(Complexes["LHCII"].shape[1])
    npsi = NPSI*(Complexes["PSI"].shape[1])
    nb6f_g = NB6F_grana*(Complexes["B6F"].shape[1])
    nb6f_sl = NB6F_SL*(Complexes["B6F"].shape[1])
    natp = NATP*(Complexes["ATP"].shape[1])
    npsii_sl = NPSII_SL*(Complexes["PSII_mono"].shape[1])
    lhcii_grana_fraction = LHCII_Localisation_Analysis(POPULATION)
    
    Ngrana = nc2s2m2 + nc2s2 + lhcii_grana_fraction*nlhcii + nb6f_g
    Nsl = npsi + natp + nb6f_sl + npsii_sl + (1-lhcii_grana_fraction)*nlhcii
    grana_density = float(Ngrana)/grana_area 
    sl_density = float(Nsl)/stroma_area
    
    return grana_density, sl_density




class Dgraph:
    
    
    xmax_factor = 10000 # this is used in init as the factor c (c*x+y)
    def __init__(self,coords,threshold):
        """creates the nx graph if the distance is less than or equal to threshold"""
        
        # import relevant modules
        import numpy as np
        import networkx as nx
        import matplotlib.pyplot as plt
        import scipy.spatial.distance as dist

        
        self.G = nx.Graph()
        self.coordinates = coords
        #node_indices apply a unique number to the x,y coordinates
        """
        for xy1 in coords:
            for xy2 in coords:
                if xy1 != xy2:
                    if self.distance(xy1,xy2) <= threshold:
                        self.G.add_edge(self.node_num(xy1),self.node_num(xy2),weight=self.distance(xy1,xy2))
        """
        D = dist.cdist(np.array(coords),np.array(coords))

        Y = np.arange(len(coords))
        for x1 in Y:
            for x2 in Y:
                if x1 != x2:
                    if D[x1,x2] <= threshold:
                        self.G.add_edge(self.node_num(coords[x1]),self.node_num(coords[x2]),weight=D[x1,x2])
        
    def node_num(self, X):
        """transforms an X = [x,y] into c*x+y using xmax_factor"""
        min = 500
        xmax_factor = 10000
        return (X[0]+500)*xmax_factor+(X[1]+500)
    
    @staticmethod
    def distance(x1,x2):
        """calculates the Euclidean distance between x1 and x2 (in the form [x,y])"""
        
        return float(dist.cdist(np.matrix(x1), np.matrix(x2)))
    
    def draw(self,nodecolor='b', edge_color='k'):
        
        options = {
            'node_color': 'blue',
            'node_size': 25,
            'line_color': 'grey',
            'linewidths': 0,
            'width': 0.1,
        }
        nx.draw(self.G, **options)

    
    def Ddraw(self):
        """draws the graph using coordinates"""
        
        X = []
        Y = []
        # draw the edges
        for x in self.coordinates:
            for x2 in self.coordinates:
                if [self.node_num(x), self.node_num(x2)] in self.G.edges: 
                    X.append([x[0],x2[0]])
                    Y.append([x[1],x2[1]])
        
        plt.plot(np.array(X).T,np.array(Y).T,'k-', lw=1)
        
        #draw the nodes
        #XY_array = np.array(self.coordinates).T
        #plt.scatter(XY_array[0,:],XY_array[1,:],color='k')
        
        plt.show()
    
    def edges(self):
        return self.G.edges
    
    def has_path(self,x1,x2):
        if [self.node_num(x1),self.node_num(x2)] in self.G.edges():
            return nx.has_path(self.G, self.node_num(x1),self.node_num(x2))
        else:
            return False
        
 
# Network classes for antenna analysis 
class Chlorophyll_Network(Dgraph):
        
    def __init__(self,POPULATION,threshold=1):
        #Node_list = coordinate_list(C2S2M2_POP+C2S2_POP+LHCII_POP)
        # import relevant modules
        import numpy as np
        import networkx as nx
        import matplotlib.pyplot as plt

        import scipy.spatial.distance as dist
        
        #instantiate the underlying nx graph
        
        
        self.G = nx.Graph()
        
        self.Distance_Cutoff = 50 # Distance an exciton can travel before fluorescence
        
        
        
        #C2S2M2
        C2S2M2_POP = [p for p in POPULATION if p.Ptype =="C2S2M2"]
        
        #C2S2
        C2S2_POP = [p for p in POPULATION if p.Ptype =="C2S2"]
        
        #LHCII
        LHCII_POP = [p for p in POPULATION if p.Ptype =="LHCII"]
        
        self.PSII_coordinates = [p.location.T.tolist()[0] for p in C2S2M2_POP+C2S2_POP]
        self.coordinates = [p.location.T.tolist()[0] for p in C2S2M2_POP+C2S2_POP+LHCII_POP]
        
        
        #go through the combinations of C2S2M2 and other complexes
        for c2s2m2 in C2S2M2_POP:
            
            C2S2M2_Chlorophylls = c2s2m2.bound_Chorophylls()
            for c2s2m2_2 in C2S2M2_POP:
                
                if c2s2m2.location.T.tolist()[0] != c2s2m2_2.location.T.tolist()[0]:
                    C2S2M2_Chlorophylls_2 = c2s2m2_2.bound_Chorophylls()
                    d = np.min(self.distances(C2S2M2_Chlorophylls,C2S2M2_Chlorophylls_2))
                    if d <=threshold:
                        
                        self.G.add_edge(self.node_num(c2s2m2.location.T.tolist()[0]),self.node_num(c2s2m2_2.location.T.tolist()[0]),weight=d)
                
            for c2s2 in C2S2_POP:
                C2S2_Chlorophylls = c2s2.bound_Chorophylls()

                d = np.min(self.distances(C2S2M2_Chlorophylls,C2S2_Chlorophylls))
                if d <=threshold:
                    self.G.add_edge(self.node_num(c2s2m2.location.T.tolist()[0]),self.node_num(c2s2.location.T.tolist()[0]),weight=d)

            for lhcii in LHCII_POP:
                LHCII_Chlorophylls = lhcii.bound_Chorophylls()

                d = np.min(self.distances(C2S2M2_Chlorophylls,LHCII_Chlorophylls))
                if d <=threshold:
                    self.G.add_edge(self.node_num(c2s2m2.location.T.tolist()[0]),self.node_num(lhcii.location.T.tolist()[0]),weight=d)
        
        #C2S2
        for c2s2 in C2S2_POP:
            C2S2_Chlorophylls = c2s2.bound_Chorophylls()
                
            for c2s2m2 in C2S2M2_POP:
                    C2S2M2_Chlorophylls = c2s2m2.bound_Chorophylls()
                    d = np.min(self.distances(C2S2_Chlorophylls,C2S2M2_Chlorophylls))
                    if d <=threshold:
                        self.G.add_edge(self.node_num(c2s2.location.T.tolist()[0]),self.node_num(c2s2m2.location.T.tolist()[0]),weight=d)

            for c2s2_2 in C2S2_POP:
                if c2s2.location.T.tolist()[0] != c2s2_2.location.T.tolist()[0]:
                    C2S2_Chlorophylls_2 = c2s2_2.bound_Chorophylls()

                    d = np.min(self.distances(C2S2_Chlorophylls,C2S2_Chlorophylls_2))
                    if d <=threshold:
                        self.G.add_edge(self.node_num(c2s2.location.T.tolist()[0]),self.node_num(c2s2_2.location.T.tolist()[0]),weight=d)

            for lhcii in LHCII_POP:
                LHCII_Chlorophylls = lhcii.bound_Chorophylls()

                d = np.min(self.distances(C2S2_Chlorophylls,LHCII_Chlorophylls))
                if d <=threshold:
                    self.G.add_edge(self.node_num(c2s2.location.T.tolist()[0]),self.node_num(lhcii.location.T.tolist()[0]),weight=d)

        #LHCII
        for lhcii in LHCII_POP:
            LHCII_Chlorophylls = lhcii.bound_Chorophylls()
            for lhcii_2 in LHCII_POP:
                if lhcii.location.T.tolist()[0] != lhcii_2.location.T.tolist()[0]:
                    LHCII_Chlorophylls_2 = lhcii_2.bound_Chorophylls()

                    d = np.min(self.distances(LHCII_Chlorophylls,LHCII_Chlorophylls_2))
                    if d <=threshold:
                        self.G.add_edge(self.node_num(lhcii.location.T.tolist()[0]),self.node_num(lhcii_2.location.T.tolist()[0]),weight=d)
            
            for c2s2m2 in C2S2M2_POP:
                    C2S2M2_Chlorophylls = c2s2m2.bound_Chorophylls()
                    d = np.min(self.distances(LHCII_Chlorophylls,C2S2M2_Chlorophylls))
                    if d <=threshold:
                        self.G.add_edge(self.node_num(lhcii.location.T.tolist()[0]),self.node_num(c2s2m2.location.T.tolist()[0]),weight=d)

            for c2s2 in C2S2_POP:
                C2S2_Chlorophylls = c2s2.bound_Chorophylls()

                d = np.min(self.distances(LHCII_Chlorophylls,C2S2_Chlorophylls))
                if d <=threshold:
                    self.G.add_edge(self.node_num(lhcii.location.T.tolist()[0]),self.node_num(c2s2.location.T.tolist()[0]),weight=d)
    
        
    @staticmethod
    def distances(M1,M2):
        """returns a 1D matrix of distances between M1 and M2"""
        D = dist.cdist(M1.T,M2.T)
        return np.absolute(np.ndarray.flatten(D))


# Antenna network for PSI
class PSI_Chlorophyll_Network(Dgraph):
        
    def __init__(self,POPULATION,threshold=1):
        #Node_list = coordinate_list(C2S2M2_POP+C2S2_POP+LHCII_POP)
        # import relevant modules
        import numpy as np
        import networkx as nx
        import matplotlib.pyplot as plt

        import scipy.spatial.distance as dist
        
        #instantiate the underlying nx graph
        
        
        self.G = nx.Graph()
        
        self.Distance_Cutoff = 50 # Distance an exciton can travel before fluorescence
        
        
        
        #PSI
        PSI_POP = [p for p in POPULATION if p.Ptype =="PSI"]
        
        #LHCII
        LHCII_POP = [p for p in POPULATION if p.Ptype =="LHCII"]
        
        self.PSI_coordinates = [p.location.T.tolist()[0] for p in PSI_POP]

        
        
        #go through the combinations of C2S2M2 and other complexes
        for psi in PSI_POP:
            
            PSI_Chlorophylls = psi.bound_Chorophylls()
 
            for lhcii in LHCII_POP:
                LHCII_Chlorophylls = lhcii.bound_Chorophylls()

                d = np.min(self.distances(PSI_Chlorophylls,LHCII_Chlorophylls))
                if d <=threshold:
                    self.G.add_edge(self.node_num(psi.location.T.tolist()[0]),self.node_num(lhcii.location.T.tolist()[0]),weight=d)

        #LHCII
        for lhcii in LHCII_POP:
            LHCII_Chlorophylls = lhcii.bound_Chorophylls()
            for lhcii_2 in LHCII_POP:
                if lhcii.location.T.tolist()[0] != lhcii_2.location.T.tolist()[0]:
                    LHCII_Chlorophylls_2 = lhcii_2.bound_Chorophylls()

                    d = np.min(self.distances(LHCII_Chlorophylls,LHCII_Chlorophylls_2))
                    if d <=threshold:
                        self.G.add_edge(self.node_num(lhcii.location.T.tolist()[0]),self.node_num(lhcii_2.location.T.tolist()[0]),weight=d)
            
            for psi in PSI_POP:
            
                PSI_Chlorophylls = psi.bound_Chorophylls()
                d = np.min(self.distances(LHCII_Chlorophylls,PSI_Chlorophylls))
                if d <=threshold:
                    self.G.add_edge(self.node_num(lhcii.location.T.tolist()[0]),self.node_num(psi.location.T.tolist()[0]),weight=d)

    
        
    @staticmethod
    def distances(M1,M2):
        """returns a 1D matrix of distances between M1 and M2"""
        D = dist.cdist(M1.T,M2.T)
        return np.absolute(np.ndarray.flatten(D))



## graph analysis functions
def PSII_Connectivity(Antenna_Graph):
        """Returns the average number PSII centres per cluster"""
        DATA = {}
        for i in Antenna_Graph.PSII_coordinates:
            N = 1 # default if no other PSIIs connected
            for j in Antenna_Graph.PSII_coordinates:
                if (i != j): # and ((i[0]-j[0])**2 + (i[1]-j[1])**2 <= Antenna_Graph.Distance_Cutoff**2):
                    # I've removed the distance dependence to focus on topology
                    if Antenna_Graph.has_path(i,j):
                        N += 1 # ie i and j are connected
            if N not in DATA.keys():
                DATA[N] = 0
            DATA[N] += 1 # add one cluster of size N
        
        # calculate average cluster size
        NClusterstimesSize = 0
        Nclusters = 0
        for n in DATA.keys():
            NClusterstimesSize += DATA[n]
            Nclusters += DATA[n]/float(n)
        Avg_Cluster_Size = NClusterstimesSize/float(Nclusters)
        return Avg_Cluster_Size

def PSII_Antenna_size(Antenna_Graph,Population):
    """Returns number of LHCIIs for each reaction centre"""
    
    Antenna_Sizes = []
    for i in Population:
        if i.Ptype == "C2S2M2" or i.Ptype == "C2S2":
            i_coord = i.location.T.tolist()[0] # location of j
            if i.Ptype == "C2S2M2":
                N_C2S2M2 = 1
                N_C2S2 = 0
            else:
                N_C2S2M2 = 0
                N_C2S2 = 1
                
            N_LHCII = 0
            for j in Population:
                j_coord = j.location.T.tolist()[0] # location of j
                if i_coord != j_coord: # and ((i_coord[0]-j_coord[0])**2 + (i_coord[1]-j_coord[1])**2 <= Antenna_Graph.Distance_Cutoff**2):
                    # I've removed the distance dependence to focus on topology
                    if Antenna_Graph.has_path(i_coord,j_coord) or Antenna_Graph.has_path(j_coord,i_coord):
                        if j.Ptype == "C2S2M2":
                            N_C2S2M2 += 1
                        elif j.Ptype == "C2S2":
                            N_C2S2 += 1
                        elif j.Ptype == "LHCII":
                            N_LHCII += 1
            Antenna_Sizes.append(4*N_C2S2M2 + 2*N_C2S2 + N_LHCII)
        
    return np.mean(np.array(Antenna_Sizes))       
                        
def PSI_Antenna_size(PSI_Antenna_Graph,Population):
    """Returns number of LHCIIs for each reaction centre"""
    
    Antenna_Sizes = []
    for i in Population:
        if i.Ptype == "PSI":
            i_coord = i.location.T.tolist()[0] # location of j
            N_PSI = 1
                
            N_LHCII = 0
            for j in Population:
                j_coord = j.location.T.tolist()[0] # location of j
                if i_coord != j_coord: # and ((i_coord[0]-j_coord[0])**2 + (i_coord[1]-j_coord[1])**2 <= Antenna_Graph.Distance_Cutoff**2):
                    # I've removed the distance dependence to focus on topology
                    if PSI_Antenna_Graph.has_path(i_coord,j_coord) or PSI_Antenna_Graph.has_path(j_coord,i_coord):
                        if j.Ptype == "PSI":
                            N_PSI += 1
                        elif j.Ptype == "LHCII":
                            N_LHCII += 1
            Antenna_Sizes.append(N_LHCII)
        
    return np.mean(np.array(Antenna_Sizes))  

def Run_analysis(GRANA_RADIUS,DATE,EXPERIMENT):
    ##parameters which describe the model
    ##see Wood et al., 2019/2020 for more details
    
    
    #the following are made global for use in analysis
    global grana_radius
    global NPSII
    global NPSI
    global grana_area
    global stroma_area
    global NC2S2M2
    global NC2S2
    global NB6F_grana
    global NB6F_SL
    global NATP
    global NPSII_SL
    global LHCII_per_RC
    global Complexes
    global Complexes_Chl 
    
    grana_radius = GRANA_RADIUS
    PSI_PSII_ratio = 1
    
    NPSII = round(1.14*(10**(-3))*np.pi*(grana_radius**2))
    NB6F = round(0.5*NPSII)
    B6F_GRANA = 0.5 #fraction of b6f allocated to grana
    NPSII_grana = round(NPSII-(NB6F*B6F_GRANA))
    NPSII = 1.2*NPSII_grana
    
    NPSI= round(PSI_PSII_ratio*2*NPSII)
    stroma_radius =  round(math.sqrt((NPSI/(math.pi*2.45*(10**(-3))))+(grana_radius**2)))
    
    grana_area = math.pi*(grana_radius**2)
    stroma_area = math.pi*(stroma_radius**2)-grana_area
    
    NC2S2 = round(0.5*NPSII_grana) # % total PSII
    NC2S2M2 = round(0.5*NPSII_grana) # % total PSII
    LHCII_per_RC = 5
    LHCII_grana = 0.7
    NLHCIIg = round(LHCII_grana*(NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4)))
    NB6F_grana = round(B6F_GRANA*NB6F)
    
    
    NB6F_SL = round((1-B6F_GRANA)*NB6F)
    NATP = round(0.7*2*NPSII)
    NPSII_SL = round(0.2*NPSII_grana) #2x10% NPSII or 10% reaction centres
    NLHCII_SL =  round((1-LHCII_grana)*(NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4)))
    
    
    Complexes_Chl = pickle.load(open("Complexes_Chl.p",'rb')) # for antenna analysis
    Complexes = pickle.load(open("Complexes.p", 'rb'))
    
    Dir_Name = "./"+EXPERIMENT+"_"+DATE+"/Data/" # directory where the data is
    
    print("Running Analysis of simulation "+EXPERIMENT+"_"+DATE)
    # The dataframe Analysis_Results contains results of all analyses
    Analysis_Results = pd.DataFrame(columns=["PSII_avg_NN","LHCII_avg_NN","PSI_avg_NN","LHCII_grana_fraction","Grana_density","SL_density"])
    
    #Lists to hold the Nearest neighbour distribution data
    PSII_NN = []
    PSII_B6F_NN = []
    LHCII_NN = []
    PSI_NN = []
    
    for t in range(10000000,11000000,2000): # data taken between 10 and 11 M in steps of 2k
         ## Population1##
         #load the data from file
         POPULATION1 = Load_Population(Dir_Name+"POPULATION1_"+str(t))
         
         #run the analysis functions defined elsewhere
         NN_PSII_1 = np.mean(np.array(Nearest_Neighbour_Distances(POPULATION1, ["C2S2M2","C2S2"])))
         NN_LHCII_1 = np.mean(np.array(Nearest_Neighbour_Distances(POPULATION1, ["LHCII"])))
         NN_PSI_1 = np.mean(np.array(Nearest_Neighbour_Distances(POPULATION1, ["PSI"])))
         LHCII_grana_fraction_1 = LHCII_Localisation_Analysis(POPULATION1)
         grana_density_1, sl_density_1 = Density_Analysis(POPULATION1)
         
         #add the results to the dataframe 
         #Analysis_Results = Analysis_Results.append({"PSII_avg_NN":NN_PSII,"LHCII_avg_NN":NN_LHCII,"PSI_avg_NN":NN_PSI,"LHCII_grana_fraction":LHCII_grana_fraction,"Grana_density":grana_density, "SL_density":sl_density},ignore_index=True)
    
         ## Population2##
         #load the data from file
         POPULATION2 = Load_Population(Dir_Name+"POPULATION2_"+str(t))
         
         #run the analysis functions defined elsewhere
         NN_PSII_2 = np.mean(np.array(Nearest_Neighbour_Distances(POPULATION2, ["C2S2M2","C2S2"])))
         NN_LHCII_2 = np.mean(np.array(Nearest_Neighbour_Distances(POPULATION2, ["LHCII"])))
         NN_PSI_2 = np.mean(np.array(Nearest_Neighbour_Distances(POPULATION2, ["PSI"])))
         LHCII_grana_fraction_2 = LHCII_Localisation_Analysis(POPULATION2)
         grana_density_2, sl_density_2 = Density_Analysis(POPULATION2)
         
         # average over the 2 populations
         NN_PSII = np.mean([NN_PSII_1,NN_PSII_2])
         NN_LHCII = np.mean([NN_LHCII_1, NN_LHCII_2])
         NN_PSI = np.mean([NN_PSI_1, NN_PSI_2])
         LHCII_grana_fraction = np.mean([LHCII_grana_fraction_1,LHCII_grana_fraction_2])
         grana_density = np.mean([grana_density_1,grana_density_2])
         sl_density = np.mean([sl_density_1,sl_density_2])
         
         #add the results to the dataframe 
         Analysis_Results = Analysis_Results.append({"PSII_avg_NN":NN_PSII,"LHCII_avg_NN":NN_LHCII,"PSI_avg_NN":NN_PSI,"LHCII_grana_fraction":LHCII_grana_fraction,"Grana_density":grana_density, "SL_density":sl_density},ignore_index=True)
         
         # Need a few whole NN datasets for the NN distributions
         
         if t%50000 == 0: # far too many data points to use all
             PSII_NN += Nearest_Neighbour_Distances(POPULATION1, ["C2S2M2","C2S2"])
             PSII_NN += Nearest_Neighbour_Distances(POPULATION2, ["C2S2M2","C2S2"])
             
             PSII_B6F_NN += Nearest_Neighbour_Distances(POPULATION1, ["C2S2M2","C2S2","B6F"])
             PSII_B6F_NN += Nearest_Neighbour_Distances(POPULATION2, ["C2S2M2","C2S2","B6F"])
             
             
             LHCII_NN += Nearest_Neighbour_Distances(POPULATION1, ["LHCII"])
             LHCII_NN += Nearest_Neighbour_Distances(POPULATION2, ["LHCII"])
             
             PSI_NN += Nearest_Neighbour_Distances(POPULATION1, ["PSI"])
             PSI_NN += Nearest_Neighbour_Distances(POPULATION2, ["PSI"])

    Analysis_Results.to_csv("./"+EXPERIMENT+"_"+DATE+"/Results.csv")
    
    # the nearest neighbour results need to be cut to the same length
    N_results = min([len(PSII_NN),len(LHCII_NN),len(PSII_B6F_NN),len(PSI_NN)]) - 1
    NN_results = pd.DataFrame({"PSII":PSII_NN[:N_results],"PSII_B6F":PSII_B6F_NN[:N_results],"LHCII":LHCII_NN[:N_results],"PSI":PSI_NN[:N_results]})
    NN_results.to_csv("./"+EXPERIMENT+"_"+DATE+"/Nearest_Neighbour_Results.csv")
    #print(Analysis_Results.tail())

def Run_graph_antenna_analysis(GRANA_RADIUS,DATE,EXPERIMENT, PSII=True, PSI=True):
    ##parameters which describe the model
    ##see Wood et al., 2019/2020 for more details
    
    
    #the following are made global for use in analysis
    global grana_radius
    global NPSII
    global NPSI
    global grana_area
    global stroma_area
    global NC2S2M2
    global NC2S2
    global NB6F_grana
    global NB6F_SL
    global NATP
    global NPSII_SL
    global LHCII_per_RC
    global Complexes
    global Complexes_Chl 
    
    grana_radius = GRANA_RADIUS
    PSI_PSII_ratio = 1
    
    NPSII = round(1.14*(10**(-3))*np.pi*(grana_radius**2))
    NB6F = round(0.5*NPSII)
    B6F_GRANA = 0.5 #fraction of b6f allocated to grana
    NPSII_grana = round(NPSII-(NB6F*B6F_GRANA))
    NPSII = 1.2*NPSII_grana
    
    NPSI= round(PSI_PSII_ratio*2*NPSII)
    stroma_radius =  round(math.sqrt((NPSI/(math.pi*2.45*(10**(-3))))+(grana_radius**2)))
    
    grana_area = math.pi*(grana_radius**2)
    stroma_area = math.pi*(stroma_radius**2)-grana_area
    
    NC2S2 = round(0.5*NPSII_grana) # % total PSII
    NC2S2M2 = round(0.5*NPSII_grana) # % total PSII
    LHCII_per_RC = 5
    LHCII_grana = 0.7
    NLHCIIg = round(LHCII_grana*(NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4)))
    NB6F_grana = round(B6F_GRANA*NB6F)
    
    
    NB6F_SL = round((1-B6F_GRANA)*NB6F)
    NATP = round(0.7*2*NPSII)
    NPSII_SL = round(0.2*NPSII_grana) #2x10% NPSII or 10% reaction centres
    NLHCII_SL =  round((1-LHCII_grana)*(NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4)))
    
    
    Complexes_Chl = pickle.load(open("Complexes_Chl.p",'rb')) # for antenna analysis
    Complexes = pickle.load(open("Complexes.p", 'rb'))
    
    Dir_Name = "./"+EXPERIMENT+"_"+DATE+"/Data/" # directory where the data is
    
    print("Running Graph Antenna Analysis of simulation "+EXPERIMENT+"_"+DATE)
    # The dataframe Analysis_Results contains results of all analyses
    PSII_Antenna_Results = pd.DataFrame(columns=["1","2","3","4","5","10"])
    PSII_Connectivity_Results = pd.DataFrame(columns=["1","2","3","4","5","10"])
    PSI_Antenna_Results = pd.DataFrame(columns=["1","2","3","4","5","10"])
    for t in range(10000000,11000000,20000): # data taken between 10 and 11 M in steps of 10k (less than other analysis because takes longer)
         
         ##Population1
         #load the data from file
         POPULATION1 = Load_Population(Dir_Name+"POPULATION1_"+str(t))
         
         #These distionaries store values for different thresolds before adding to data frame 
         PSII_Ant = {}
         PSII_Con = {}
         PSI_Ant = {}
         
         #run the analysis functions defined elsewhere
         
         for THRESHOLD in [1,2,3,4,5,10]:
             if PSII == True:
                 Antenna_graph = Chlorophyll_Network(POPULATION1,threshold=THRESHOLD)
                
                 AS = PSII_Antenna_size(Antenna_graph,POPULATION1)
                 PSII_Ant[str(THRESHOLD)] = AS
                
                 CT = PSII_Connectivity(Antenna_graph)
                 PSII_Con[str(THRESHOLD)] = CT
             
             if PSI == True:
                 PSI_Antenna_graph = PSI_Chlorophyll_Network(POPULATION1,threshold=THRESHOLD)
                 PSI_AS = PSI_Antenna_size(PSI_Antenna_graph,POPULATION1)
                 PSI_Ant[str(THRESHOLD)] = PSI_AS
                
         
         #add the results to the dataframe
         if PSII == True:
            PSII_Antenna_Results = PSII_Antenna_Results.append(PSII_Ant,ignore_index=True)
            PSII_Connectivity_Results = PSII_Connectivity_Results.append(PSII_Con,ignore_index=True)
         if PSI == True:
            PSI_Antenna_Results = PSI_Antenna_Results.append(PSI_Ant,ignore_index=True)
         
         ##Population2
         #load the data from file
         POPULATION2 = Load_Population(Dir_Name+"POPULATION2_"+str(t))
         
         #These distionaries store values for different thresolds before adding to data frame 
         PSII_Ant = {}
         PSII_Con = {}
         PSI_Ant = {}
         
         #run the analysis functions defined elsewhere
         
         for THRESHOLD in [1,2,3,4,5,10]:
             if PSII == True:
                 Antenna_graph = Chlorophyll_Network(POPULATION2,threshold=THRESHOLD)
                
                 AS = PSII_Antenna_size(Antenna_graph,POPULATION2)
                 PSII_Ant[str(THRESHOLD)] = AS
                
                 CT = PSII_Connectivity(Antenna_graph)
                 PSII_Con[str(THRESHOLD)] = CT
                 
             if PSI == True:
                 PSI_Antenna_graph = PSI_Chlorophyll_Network(POPULATION2,threshold=THRESHOLD)
                 PSI_AS = PSI_Antenna_size(PSI_Antenna_graph,POPULATION2)
                 PSI_Ant[str(THRESHOLD)] = PSI_AS
            
         
         #add the results to the dataframe
         if PSII == True:
            PSII_Antenna_Results = PSII_Antenna_Results.append(PSII_Ant,ignore_index=True)
            PSII_Connectivity_Results = PSII_Connectivity_Results.append(PSII_Con,ignore_index=True)
         if PSI == True:
            PSI_Antenna_Results = PSI_Antenna_Results.append(PSI_Ant,ignore_index=True)
    


    
    if PSII == True:
        PSII_Antenna_Results.to_csv("./"+EXPERIMENT+"_"+DATE+"/PSII_Antenna_Results.csv")
        PSII_Connectivity_Results.to_csv("./"+EXPERIMENT+"_"+DATE+"/PSII_Connectivity_Results.csv")
    if PSI == True:
        PSI_Antenna_Results.to_csv("./"+EXPERIMENT+"_"+DATE+"/PSI_Antenna_Results.csv")



if __name__== '__main__':
    
    
    # generally constants
    t0 = time.time()
    GRANA_SIZE = 170 # width of grana, nm.
    Number_of_iterations = 11000001 # number of Monte Carlo steps, Note that data is only collected after 10M iterations.
    DATE = "02_08_19"  # a reference date in which the simulations are run.
    
    #create_initial_population(GRANA_SIZE) # optional usually
    #print("Time elapsed ", (time.time()-t0)/3600.0, "hours")
    
    
    ###Test Simulation###
    
    EXPERIMENT = "Stacking_test"   # a reference for which experiment is being run.
    Stacking_Interaction_Energy = 4 # stacking interaction strength, kT (default = 4).
    LHCII_Binding_Interaction_Energy = 0 # intralayer LHCII interaction strength, kT (default = 2).
    PSI_interaction_energy = 0 # PSI - LHCII interaction strength, kT (default = 0, SII = 2).
    
    POP1 = create_test_population(GRANA_SIZE,Particle_Numbers= [0,0,180, 0,0,20,0, 0, 0])
    POP2 = create_test_population(GRANA_SIZE,Particle_Numbers= [0,0,180, 0,0,20,0, 0, 0])


    POP1, POP2 = Run_Simulation(GRANA_SIZE,DATE,EXPERIMENT,Number_of_iterations,Stacking_Interaction_Energy,LHCII_Binding_Interaction_Energy,PSI_interaction_energy,POPULATION1=POP1,POPULATION2=POP2)
    print("Time elapsed ", (time.time()-t0)/3600.0, "hours")
    
    
    
        
    EXPERIMENT = "SII_test"   # a reference for which experiment is being run.
    Stacking_Interaction_Energy = 0 # stacking interaction strength, kT (default = 4).
    LHCII_Binding_Interaction_Energy = 0 # intralayer LHCII interaction strength, kT (default = 2).
    PSI_interaction_energy = 4 # PSI - LHCII interaction strength, kT (default = 0, SII = 2).
    
    POP1 = create_test_population(GRANA_SIZE,Particle_Numbers= [0,0,180, 0,50,20,0, 0, 0])
    POP2 = create_test_population(GRANA_SIZE,Particle_Numbers= [0,0,180, 0,50,20,0, 0, 0])


    POP1, POP2 = Run_Simulation(GRANA_SIZE,DATE,EXPERIMENT,Number_of_iterations,Stacking_Interaction_Energy,LHCII_Binding_Interaction_Energy,PSI_interaction_energy,POPULATION1=POP1,POPULATION2=POP2)
    print("Time elapsed ", (time.time()-t0)/3600.0, "hours")
    
    
    #Run_analysis(GRANA_SIZE,DATE,EXPERIMENT)
    #print("Time elapsed ", (time.time()-t0)/3600.0, "hours")
    
    #Run_graph_antenna_analysis(GRANA_SIZE,DATE,EXPERIMENT,PSII=True,PSI=True)
      
    #print("Completed in ", (time.time()-t0)/3600.0, "hours")
    

    