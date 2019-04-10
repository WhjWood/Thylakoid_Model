import math, random, time
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as dist
import pandas as pd
import pickle

##WHJWood Thylakoid model 2.0 05/03/19
#Date of last edit: 07/03/19




#Particle class definition

class Particle(object):
    
    # class variables (all unshifted particle matrices)
    
    # offsets of PSII-bound LHCIIs
    C2S2M2_LHCIIa = np.matrix([[-10.0,-8.0]]).T
    C2S2M2_LHCIIb = np.matrix([[-10.0,1.0]]).T
    C2S2M2_LHCIIc = np.matrix([[11.0,8.0]]).T
    C2S2M2_LHCIId = np.matrix([[11.0,0.0]]).T 
    C2S2_LHCIIa = np.matrix([[-7.0,7.0]]).T
    C2S2_LHCIIb = np.matrix([[7.0,-7.0]]).T
    
    #make a pickled diction of PDB matrices
    #Unshifted_particle_matrix = {'test' : np.matrix([[-1,1],[0,1],[1,1],[-1,0],[0,0],[1,0],[-1,-1],[0,-1],[1,-1]]).T,}
    
    def __init__(self,x,y,theta, particle_type):
        try:            
            self.location = np.matrix([[int(x)],[int(y)]])
            self.rotation = float(theta)
            self.Ptype = str(particle_type)
            self.Pmatrix = self.particle_matrix()
        except:
            print("Values for particle location, rotation or type are of incorrect format")
    
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
        

    
    @staticmethod
    def rotate(theta,A, axis):
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
        try:            
            self.location = Location
            self.rotation = float(theta)
            self.Ptype = str(particle_type)
            self.Pmatrix = self.particle_matrix()
        except:
            print("Values for particle location, rotation or type are of incorrect format")

# #### Functions for overlap detection    
# 


def unique_2D_mapping(A,C=10000):
        """converts a 2D matrix to a 1d vector"""
        return np.array(A[0,:]+C*A[1,:])[0] # this maps x and y onto a unique z and returns array 

def N_unique(A):
    """returns the number of unique columns in matrix A"""
    MIN = np.min(A)
    return pd.unique(unique_2D_mapping(A-MIN)).shape[0]

    
# def remove_repeat_cols(A):
#     if A.shape[0] == 2:
#         C = 10000
#         MIN = np.min(A)
#         MIN_M = np.repeat(np.matrix([[MIN],[MIN]]),A.shape[1],axis=1)
#         unique_A = pd.unique(unique_2D_mapping(A-MIN_M))
#         X = np.matrix(unique_A%C)
#         Y = np.matrix(MIN+ (unique_A-X)/C)
#         X+=MIN
# 
#         new_A = np.concatenate((Y,X),axis=0)
#         return new_A
##
def remove_repeat_cols(A):
    return np.matrix(np.unique(A,axis=1))


def collision2(X,P):
    return np.unique(unique_2D_mapping(np.concatenate((X,P),axis=1))).shape[0] < X.shape[1]+P.shape[1] 

def collision3(X,P):
    return N_unique(np.concatenate((X,P),axis=1)) < X.shape[1]+P.shape[1] 


###utility functions
def list_to_matrix(coords):
    return np.matrix(coords).T


def rotate(theta,A, axis):
        return np.round(np.matrix([[np.cos(theta),-np.sin(theta)],[np.sin(theta), np.cos(theta)]])*np.matrix(A))


def population_matrix(X):
    """X is a list of Particle objects"""
    return np.concatenate([x.Pmatrix for x in X],axis=1)

 
def coordinate_matrix(X):
    return np.concatenate([np.matrix(x) for x in X],axis=1)
# 
# def population_tuple(X,P):
#     return [particle_matrix(x,P) for x in X]
# 
#
 
def pdb_to_matrix(file_name):
    particle_file = open(file_name,'r')
    COORDS = []
    for line in particle_file:
        if line[0:4] == "ATOM" or line[0:6] == "HETATM" or line[0:6] == "ANISOU":
            # divide by 10 to get nm or leave as is for Angstroms
            COORDS.append((int(round(float(line[30:38])/10)),int(round(float(line[38:46])/10))))#,int(float(line[46:54])/10.0)))

    X = remove_repeat_cols(list_to_matrix(COORDS))
    return (X-X.mean(1)).astype(int)
# 
# 
# def matrix_difference(X1,X2):
#     """removes columns of X2 from X1"""
#     S1 = set([tuple(row) for row in X1.T.tolist()])
#     return np.matrix([list(X) for X in (S1.difference(set([tuple(row) for row in X2.T.tolist()]))) ]).T
# 
# 
# def create_initial_population(BOUNDARIES,N,P):
#     """adds N particles (P) to the population
#     boundaries are in matrix form"""
#      
#     POPULATION = []
#     while len(POPULATION)<N:
# 
#         
#         x =  np.random.randint(0,2*grana_radius)-grana_radius
#         y =  np.random.randint(0,2*grana_radius)-grana_radius
#         
#         while np.sqrt(x**2 + y**2) >= grana_radius:
#             x =  np.random.randint(0,2*grana_radius)-grana_radius
#             y =  np.random.randint(0,2*grana_radius)-grana_radius
#             
#         coord = [np.matrix([[x],[y]]),0.0, 0.0]
#         #coord = [np.matrix([[x],[y]]),0.0, np.random.random()*2*np.pi]                                                   
#         
#         if len(POPULATION)==0:#only boundaries taken into account
#             
#             COLL = collision2(BOUNDARIES,particle_matrix(coord,P))
# 
#             while COLL == True:
#                 dxdy = np.matrix(random.choice([[0,1],[0,-1],[1,0],[-1,0]])).T
#                 
#                     
#                 d0 = 0#random.choice([np.pi/12,-1*np.pi/12])
#                 #d1 = random.choice([np.pi/12,-1*np.pi/12])
#                 coord1 = [coord[0]+dxdy,coord[1]+d0]#,coord[2]+d1]
#                 
#                 if np.sqrt(coord1[0][0,0]**2 + coord1[0][1,0]**2) < grana_radius:
#                     coord = coord1
#                 COLL = collision2(BOUNDARIES,particle_matrix(coord,P))
#                 
#         
#         else:# other particles involved
#             
#             COLL = collision2(np.concatenate((population_matrix(POPULATION,P),BOUNDARIES),axis =1),particle_matrix(coord,P))
# 
#             while COLL == True:
#                 dxdy = 2*np.matrix(random.choice([[0,1],[0,-1],[1,0],[-1,0]])).T
# 
#                     
#                 d0 = 0#random.choice([np.pi/12,-1*np.pi/12])
#                 #d1 = random.choice([np.pi/12,-1*np.pi/12])
#                 #coord = [coord[0]+dxdy,coord[1]+d0]#coord[2]+d1]
#                 coord1 = [coord[0]+dxdy,coord[1]+d0]#,coord[2]+d1]
#                 
#                 if np.sqrt(coord1[0][0,0]**2 + coord1[0][1,0]**2) < grana_radius:
#                     coord = coord1
#                 COLL = collision2(np.concatenate((population_matrix(POPULATION,P),BOUNDARIES),axis =1),particle_matrix(coord,P))
#         
#         POPULATION.append(coord)
#         print(len(POPULATION))
#         #print(pd.unique(unique_2D_mapping(population_matrix(POPULATION,P))).shape[0]/population_matrix(POPULATION,P).shape[1])
#     return POPULATION




def particle_step2(POPULAT,BOUNDARIES,PARTICLE):
    
    
    """in this version of the function, POPULATION is a matrix not a list"""
    dxdy = np.matrix(random.choice([[0,1],[0,-1],[1,0],[-1,0]])).T                 
    d0 = random.choice([np.pi/33.0,-1*np.pi/33.0])

    
    PARTICLE2 = SudoParticle(PARTICLE.location+dxdy,PARTICLE.rotation+d0,PARTICLE.Ptype)   
    #coord = [x0+np.matrix(random.choice([[0,1,0],[0,-1,0],[1,0,0],[-1,0,0],[0,0,1],[0,0,-1]])).T,(th0+random.choice([np.pi/12,-1*np.pi/12]))%2*np.pi,(th0+random.choice([np.pi/12,-1*np.pi/12]))%2*np.pi ] 
    COLL = collision3(np.concatenate((POPULAT,BOUNDARIES),axis =1),PARTICLE2.Pmatrix)
    BOUND = boundary_check(PARTICLE.Ptype,PARTICLE.location[0,0],dxdy[0,0],PARTICLE.location[1,0],dxdy[1,0])
    return (COLL or BOUND), PARTICLE2   


def time_step(Population1,Population2,L,BOUNDARIES,layer,current_energy):
    """chooses a random particle from population.Population is a list of particles. L is length of population"""
    if layer == 1:
        Population = Population1
        other_population = Population2
    else:
        Population = Population2
        other_population = Population2
    x = np.random.randint(0,L) # index of particle
    coll, new_particle = particle_step2(population_matrix(Population[0:x]+Population[x+1:]),BOUNDARIES,Population[x])
    
    if coll == False:
        new_energy = Hamiltonian(Population[0:x]+[new_particle]+Population[x+1:],other_population)
        if ((new_particle.Ptype =="C2S2M2") or (new_particle.Ptype =="C2S2") or (new_particle.Ptype =="LHCII")):
            dE = new_energy - current_energy # change in energy
        else:
            dE=0
        
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


# def time_step2(T,POPULATION,BOUNDARIES,P):
#     """ new version uses population tuple and concetenate, but works with a list"""
# 
#     for t in np.arange(T):
#         POP = population_tuple(POPULATION,P)
#         order = np.random.permutation(len(POPULATION)) # goes only up to n -1
#         
#         for i in order:
#             
#             POP2 = tuple(POP[0:i]+POP[i+1:])
#             COLL,coord = particle_step2(np.concatenate(POP2,axis=1),BOUNDARIES,POPULATION[i],P)
# 
#             if COLL == False:
#                 POP[i] = particle_matrix(coord,P)
#                 POPULATION[i] = coord
# 
#         
#     return POPULATION        
#          
# 
# def NN_vals(pop):
#     T = [np.array(x[0]) for x in pop]
#     for t1 in T:
#         vals = []
#         for t2 in T:
#             if list(t1)!=list(t2):
#                 vals.append(np.sqrt((t1[0][0]-t2[0][0])**2+(t1[1][0]-t2[1][0])**2))
#         print(min(vals))                  
# 
# def display_coordinates(pop):
#     T = [np.array(x[0]) for x in pop]
#     for t in T:
#         print(str(t[0][0])+" "+str(t[1][0]))
# 
# def save_coordinates(pop,filename):
#     
#     file1 = open(filename, 'w')
#     for p in pop:
#         t = np.array(p[0])
#         file1.write(str(t[0][0])+" "+str(t[1][0])+" "+str(p[1])+"\n")
# 
#     file1.close()
# 
# def load_coordinates(filename):
#     """loads from the files: line: x y theta"""
#     coords = open(filename,'r').readlines()
#     POP = []
#     for coord in coords:
#         f=coord.split()
#         POP.append([np.matrix([[float(f[0])],[float(f[1])]]),float(f[2])])
#     
#     return POP 


def load_coordinates(filename, ptype):
    """loads from the files: line: x y theta"""
    coords = open(filename,'r').readlines()
    POP = []
    for coord in coords:
        f=coord.split()
        POP.append(Particle(float(f[0]),float(f[1]),float(f[2]),ptype))
    
    return POP 

#         
#    
# def save_matrix(filename,m):
#     file1=open(filename,'w')
#     if m.shape[0] == 2:
#         
#         for i in np.arange(m.shape[1]):
#             file1.write(str(m[0,i])+" "+str(m[1,i])+"\n")
#     if m.shape[0] == 3:
#         
#         for i in np.arange(m.shape[1]):
#             file1.write(str(m[0,i])+" "+str(m[1,i])+" "+str(m[2,i])+"\n")
#     file1.close()
# 
# 
#


def distances(M1,M2):
    """returns a 1D matrix of distances between M1 and M2"""
    D = dist.cdist(M1.T,M2.T)
    return np.absolute(np.ndarray.flatten(D))

#jit(nopython=True)
def intra_layer(D):
    """calculates intralayer potentials from 1D matrix of distances"""
    return E_intra*(D[(2*LHCII_radius<=D)&(D<=(2*LHCII_radius+3))].shape[0]) # may need to tighten the threshold from 6+3 nm m

#@jit(nopython=True)
def inter_layer(X1,X2):
    """X1 and X2 are the LHCII locations in layers 1 and 2"""
    D = distances(X1,X2)
    return np.sum(E_inter*np.square((D[D<LHCII_radius]-LHCII_radius)/LHCII_radius))


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

# This function cannot be jit compiled right now
def Hamiltonian(Population1,Population2):
    LHCIIs_Layer1 = np.concatenate(tuple(p.bound_LHCII() for p in Population1 if (p.Ptype=="C2S2M2") or (p.Ptype=="C2S2") or (p.Ptype=="LHCII")), axis=1)
    LHCIIs_Layer2 = np.concatenate(tuple(p.bound_LHCII() for p in Population1 if (p.Ptype=="C2S2M2") or (p.Ptype=="C2S2") or (p.Ptype=="LHCII")), axis=1)
    
    # Distances between LHCIIs in each layer
    D1 = distances(LHCIIs_Layer1,LHCIIs_Layer1)
    D2 = distances(LHCIIs_Layer2,LHCIIs_Layer2)
    
    # intra particle interaction energy in each layer
    Eintra1 = intra_layer(D1)
    Eintra2 = intra_layer(D2)
    
    # interlayer interaction energy (stacking interaction)
    Einter1 = inter_layer(LHCIIs_Layer1,LHCIIs_Layer2)
    
    return Eintra1 + Eintra2 + Einter1


# t0 = time.time()
# 
# 
# """
#PDB files
# C2S2_PDB = pdb_to_matrix("3JCU.pdb")
# C2S2M2_PDB = pdb_to_matrix("5MDX.pdb")
# B6F_PDB = pdb_to_matrix("1UM3_FULL.pdb")
# LHCII_PDB =  pdb_to_matrix("2BHW.pdb")
# PSI_PDB = pdb_to_matrix("4Y28_2.pdb")
# ATP_PDB = pdb_to_matrix("5ARA_memb.pdb")
# PSII_mono_PDB = pdb_to_matrix("3JCU_MONO.pdb")
# 
# Complexes = {"C2S2M2": C2S2M2_PDB,"C2S2": C2S2_PDB, "B6F" : B6F_PDB, "LHCII" : LHCII_PDB, "PSI" : PSI_PDB, "ATP" : ATP_PDB, "PSII_mono" : PSII_mono_PDB}
# pickle.dump(Complexes, open("Complexes.p", 'wb'))


def Run_Simulation(GRANA_RADIUS=170,DATE="11_05_18",EXPERIMENT="FREE_LHCII",TMAX=100):
    #main
    ##parameters
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
    
    grana_area = math.pi*(grana_radius**2)
    stroma_area = math.pi*(stroma_radius**2)-grana_area
    
    NC2S2 = round(0.5*NPSII_grana) # % total PSII
    NC2S2M2 = round(0.5*NPSII_grana) # % total PSII
    LHCII_per_RC = 5
    LHCII_grana = 0.7
    NLHCIIg = round(LHCII_grana*(NPSII*2*LHCII_per_RC-(NC2S2*2+NC2S2M2*4)))
    NB6F_grana = round(B6F_GRANA*NB6F)
    
    #see lab book for details
    #stroma_radius =  round(math.sqrt((NPSI/(math.pi*2.45*(10**(-3))))+(grana_radius**2)))
    
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
    
    LHCII_radius = 3.5 # nm
    interaction_radius = 5.5 # nm
    E_intra = 2.0 #KT , intralayer interactions
    E_inter = 4.0 #KT, interactions over the stromal gap
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
    
    # # first: initialise the populations and M matrices
    date=DATE
    exp=EXPERIMENT
    start= "_"+str(0)+"us"
    
    # Old way of loading in populations is replaced with
    POPULATION1 = pickle.load(open("Population1_initial", 'rb'))
    POPULATION2 = pickle.load(open("Population2_initial", 'rb'))
    """
    ### LAYER 1
    #Grana
    #C2S2M2
    POP_C2S2M2 = load_coordinates("./C2S2M2/"+date+"_C2S2M2_Rg_"+str(grana_radius)+start+exp,"C2S2M2")
    POP_C2S2 = load_coordinates("./C2S2/"+date+"_C2S2_Rg_"+str(grana_radius)+start+exp,"C2S2")
    POP_B6F_g = load_coordinates("./B6Fg/"+date+"_B6F_Rg_"+str(grana_radius)+start+exp,"B6F")
    POP_LHCII_ALL = load_coordinates("./LHCII/"+date+"_LHCII_Rg_"+str(grana_radius)+start+exp, "LHCII")
    POP_PSI = load_coordinates("./PSI/"+date+"_PSI_Rg_"+str(grana_radius)+start+exp,"PSI")
    POP_B6F_s = load_coordinates("./B6Fs/"+date+"_B6Fs_Rg_"+str(grana_radius)+start+exp,"B6F")
    POP_ATP = load_coordinates("./ATP/"+date+"_ATP_Rg_"+str(grana_radius)+start+exp,"ATP")
    POP_PSII_SL = load_coordinates("./PSIIs/"+date+"_PSIIs_Rg_"+str(grana_radius)+start+exp,"PSII_mono")
    
    ### LAYER 2
    POP_C2S2M2_2 = load_coordinates("./C2S2M2/"+date+"_C2S2M2_2_Rg_"+str(grana_radius)+start+exp,"C2S2M2")
    POP_C2S2_2 = load_coordinates("./C2S2/"+date+"_C2S2_2_Rg_"+str(grana_radius)+start+exp,"C2S2")
    POP_B6F_g_2 = load_coordinates("./B6Fg/"+date+"_B6F_2_Rg_"+str(grana_radius)+start+exp,"B6F")
    POP_LHCII_ALL_2 = load_coordinates("./LHCII/"+date+"_LHCII_2_Rg_"+str(grana_radius)+start+exp,"LHCII")
    POP_PSI_2 = load_coordinates("./PSI/"+date+"_PSI_2_Rg_"+str(grana_radius)+start+exp,"PSI")
    POP_B6F_s_2 = load_coordinates("./B6Fs/"+date+"_B6Fs_2_Rg_"+str(grana_radius)+start+exp,"B6F")
    POP_ATP_2 = load_coordinates("./ATP/"+date+"_ATP_2_Rg_"+str(grana_radius)+start+exp,"ATP")
    POP_PSII_SL_2 = load_coordinates("./PSIIs/"+date+"_PSIIs_2_Rg_"+str(grana_radius)+start+exp,"PSII_mono")
    
    
    print("Initialised")
    POPULATION1 = POP_C2S2M2+POP_C2S2+POP_LHCII_ALL+POP_B6F_g+POP_B6F_s+POP_PSII_SL+POP_ATP+POP_PSI
    POPULATION2 = POP_C2S2M2_2+POP_C2S2_2+POP_LHCII_ALL_2+POP_B6F_g_2+POP_B6F_s_2+POP_PSII_SL_2+POP_ATP_2+POP_PSI_2
    """
    
    Nparticles = len(POPULATION1)
    
    t0 = time.time()
    tstart = 0
    Tmax= TMAX
    ENERGY = [Hamiltonian(POPULATION1,POPULATION2)]
    
    # M = population_matrix(POPULATION1)
    # fig=plt.figure()
    # ax=fig.add_subplot(111)
    # ax.set_xlabel("Distance (nm)")
    # plt.axis([-300,300,-300,300])
    # plt.axis('equal')
    # plt.axes().set_aspect('equal')
    # plt.scatter(np.array(M[0,:]),np.array(M[1,:]),s=0.1)


    X_array = np.random.random(Tmax)
    for i in np.arange(Tmax-tstart):
        t=i+tstart
    
        
        # the following is a decision tree that leads to the perturbation
        X = X_array[i]
        if X<0.5: # Layer one
            layer = 1
            POPULATION1, new_energy = time_step(POPULATION1,POPULATION2, Nparticles, BOUNDARIES_LHCII, layer, ENERGY[-1])
            ENERGY.append(new_energy)
        else: # Layer one
            layer = 2
            POPULATION2 , new_energy = time_step(POPULATION1, POPULATION2,Nparticles,BOUNDARIES_LHCII, layer, ENERGY[-1])
            ENERGY.append(new_energy)
    
            
    # M = population_matrix(POPULATION1)
    # plt.scatter(np.array(M[0,:]),np.array(M[1,:]),s=0.1)
    # plt.show()
    plt.plot(np.arange(len(ENERGY)),np.array(ENERGY))
    return POPULATION1, POPULATION2
              

