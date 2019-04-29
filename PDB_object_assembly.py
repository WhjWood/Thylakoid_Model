import bpy
import math
## WhjWood 26/02/19 
## This program is used to construct models in Blender for publication images.

def load_population(filename):
    """ loads in data from a file: x y theta"""
    POP = []
    coord_file = open(filename,'r').readlines()
    for l in coord_file:
        f = l.split(" ")
        POP.append([float(f[0]), float(f[1]),float(f[2])])
    return POP

#paths
#Linux
#path = r'/media/whjw/My Passport/Protein_dynamics/MC_data/LHCII_25_09_18/''
#windows
#path = r'D:\WillW\Thylakoid_Protein_Dynamics\Part2_MC\LHCII_25_09_18\''

#filenames
c2s2m2_file = r'D:\WillW\Thylakoid_Protein_Dynamics\Part2_MC\LHCII_25_09_18\C2S2M2\11_05_18_C2S2M2_Rg_170_9990000usFREE_LHCII'
c2s2_file = r'D:\WillW\Thylakoid_Protein_Dynamics\Part2_MC\LHCII_25_09_18\C2S2\11_05_18_C2S2_Rg_170_9990000usFREE_LHCII'
lhcii_file = r'D:\WillW\Thylakoid_Protein_Dynamics\Part2_MC\LHCII_25_09_18\LHCII\11_05_18_LHCII_Rg_170_9990000usFREE_LHCII'
psi_file = r'D:\WillW\Thylakoid_Protein_Dynamics\Part2_MC\LHCII_25_09_18\PSI\11_05_18_PSI_Rg_170_9990000usFREE_LHCII'
b6fg_file = r'D:\WillW\Thylakoid_Protein_Dynamics\Part2_MC\LHCII_25_09_18\B6Fg\11_05_18_B6F_Rg_170_9990000usFREE_LHCII'
b6fs_file = r'D:\WillW\Thylakoid_Protein_Dynamics\Part2_MC\LHCII_25_09_18\B6Fs\11_05_18_B6Fs_Rg_170_9990000usFREE_LHCII'
atp_file = r'D:\WillW\Thylakoid_Protein_Dynamics\Part2_MC\LHCII_25_09_18\ATP\11_05_18_ATP_Rg_170_9990000usFREE_LHCII'
psiis_file = r'D:\WillW\Thylakoid_Protein_Dynamics\Part2_MC\LHCII_25_09_18\PSIIs\11_05_18_PSIIs_Rg_170_9990000usFREE_LHCII'

#load in data
C2S2M2_POP = load_population(c2s2m2_file)
C2S2_POP = load_population(c2s2_file)
LHCII_POP = load_population(lhcii_file)
PSI_POP = load_population(psi_file)
B6Fg_POP = load_population(b6fg_file)
B6Fs_POP = load_population(b6fs_file)
B6F_POP = B6Fg_POP + B6Fs_POP
ATP_POP = load_population(atp_file)
PSIIs_POP = load_population(psiis_file)



z_offset = 0
scale = 0.1
theta_offset = 0 # math.pi/4.0

for x in C2S2M2_POP:
    bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\C2S2M2_mod.obj")
    bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
    bpy.context.selected_objects[0].location = (x[0],x[1],z_offset)
    bpy.context.selected_objects[0].rotation_euler[2] = x[2]+theta_offset + math.pi/2.0
    bpy.context.selected_objects[0].active_material.diffuse_color = (0, 0.273831, 0)

for x in C2S2_POP:
    bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\C2S2_mod.obj")
    bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
    bpy.context.selected_objects[0].location = (x[0],x[1],z_offset)
    bpy.context.selected_objects[0].rotation_euler[2] = x[2]+theta_offset 
    bpy.context.selected_objects[0].active_material.diffuse_color = (0, 0.273831, 0)

for x in LHCII_POP:
    bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\LHCII_mod2.obj")
    bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
    bpy.context.selected_objects[0].location = (x[0],x[1],z_offset)
    bpy.context.selected_objects[0].rotation_euler[2] = x[2]%(2*math.pi)+theta_offset 
    bpy.context.selected_objects[0].active_material.diffuse_color = (0, 0.8, 0.8)

for x in PSI_POP:
    bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\PSI_mod2.obj")
    bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
    bpy.context.selected_objects[0].location = (x[0],x[1],z_offset)
    bpy.context.selected_objects[0].rotation_euler[2] = x[1]%(2*math.pi)+theta_offset 
    bpy.context.selected_objects[0].active_material.diffuse_color = (0.8, 0, 0)

for x in B6F_POP:
    bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\B6F_mod.obj")
    bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
    bpy.context.selected_objects[0].location = (x[0],x[1],z_offset)
    bpy.context.selected_objects[0].rotation_euler[2] = x[1]%(2*math.pi)+theta_offset 
    bpy.context.selected_objects[0].active_material.diffuse_color = (0.8, 0, 0.8)

for x in ATP_POP:
    bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\ATP_mod.obj")
    bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
    bpy.context.selected_objects[0].location = (x[0],x[1],z_offset+9.2)
    bpy.context.selected_objects[0].rotation_euler[2] = x[1]%(2*math.pi)+theta_offset

for x in PSIIs_POP:
    bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\PSIImono_mod.obj")
    bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
    bpy.context.selected_objects[0].location = (x[0],x[1],z_offset)
    bpy.context.selected_objects[0].rotation_euler[2] = x[1]%(2*math.pi)+theta_offset

