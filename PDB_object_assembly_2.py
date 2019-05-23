import bpy
import math
import csv


Population1_path = r"D:\Thylakoid_model\RANDOM_19_05_19\Data\POPULATION1_11000000"
Population2_path = r"D:\Thylakoid_model\RANDOM_19_05_19\Data\POPULATION2_11000000"

bpy.context.scene.layers[2] = True
bpy.ops.mesh.primitive_cylinder_add(radius=230, depth=3, view_align=False, enter_editmode=False, location=(0, 0, 0), layers=(True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False))
bpy.context.scene.layers[0] = True

if(len(bpy.data.cameras) == 1):
    obj = bpy.data.objects['Camera'] # bpy.types.Camera
    obj.location.x = 0.0
    obj.location.y = -700.0
    obj.location.z = 500.0
    obj.rotation_euler[0] = 50.0
    obj.rotation_euler[1] = 0.0
    obj.rotation_euler[2]  = 5

POPULATION1 = []
with open(Population1_path, 'r') as f1:
    POP1 = csv.reader(f1)
    for row in POP1:
        x = row[1]
        y = row[2]
        theta = row[3]
        ptype = row[4]
        POPULATION1.append([x,y,theta,ptype])
    

POPULATION2 = []
with open(Population2_path, 'r') as f2:
    POP2 = csv.reader(f2)
    for row in POP2:
        x = row[1]
        y = row[2]
        theta = row[3]
        ptype = row[4]
        POPULATION2.append([x,y,theta,ptype])
    
bpy.context.scene.layers[0] = True
    
for row in POPULATION1[1:]:
    x = float(row[0])
    y = float(row[1])
    theta = float(row[2])
    ptype = row[3]    
    z_offset = 0
    scale = 0.1
    theta_offset = 0 # math.pi/4.0
    
    ATP_offset = 10

    if ptype == "C2S2M2":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\C2S2M2_mod.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
        bpy.context.selected_objects[0].rotation_euler[2] = theta+theta_offset + math.pi/2.0
        bpy.context.selected_objects[0].active_material.diffuse_color = (0, 0.273831, 0)
    
    elif ptype == "C2S2":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\C2S2_mod.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
        bpy.context.selected_objects[0].rotation_euler[2] = theta+theta_offset 
        bpy.context.selected_objects[0].active_material.diffuse_color = (0, 0.273831, 0)
    
    elif ptype == "LHCII":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\LHCII_mod2.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
        bpy.context.selected_objects[0].rotation_euler[2] = theta%(2*math.pi)+theta_offset 
        bpy.context.selected_objects[0].active_material.diffuse_color = (0, 0.8, 0.8)
    
    elif ptype == "PSI":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\PSI_mod2.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
        bpy.context.selected_objects[0].rotation_euler[2] = theta%(2*math.pi)+theta_offset  
        bpy.context.selected_objects[0].active_material.diffuse_color = (0.8, 0, 0)

    elif ptype == "B6F":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\B6F_mod.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
        bpy.context.selected_objects[0].rotation_euler[2] = theta%(2*math.pi)+theta_offset 
        bpy.context.selected_objects[0].active_material.diffuse_color = (0.8, 0, 0.8)
    
    elif ptype == "ATP":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\ATP_mod.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset+ATP_offset )
        bpy.context.selected_objects[0].rotation_euler[2] = theta%(2*math.pi)+theta_offset 
    
    elif ptype == "PSII_mono":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\PSIImono_mod.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
    bpy.context.selected_objects[0].rotation_euler[2] = theta%(2*math.pi)+theta_offset 



### Population 2 (upside down) ###


bpy.context.scene.layers[2] = True
bpy.ops.mesh.primitive_cylinder_add(radius=230, depth=3, view_align=False, enter_editmode=False, location=(0, 0, -20.0), layers=(True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False))
bpy.context.scene.layers[1] = True

for row in POPULATION2[1:]:
    x = float(row[0])
    y = float(row[1])
    theta = float(row[2])
    ptype = row[3]    
    z_offset = -20.0
    scale = 0.1
    theta_offset = 0 # math.pi/4.0
    ATP_offset = 10

    if ptype == "C2S2M2":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\C2S2M2_mod.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
        bpy.context.selected_objects[0].rotation_euler[1] = math.pi
        bpy.context.selected_objects[0].rotation_euler[2] = theta+theta_offset + math.pi/2.0
        bpy.context.selected_objects[0].active_material.diffuse_color = (0, 0.273831, 0)
    
    elif ptype == "C2S2":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\C2S2_mod.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
        bpy.context.selected_objects[0].rotation_euler[1] = math.pi
        bpy.context.selected_objects[0].rotation_euler[2] = theta+theta_offset 
        bpy.context.selected_objects[0].active_material.diffuse_color = (0, 0.273831, 0)
    
    elif ptype == "LHCII":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\LHCII_mod2.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
        bpy.context.selected_objects[0].rotation_euler[1] = math.pi
        bpy.context.selected_objects[0].rotation_euler[2] = theta%(2*math.pi)+theta_offset 
        bpy.context.selected_objects[0].active_material.diffuse_color = (0, 0.8, 0.8)
    
    elif ptype == "PSI":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\PSI_mod2.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
        bpy.context.selected_objects[0].rotation_euler[1] = math.pi
        bpy.context.selected_objects[0].rotation_euler[2] = theta%(2*math.pi)+theta_offset  
        bpy.context.selected_objects[0].active_material.diffuse_color = (0.8, 0, 0)
    
    elif ptype == "B6F":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\B6F_mod.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
        bpy.context.selected_objects[0].rotation_euler[1] = math.pi
        bpy.context.selected_objects[0].rotation_euler[2] = theta%(2*math.pi)+theta_offset 
        bpy.context.selected_objects[0].active_material.diffuse_color = (0.8, 0, 0.8)
    
    elif ptype == "ATP":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\ATP_mod.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset-ATP_offset)
        bpy.context.selected_objects[0].rotation_euler[1] = math.pi
        bpy.context.selected_objects[0].rotation_euler[2] = theta%(2*math.pi)+theta_offset 
    
    elif ptype == "PSII_mono":
        bpy.ops.import_scene.obj(filepath=r"F:\programs\Blender_scripting\PSIImono_mod.obj")
        bpy.ops.transform.resize(value=(scale,scale,scale), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.context.selected_objects[0].location = (x,y,z_offset)
        bpy.context.selected_objects[0].rotation_euler[1] = math.pi
        bpy.context.selected_objects[0].rotation_euler[2] = theta%(2*math.pi)+theta_offset 


bpy.ops.render.render()