import pylab as py
from os import path
from sys import platform

if "linux" == platform or "linux2" == platform:
    win_home = path.join("/media","kloenbaek","Windows")
else:
    win_home = "C:"

fiber_layup_filename_default = path.join(win_home,"Users","kloenbaek","Desktop","DoLCS",r"DTU_10MW","structural_models_v1.0","structural_models","composite_layup","composite_layup_Caps.txt")
material_prop_filename_default = path.join(win_home,"Users","kloenbaek","Desktop","DoLCS",r"DTU_10MW","structural_models_v1.0","structural_models","ABAQUS","refblade_materials.inp")

def get_sectional_data():
    # Reading initial sectional data
    sections = read_initial_section_fiber_layup()
    sections = change_r_start_to_be_zero(sections)
    # Reading fiber material properties
    mat_prop =  read_material_properties()

    # Setting material properties from fiber type
    sections = set_material_prop_for_lamina(sections, mat_prop)

    # Setting fiber layer coordinate
    sections = set_z_start_z_end(sections)

    # Setting fiber angle
    sections = set_fiber_angle(sections)
    return sections

def set_fiber_angle(sections):
    for i_sec in range(1, sections["n_sec"]):
        for i_fib in range(1, sections[i_sec]["fiber_nr"] + 1):
            sections[i_sec][i_fib]["angle"] = 0.0
    return sections

def set_z_start_z_end(sections):
    for i_sec in range(1, sections["n_sec"]):
        thickness = sections[i_sec]["thickness"]
        z_start = -thickness/2.0
        for i_fib in range(1, sections[i_sec]["fiber_nr"] + 1):
            z_end = z_start + sections[i_sec][i_fib]["thickness"]
            sections[i_sec][i_fib]["z_start"] = z_start
            sections[i_sec][i_fib]["z_end"] = z_end
            z_start = z_end
    return sections

def set_material_prop_for_lamina(sections, mat_prop):
    for i_sec in range(1,sections["n_sec"]):
        for i_fib in range(1,sections[i_sec]["fiber_nr"]+1):
            fib_type = sections[i_sec][i_fib]["fiber_type"]
            for name in mat_prop[fib_type]:
                sections[i_sec][i_fib][name] = mat_prop[fib_type][name]
    return(sections)

def read_material_properties(material_prop_filename=material_prop_filename_default):
    material_prop = {}
    with open(material_prop_filename, "r") as file:
        [file.readline() for ii in range(4)]
        mat_name = [name.strip() for name in file.readline().replace("**","").strip().split(",")]
        mat_name.append(file.readline().replace("**","").strip())

        read_density = False
        read_mat_prop = 0
        for line in file:
            if "MATERIAL" in line:
                lamina = line.strip().split("=")[-1].strip()
                if lamina == "TE_GLUE_MAT":
                    break
                material_prop[lamina] = {}

            if "DENSITY" in line:
                read_density = True
            elif read_density:
                read_density = False
                material_prop[lamina]["density"] = float(line.strip())

            if "ELASTIC" in line:
                read_mat_prop = 1
            elif read_mat_prop == 1:
                values = [float(num) for num in line.strip().replace(",","").split()]
                for mat_ind in range(len(values)):
                    material_prop[lamina][mat_name[mat_ind]] = values[mat_ind]
                read_mat_prop = 2
            elif read_mat_prop == 2:
                mat_ind += 1
                line = float(line.strip())
                material_prop[lamina][mat_name[mat_ind]] = line
                read_mat_prop = 0
    return material_prop

def read_initial_section_fiber_layup(fiber_layup_filename=fiber_layup_filename_default):
    # Reading data
    data = py.loadtxt(fiber_layup_filename,skiprows=5)
    # Reading header
    with open(fiber_layup_filename,"r") as file:
        [file.readline() for ii in range(3)]
        header = file.readline().strip().split()

    sections = {}
    for ii in range(len(data)):
        sec_nr = int(data[ii,0])
        sections[sec_nr] = {}

        fiber_nr = 0
        thickness = 0
        for iii, name in enumerate(header):
            if (name in ["Number","r_start","r_end"]):
                if not(name == "Number"):
                    sections[sec_nr][name] = data[ii,iii]
            else:
                fiber_nr += 1
                sections[sec_nr][fiber_nr] = {}
                sections[sec_nr][fiber_nr]["fiber_type"] = name
                sections[sec_nr][fiber_nr]["thickness"] = data[ii, iii]*1e-3
                thickness += sections[sec_nr][fiber_nr]["thickness"]

        sections[sec_nr]["thickness"] = thickness
        sections[sec_nr]["fiber_nr"]  = fiber_nr
    sections["n_sec"] = len(data)+1
    return(sections)

def change_r_start_to_be_zero(sections):
    r_start = sections[1]["r_start"]
    for i_sec in range(1, sections["n_sec"]):
        sections[i_sec]["r_start"] = sections[i_sec]["r_start"]-r_start
        sections[i_sec]["r_end"] = sections[i_sec]["r_end"]-r_start
    return(sections)

