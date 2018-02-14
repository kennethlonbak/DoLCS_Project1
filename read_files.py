import pylab as py

fiber_layup_filename_default = r"..\DTU 10MW\structural_models_v1.0\structural_models\composite_layup\composite_layup_Caps.txt"
material_prop_filename_default = r"..\DTU 10MW\structural_models_v1.0\structural_models\ABAQUS\refblade_materials.inp"

# Reading files (fiber layup, material properties)
def read_fiber_layup(fiber_layup_filename=fiber_layup_filename_default):
    # Reading data
    data = py.loadtxt(fiber_layup_filename,skiprows=5)
    # Reading header
    with open(fiber_layup_filename,"r") as file:
        [file.readline() for ii in range(3)]
        header = file.readline().strip().split()


    fiber_layup = []
    for ii in range(len(data)):
        fiber_layup.append({})
        for iii, name in enumerate(header):
            fiber_layup[ii][name] = data[ii,iii]
    return(fiber_layup)

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
