import pylab as py
import read_files

def test_function():
    # Read fiber_layup
    fiber_layup = read_files.get_sectional_data()[1]

    # Calculate ABD matrix
    fib2ABD(fiber_layup)
    return

def fib2ABD(fiber_layup):
    '''
    :param fiber_layup: type(fiber_layup)=dict, Mandatory: numbers from 1 to n_layer where n_layer is the number of layers. n_layer
    :subparam fiber_layup[i_layer]: type(fiber_layup[i_layer])=dict, Mandatory fields: E1, E2, nu12, G12, angle, thickness, z_start, z_end. Optional fields: fiber_type,
    :return fiber_layup["ABD"] = ABD: ["thickness"] = sum(fiber_layup[i_layer]["thickness"])
    :subfield fiber_layup[i_layer]: ["S_l"] = S_l, ["Q_l"] = Q_l, ["Q_G"] = Q_G, ["A"] = A, ["B"] = B, ["D"] = D
    '''

    A = py.zeros((3,3))
    B = py.zeros((3,3))
    D = py.zeros((3,3))

    for i_fib in range(1, fiber_layup["fiber_nr"] + 1):
        # Setting up local compliance matrix (S) for each fiber
        fiber_layup[i_fib]["S_l"] = get_compliance_matrix(fiber_layup[i_fib])

        # Getting local stiffness matrix (Q) - (Q = inv(S))
        fiber_layup[i_fib]["Q_l"] = py.inv(fiber_layup[i_fib]["S_l"])

        # Make transformation matrix for each fiber (T_lg, local->global - inv(T_lg)=T_gl, Global->local)
        fiber_layup[i_fib]["T_L2G"] = get_transform_Local2Global(fiber_layup[i_fib]["angle"])
        fiber_layup[i_fib]["T_G2L"] = py.inv(fiber_layup[i_fib]["T_L2G"])

        # Make global stiffness matrix
        fiber_layup[i_fib]["Q_G"] = py.dot(fiber_layup[i_fib]["T_L2G"], py.dot(fiber_layup[i_fib]["Q_l"], fiber_layup[i_fib]["T_L2G"].T))

        # Make global A, B, D
        A += fiber_layup[i_fib]["Q_G"]*(fiber_layup[i_fib]["z_end"]-fiber_layup[i_fib]["z_start"])
        B += fiber_layup[i_fib]["Q_G"] *1.0/2.0* (fiber_layup[i_fib]["z_end"]**2 - fiber_layup[i_fib]["z_start"]**2)
        D += fiber_layup[i_fib]["Q_G"] *1.0/3.0* (fiber_layup[i_fib]["z_end"]**3 - fiber_layup[i_fib]["z_start"]**3)

    # Collect ABD matrix
    ABD = py.zeros((6,6))
    ABD[:3,:3] = A
    ABD[:3, 3:] = B
    ABD[3:, :3] = B
    ABD[3:, 3:] = D
    fiber_layup["ABD"] = ABD
    fiber_layup["abd"] = py.inv(fiber_layup["ABD"])
    fiber_layup["E_x"] = 1.0/(fiber_layup["abd"][0,0]*fiber_layup["thickness"])
    fiber_layup["E_y"] = 1.0/(fiber_layup["abd"][1,1]*fiber_layup["thickness"])
    fiber_layup["G_xy"] = 1.0/(fiber_layup["abd"][2,2]*fiber_layup["thickness"])
    return(fiber_layup)


# Compliance and stiffness matrices ---------------------------------------------------------------------------------- #
def get_compliance_matrix(fiber_layup):
    S_11 = 1.0/fiber_layup["E1"]
    S_22 = 1.0/fiber_layup["E2"]
    S_12 = S_21 = -fiber_layup["nu12"]/fiber_layup["E1"]
    S_33 = 1.0/fiber_layup["G12"]
    S_l = py.array([[S_11,S_12,0],[S_21,S_22,0],[0,0,S_33]])
    return(S_l)

def get_transform_Local2Global(angle):
    c = py.cos(angle)
    s = py.sin(angle)
    T_L2G = py.array([[c**2 ,s**2   ,-2*c*s],
                      [s**2 ,c**2   ,2*c*s],
                      [c*s  ,-c*s   ,c**2-s**2]])
    return(T_L2G)

if __name__ == "__main__":
    test_function()


