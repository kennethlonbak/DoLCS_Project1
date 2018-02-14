import pylab as py
import read_files

def test_function():
    # Read fiber_layup
    fiber_layup = read_files.read_fiber_layup()[0]
    # Read Material properties
    material_prop = read_files.read_material_properties()

    # Calculate ABD matrix
    fib_mat2ABD(fiber_layup, material_prop)
    return

def fib_mat2ABD(fiber_layup, material_prop):
    # Setting up local compliance matrix (S) for each fiber
    S = get_compliance_matrix(material_prop)

    # Getting local stiffness matrix (Q) - (Q = inv(S))
    Q = get_stiffness_matrix(S)

    # Make transformation matrix for each fiber (T_lg, local->global - inv(T_lg)=T_gl, Global->local)
    pass
    # Make global stiffness matrix

    # Make global A, B, D

    # Collect ABD matrix

# Compliance and stiffness matrices ---------------------------------------------------------------------------------- #
def get_compliance_matrix(material_prop):
    S = {}
    for name in material_prop: # Lecture 3, slide 13
        S_11 = 1.0/material_prop[name]["E1"]
        S_22 = 1.0/material_prop[name]["E2"]
        S_12 = S_21 = -material_prop[name]["nu12"]/material_prop[name]["E1"]
        S_33 = 1.0/material_prop[name]["G12"]
        S[name] = py.array([[S_11,S_12,0],[S_21,S_22,0],[0,0,S_33]])
    return(S)

def get_stiffness_matrix(S):
    Q = {}
    for name in S:
        Q[name] = py.inv(S[name])
    return(Q)


if __name__ == "__main__":
    test_function()


