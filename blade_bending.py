import pylab as py
import read_files
import ABD_matrix
import load_blade_shape

def calculate_blade_bending():

    # Loading blade section data
    sections = read_files.get_sectional_data()

    # Load blade shape
    c_fun, th_fun, tc_fun, w_fun = load_blade_shape.get_shape_functions()

    # Load blade loads
    r_i = py.array([16.68, 33.46, 51.35, 68.28, 85.43])
    F_i = py.array([210.5, 252.8, 289.0, 316.8, 191.0])*1e3

    # Initial displacement(w), angle(theta) and curvature(kappa)
    w = 0.0
    theta = 0.0
    kappa = 0.0

    r_vec = []
    w_vec = []
    # Calculating section bending
    for i_sec in range(1, sections["n_sec"]):
        # Getting ABD matrix for each section
        sections[i_sec] = ABD_matrix.fib2ABD(sections[i_sec])

        # Calculating r_mean
        sections[i_sec]["r_mean"] = (sections[i_sec]["r_start"]+sections[i_sec]["r_end"])/2.0

        # Calculating Inertia moment (I_xx)
        #Ixx_i = 1 / 12 * w. * t_cap. ^ 3
        sections[i_sec]["I_xx"] = w_fun(sections[i_sec]["r_mean"])*th_fun(sections[i_sec]["r_mean"])**3/12.0

        # Calculating EI
        sections[i_sec]["EI_x"] = sections[i_sec]["E_x"]*sections[i_sec]["I_xx"]

        # Setting start values
        sections[i_sec]["w_start"] = w
        sections[i_sec]["theta_start"] = theta
        sections[i_sec]["kappa_start"] = kappa

        # Assign r_end to a short variable
        r = sections[i_sec]["r_end"]-sections[i_sec]["r_start"]

        # Get equvialent moment
        sections[i_sec]["M_w"], sections[i_sec]["M_theta"],sections[i_sec]["M_kappa"] = get_moment(sections[i_sec]["r_end"], r_i, F_i)

        # Calculate end displacement, angle and curvature
        sections[i_sec]["w_end"] = r**2/(6*sections[i_sec]["EI_x"])*sections[i_sec]["M_w"] + kappa**2/2*r**2+theta*r+w
        sections[i_sec]["theta_end"] = r/(2*sections[i_sec]["EI_x"])*sections[i_sec]["M_theta"]+ kappa*r**2+theta
        sections[i_sec]["kappa_end"] = 1.0/(sections[i_sec]["EI_x"])*sections[i_sec]["M_kappa"]+kappa

        # Assigning temp var for next iteration
        w = sections[i_sec]["w_end"]
        theta = sections[i_sec]["theta_end"]
        kappa = sections[i_sec]["kappa_end"]

        r_vec.append(sections[i_sec]["r_mean"])
        w_vec.append(w)

    py.plot(r_vec,w_vec)
    py.show()
def get_moment(r_in, r_i, F_i):
    M_w = 0.0
    M_theta = 0.0
    M_kappa = 0.0
    for i, r in enumerate(r_i):
        if (r>r_in):
            M_w += (r_in-3*r)*F_i[i]
            M_theta += (r_in-2*r)*F_i[i]
            M_kappa += (r_in-r)*F_i[i]
    return(M_w, M_theta, M_kappa)

if __name__ == "__main__":
    calculate_blade_bending()