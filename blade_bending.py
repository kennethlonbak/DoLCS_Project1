import pylab as py
import read_files
import ABD_matrix
import load_blade_shape

fig_path = r"C:\Users\kloenbaek\Desktop\DoLCS\Project1\Report\Figures\\"
baseline_filename =r"C:\Users\kloenbaek\Desktop\DoLCS\DTU 10MW\aeroelastic_models_v1.0\aeroelastic_models\hawc2\data\DTU_10MW_RWT_Blade_st.dat"

py.rc("font", size=9)

def plot_blade_stiffness():
    sections = calculate_blade_bending()

    # Get r
    r = sections2value(sections, "r_start")

    # Get Ex
    Ex = sections2value(sections, "E_x")

    # Get Ixx
    Ixx = sections2value(sections, "I_xx")

    # Get Ixx
    EI_x = Ex*Ixx#sections2value(sections, "EI_x")

    # Load baseline values
    data = read_baseline_data()
    #data = py.loadtxt("bending_Stiffness.dat")

    fig, ax = py.subplots(3,1,figsize=(6,4), gridspec_kw = {'hspace':0})

    for temp in ax:
        temp.grid("on",which='both')

    ax[0].semilogy(r,Ex,'.-',label="Own")
    ax[0].semilogy(data["r"], data["E"],label="HAWC2 input")
    ax[0].set_xticklabels([])
    ax[0].set_ylabel(r"$E_x$ [Pa]")

    ax[1].semilogy(r,Ixx,'.-',label="Own")
    ax[1].semilogy(data["r"], data["I_x"],label="HAWC2 input")
    ax[1].set_xticklabels([])
    ax[1].set_ylabel(r"$I_{xx}$ [m$^4$]")

    ax[2].semilogy(r,EI_x,'.-',label="Own")
    ax[2].semilogy(data["r"], data["E"]*data["I_x"],label="HAWC2 input")
    ax[2].set_ylabel(r"$EI_{x}$ [Nm$^2$]")
    ax[2].set_xlabel(r"Radius [m]")
    ax[1].legend(loc=0)
    py.tight_layout()
    fig.savefig(fig_path + "Bending_Stiffness.png")
    py.show()

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

    # Calculating section bending
    for i_sec in range(1, sections["n_sec"]):
        # Getting ABD matrix for each section
        sections[i_sec] = ABD_matrix.fib2ABD(sections[i_sec])

        # Calculating r_mean
        sections[i_sec]["r_mean"] = (sections[i_sec]["r_start"]+sections[i_sec]["r_end"])/2.0

        # Calculating Inertia moment (I_xx)
        sections[i_sec]["I_xx"] = 1.0/12.0*w_fun(sections[i_sec]["r_mean"])*(th_fun(sections[i_sec]["r_mean"])**3-(th_fun(sections[i_sec]["r_mean"])-2*sections[i_sec]["thickness"])**3)
        #sections[i_sec]["I_xx"] = 1.0/12.0*w_fun(sections[i_sec]["r_mean"])*th_fun(sections[i_sec]["r_mean"])**3

        # Calculating EI
        sections[i_sec]["EI_x"] = sections[i_sec]["E_x"]*sections[i_sec]["I_xx"]

        # Setting start values
        sections[i_sec]["w_start"] = w
        sections[i_sec]["theta_start"] = theta
        sections[i_sec]["kappa_start"] = kappa

        # Assign r_end to a short variable
        r = sections[i_sec]["r_end"]-sections[i_sec]["r_start"]

        # Get equivalent moment
        sections[i_sec]["M_w"], sections[i_sec]["M_theta"],sections[i_sec]["M_kappa"] = get_moment(sections[i_sec]["r_start"], r_i, F_i)

        # Calculate end displacement, angle and curvature
        sections[i_sec]["w_end"] = r**2/(6*sections[i_sec]["EI_x"])*sections[i_sec]["M_w"] +w +theta*r #+ kappa/2*r**2
        sections[i_sec]["theta_end"] = r/(2*sections[i_sec]["EI_x"])*sections[i_sec]["M_theta"]+theta #+ kappa*r
        sections[i_sec]["kappa_end"] = 1.0/(sections[i_sec]["EI_x"])*sections[i_sec]["M_kappa"]#+kappa
        #sections[i_sec]["w_end"] = kappa_integration(sections)[1][i_sec-1]

        # Assigning temp var for next iteration
        w = sections[i_sec]["w_end"]
        theta = sections[i_sec]["theta_end"]
        kappa = sections[i_sec]["kappa_end"]
    w_new = kappa_integration(sections)[1]
    for i_sec in range(1, sections["n_sec"]):
        sections[i_sec]["w_start"] = w_new[i_sec-1]
        sections[i_sec]["w_end"] = w_new[i_sec]
    return(sections)

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

def kappa_integration(sections):
    r = []
    kappa = []
    for i_sec in range(1, sections["n_sec"]):
        r.append(sections[i_sec]["r_start"])
        kappa.append(sections[i_sec]["kappa_start"])
    r.append(sections[i_sec]["r_end"])
    kappa.append(sections[i_sec]["kappa_end"])

    r = py.array(r)
    r = r-r[0]
    kappa = py.array(kappa)
    theta_fun = lambda i_in: py.trapz(kappa[:i_in],r[:i_in])
    w_new = []
    theta_new = []
    for i_in in range(1,len(r)+1):
        theta_new.append(theta_fun(i_in))
        w_new.append(py.trapz(theta_new,r[:i_in]))
    return(r, w_new)

def sections2value(sections, name):
    values = []
    for i_sec in range(1, sections["n_sec"]):
        values.append(sections[i_sec][name])
    return(py.array(values))

def read_baseline_data(filename=baseline_filename):
    with open(filename) as file:
        [file.readline() for i in range(3)]
        names = [name.strip() for name in file.readline().strip().split()]
        n_sec = int(file.readline().strip().split()[-1])
        data = {}
        for name in names:
            data[name] = []
        for i_sec in range(n_sec):
            values = [float(value) for value in file.readline().strip().split()]
            for i_name, name in enumerate(names):
                data[name].append(values[i_name])

        for name in names:
            data[name] = py.array(data[name])
    return(data)
if __name__ == "__main__":
    plot_blade_stiffness()