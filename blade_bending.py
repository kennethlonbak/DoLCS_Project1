import pylab as py
import read_files
import ABD_matrix
import load_blade_shape
from os import path

fig_path = ABD_matrix.fig_path
baseline_filename = path.join(read_files.win_home,"Users","kloenbaek","Desktop","DoLCS",r"DTU_10MW","aeroelastic_models_v1.0","aeroelastic_models","hawc2","data","DTU_10MW_RWT_Blade_st.dat")

py.rc("font", size=9)

def plot_blade_stiffness():
    sections = calculate_blade_bending()

    # Get r
    r = sections2value(sections, "r_start")

    # Get Ex
    Ex = sections2value(sections, "E_x")

    # Get Ixx
    Ixx = sections2value(sections, "I_xx")

    # Get EI_x
    EI_x = sections2value(sections, "EI_x")

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
    fig.savefig(path.join(fig_path,"Bending_Stiffness.png"))
    py.show()

def calculate_blade_bending(use_HAWC = False):

    # Loading blade section data
    sections = read_files.get_sectional_data()

    if (use_HAWC):
        HAWC_data = read_baseline_data()
        EI_X_intap = lambda r: py.interp(r,HAWC_data["r"],HAWC_data["E"]*HAWC_data["I_x"])

    # Load blade shape
    c_fun, th_fun, tc_fun, w_fun = load_blade_shape.get_shape_functions()

    # Load blade loads
    r_i = py.array([16.68, 33.46, 51.35, 68.28, 85.43]); sections["r_i"] = r_i
    F_i = py.array([210.5, 252.8, 289.0, 316.8, 191.0])*1e3; sections["F_i"] = F_i

    # Calculating section bending
    for i_sec in range(1, sections["n_sec"]):
        # Getting ABD matrix for each section
        sections[i_sec] = ABD_matrix.fib2ABD(sections[i_sec])

        # Section shape
        sections[i_sec]["cap_width"] = w_fun(sections[i_sec]["r_start"])
        sections[i_sec]["height"] = th_fun(sections[i_sec]["r_start"])
        sections[i_sec]["sec_length"] = sections[i_sec]["r_end"]-sections[i_sec]["r_start"]

        # Calculating Inertia moment (I_xx)
        sections[i_sec]["I_xx"] = 1.0/12.0*sections[i_sec]["cap_width"]*(sections[i_sec]["height"]**3-(sections[i_sec]["height"]-2*sections[i_sec]["thickness"])**3)

        # Calculating EI
        if (use_HAWC):
            sections[i_sec]["EI_x"] = EI_X_intap(sections[i_sec]["r_start"])
        else:
            sections[i_sec]["EI_x"] = sections[i_sec]["E_x"]*sections[i_sec]["I_xx"]

        # Calculating moment
        sections[i_sec]["M"] = sum((r_i[r_i > sections[i_sec]["r_start"]] - sections[i_sec]["r_start"]) * F_i[r_i > sections[i_sec]["r_start"]])

        # Setting start values
        sections[i_sec]["kappa"] = -sections[i_sec]["M"] / sections[i_sec]["EI_x"]

    delta_new = kappa_integration(sections)[1]
    for i_sec in range(1, sections["n_sec"]):
        sections[i_sec]["delta_start"] = delta_new[i_sec-1]
        sections[i_sec]["delta_end"] = delta_new[i_sec]
    return(sections)

def kappa_integration(sections):
    r = []
    kappa = []
    for i_sec in range(1, sections["n_sec"]):
        r.append(sections[i_sec]["r_start"])
        kappa.append(sections[i_sec]["kappa"])
    r.append(sections[i_sec]["r_end"])
    kappa.append(0)

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

def sections2value(sections, name, add_end = False):
    values = []
    for i_sec in range(1, sections["n_sec"]):
        values.append(sections[i_sec][name])

    if (add_end):
        name_end = name.replace("start","end")
        values.append(sections[i_sec][name_end])
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