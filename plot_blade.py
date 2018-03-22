import pylab as py
import blade_bending
from os import path
import ABD_matrix
import plotly

fig_path = ABD_matrix.fig_path

def plot_blade_deflection(use_plotly=True):
    sections = blade_bending.calculate_blade_bending()
    sections_HAWC = blade_bending.calculate_blade_bending(True)

    # Radial vector
    r = blade_bending.sections2value(sections, "r_start", True)

    # Deflection vector
    delta = blade_bending.sections2value(sections, "delta_start", True)
    delta_HAWC = blade_bending.sections2value(sections_HAWC, "delta_start", True)

    # Moment vector
    M = blade_bending.sections2value(sections, "M")

    # Curvature vector
    kappa = blade_bending.sections2value(sections, "kappa")

    r2M = lambda r_in: py.interp(r_in,r[:-1],M)

    fig, ax = py.subplots(2,1, gridspec_kw = {'hspace':0},figsize=(6,4))
    ax[0].plot(r[:-1],M)
    ax[0].quiver(sections["r_i"],r2M(sections["r_i"]),[0]*len(sections["F_i"]),-sections["F_i"],label="Force arrows")
    ax[0].set_ylabel("Moment [Nm]")
    ax[0].grid("on")
    ax[0].legend()
    ax[0].set_ylim([-2e7,7e7])
    ax[1].grid("on")
    ax[1].plot(r,delta,label="Only Spar Caps")
    ax[1].plot(r, delta_HAWC, label="Using HAWC $EI_x$")
    ax[1].plot(r,[-18.26]*len(r),label="Tower clearance",ls="--",lw=0.9)
    ax[1].set_xlabel("Radius [m]")
    ax[1].set_ylabel("Deflection [m]")
    ax[1].legend(loc=0)
    py.tight_layout()
    fig.savefig(path.join(fig_path, "Deflection.png"))

    if use_plotly:
        plotly.offline.plot_mpl(fig)
    else:
        py.show()

def plot_undef_blade(ax=None):
    if ax==None:
        fig, ax = py.subplots(1,1)

    blade = get_undef_blade()

    ax.plot(*blade["tower"],c="k")
    ax.plot(*blade["shaft"],c="k")
    ax.plot(*blade["hub"],c="k",label="Tower")
    ax.plot(*blade["blade"],label="undef. blade")
    ax.axis("equal")
    #py.show()


def get_undef_blade():
    blade = {}
    blade["tower"] = py.array([[0.0, 4.15 / 2, 4.15 / 2, -4.15 / 2, -4.15 / 2, 0.0], [0.0, 0.0, 115.63, 115.63, 0.0, 0.0]])
    blade["shaft"] = py.array([[blade["tower"][0,1],blade["tower"][0,1]-7.1*py.cos(5*py.pi/180)],[blade["tower"][1,2]+2.75,blade["tower"][1,2]+2.75+abs(7.1)*py.sin(5*py.pi/180)]])
    shaft_tan = py.diff(blade["shaft"])
    shaft_tan = shaft_tan[0] +1j*shaft_tan[1]
    shaft_tan /= abs(shaft_tan)
    shaft_normal = shaft_tan*1j

    blade["hub_fun"] = lambda r: blade["shaft"][0,-1]+1j*blade["shaft"][1,-1] +r*shaft_normal

    blade["hub"] = py.array([[py.real(blade["hub_fun"](0)),py.real(blade["hub_fun"](2.8))],[py.imag(blade["hub_fun"](0)),py.imag(blade["hub_fun"](2.8))]])
    cone = -2.5*py.pi/180 # Cone angle
    blade_normal = (py.cos(cone)+1j*py.sin(cone))*shaft_normal
    blade["blade_fun"] = lambda r, R, defl: blade["hub"][0,-1]+1j*blade["hub"][1,-1] +r*blade_normal + r/R*2.332*blade_normal/1j+ defl*blade_normal/1j
    R = 86.366
    blade["blade"] = py.array([[py.real(blade["blade_fun"](0,R,0)),py.real(blade["blade_fun"](R,R,0))],[py.imag(blade["blade_fun"](0,R,0)),py.imag(blade["blade_fun"](R,R,0))]])
    #print(py.angle(blade_normal)*180/py.pi,py.angle(shaft_normal)*180/py.pi)
    return(blade)

def write_baseline_values(i_sec= 50, filename = "info_sec"):
    section = blade_bending.calculate_blade_bending()[i_sec]

    with open("%s%d.dat"%(filename,i_sec),"w") as file:
        for name in section:
            file.write(str(name)+" =")
            file.write(str(section[name]))
            file.write("\n")




if __name__ == "__main__":
    write_baseline_values()
    plot_blade_deflection(False)