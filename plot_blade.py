import pylab as py
import blade_bending
from scipy.integrate import quad

def plot_blade_deflection():
    sections = blade_bending.calculate_blade_bending()

    # Extraction of defelection
    w = []
    r = []
    kappa = []
    for i_sec in range(1, sections["n_sec"]):
        w.append(sections[i_sec]["w_start"])
        r.append(sections[i_sec]["r_start"])
    w.append(sections[i_sec]["w_end"])
    r.append(sections[i_sec]["r_end"])

    w = py.array(w)
    r = py.array(r)
    r = r-r[0]

    turbine_undef = get_undef_blade()
    turbine_def = turbine_undef.copy()
    blade_cmpx =  turbine_def["blade_fun"](r,r[-1],w)
    turbine_def["blade"] = py.array([py.real(blade_cmpx),py.imag(blade_cmpx)])

    fig, ax = py.subplots(2,1)
    plot_undef_blade(ax[0])
    ax[0].plot(*turbine_def["blade"],label="def. blade")
    ax[0].grid("on")
    ax[0].legend(loc=0)
    ax[1].grid("on")
    ax[1].plot(r,w,label="Deformed blade centerline")
    ax[1].plot(r,[-18.26]*len(r),label="Tower clearance",ls="--",lw=0.9)
    ax[1].set_xlabel("Radius [m]")
    ax[1].set_ylabel("Deflection [m]")
    ax[1].legend(loc=0)
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



if __name__ == "__main__":
    plot_blade_deflection()