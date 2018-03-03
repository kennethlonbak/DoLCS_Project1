import pylab as py


def test_blade_shape():
    c_fun, th_fun, tc_fun, tw_fun = get_shape_functions()

    r_new = py.linspace(0,89,1000)

    fig, ax = py.subplots(4,1)
    for ax_temp in ax:
        ax_temp.grid("on")
    ax[0].plot(r_new, c_fun(r_new))
    ax[0].set_ylabel("Chord [m]")

    ax[1].plot(r_new, tc_fun(r_new))
    ax[1].set_ylabel("Relative thickness [\%]")

    ax[2].plot(r_new, th_fun(r_new))
    ax[2].set_ylabel("Thickness [m]")

    ax[3].plot(r_new, tw_fun(r_new))
    ax[3].set_ylabel("Cap width [m]")
    ax[3].set_xlabel("Radius [m]")
    py.show()

def get_shape_functions():
    tck_chord = create_tck_from_file("radius_vs_chord.dat")
    c_fun = lambda r: spline_interpolate(tck_chord, r)

    # Relative thickness (tc)
    tck_tc = create_tck_from_file("radius_vs_relative_thickness.dat")
    tc_fun = lambda r: spline_interpolate(tck_tc, r)

    # thickness
    th_fun = lambda r: tc_fun(r) * c_fun(r)

    # Cap width(tw)
    tck_w = create_tck_from_file("radius_vs_cap_width.dat")
    w_fun = lambda r: spline_interpolate(tck_w, r)
    return c_fun, th_fun, tc_fun, w_fun

def create_tck_from_file(filename):
    data = py.loadtxt(filename,skiprows=1)
    tck_0 = data[:, :2]
    tck_1 = data[:,2:]

    spline_deg = len(tck_1[0])
    tck = (tck_0,tck_1,spline_deg)
    return tck

def spline_interpolate(tck, r_in):
    if not(type(r_in) in [list, py.ndarray]):
        r_in = [r_in]

    out_list = []
    for r in r_in:
        i_intp = py.argmin(abs(r-tck[0][:,0]))

        if (r < tck[0][0,0]):
            out = tck[1][0,-1]
        elif (r>tck[0][-1,-1]):
            out = 0.0#tck[1][-1,-1]
        else:
            for i_intp,[r_min, r_max] in enumerate(tck[0]):
                if (r > r_min and r <= r_max):
                    break
            out = 0
            for ii in range(tck[2]):
                out += (r-tck[0][i_intp,0])**(tck[2]-(ii+1))*tck[1][i_intp,ii]
        out_list.append(out)
    return py.array(out_list)

if __name__ == "__main__":
    test_blade_shape()