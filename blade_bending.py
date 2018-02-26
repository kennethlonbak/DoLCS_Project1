import pylab as py
import read_files
import Sectional_ABD_matrix

def calculate_blade_bending():

    # Loading blade section data
    sections = read_files.get_sectional_data()

    # Load blade shape

    # Load blade loads



    # Calculating ABD
    for i_sec in range(1, sections["n_sec"]):
        sections[i_sec] = Sectional_ABD_matrix.fib_mat2ABD(sections[i_sec])