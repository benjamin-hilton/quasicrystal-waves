import matplotlib.pyplot as plt
import scipy as sp

import lattice.latticegeneration as lg
import conf.configloader as conf
import vectorise.pseudovectorised as pv

import math

# Load parameters from configuration
spring_constant = conf.get_configuration_item("latticeGeneration", "springConstant")
mass_default = conf.get_configuration_item("latticeGeneration", "massDefault")
mass_ratio = conf.get_configuration_item("latticeGeneration", "massRatio")
generation = 10 #conf.get_configuration_item("latticeGeneration", "generationCount")
generation_function = str(conf.get_configuration_item("latticeGeneration", "generationFunction"))
method_name = str(conf.get_configuration_item("latticeGeneration", "methodName"))
graph_width = conf.get_configuration_item("latticeGeneration", "size", "width")
graph_height = conf.get_configuration_item("latticeGeneration", "size", "height")


# Generate mass crystal
masses = sp.array(lg.generate_by_name(generation_function, generation, mass_ratio, mass_default))

##def matrix_multiplication(current, matrix):
##    a = current[0] * matrix[0,0] + current[1] * matrix[1,0]
##    b = current[0] * matrix[0,1] + current[1] * matrix[1,1]
##    c = current[2] * matrix[0,0] + current[3] * matrix[1,0]
##    d = current[2] * matrix[0,1] + current[3] * matrix[1,1]
##    return [a,b,c,d]
##
##def find_product(spring_constant, masses, angular_freq):
##    """Evaluate the product of all of the tranfer matrices"""
##
##    # Evaluate all of the a_i coefficients
##    a = 2 - ((masses * angular_freq * angular_freq) / spring_constant)
##
##    # Find all of the top rows of the matrices
##    top_rows = sp.zeros(len(a) * 2)
##    top_rows[::2] = a
##    top_rows[1::2] = 1
##
##    # Find all of the bottom rows of the matrices
##    bottom_rows = sp.zeros(len(a) * 2)
##    bottom_rows[::2] = -1
##
##    # Split the large array into len(a) matrices
##    matrices = sp.split(
##        sp.matrix(
##            sp.vstack((top_rows, bottom_rows))
##        ), len(a), axis=1
##    )
##
##    # Start with the identity matrix
##    current = [1,0,0,1]
##
##    # Can't find a way to vectorise out this loop. Maybe I should
##    # restructure the entire function.
##    for matrix in matrices:
##        current = matrix_multiplication(current, matrix)
##
##        determinant = current[0] * current[3] - current[1] * current[2]
##
##        assert sp.fabs(determinant - 1.0) < 1e-1, "det(T) = " + str(determinant) + " which is invalid (f = " + str(angular_freq) + ") \n\n" + str(current) + "\n \n" + str(current * matrix.I) 
##
##    matrix = sp.matrix([[current[0], current[1]], [current[2], current[3]]])
##    
##    return current

def transfer_matrices_product(masses, spring_constant, angular_freq):
    a_values = map(float, 2 - ((masses * angular_freq ** 2) / spring_constant))

    # Output variables
    tl = 1.0
    tr = 0.0
    bl = 0.0
    br = 1.0

    b = False

    for a in a_values:
        prev = tl, tr, bl, br
        
        tl, tr, bl, br = tl * a - tr, tl, bl * a - br,  bl

        if abs(tl * br - tr * bl) - 1 > 1e-2:
            b = True
            
            #print "DET", tl * br - tr * bl
            #print "PREVIOUS", prev
            #print "TRANSFER", a, 1, -1, 0, "\n\n\n"
            #print "PREV_DET", prev[0] * prev[1] - prev[2] * prev[3]

    #if b:
        #print "\n\n\n\n\n\n\n\n\n"

    return [[tl, tr], [bl, br]]

def phase_calculation(angular_freq, spring_constant, masses):
    """Calculates the phase for a given frequency"""

    # Determine the product of all of the transfer matrices
    product_matrix = transfer_matrices_product(masses, spring_constant, angular_freq)

    trace = product_matrix[0][0] + product_matrix[1][1]
    det = product_matrix[0][0] * product_matrix[1][1] - product_matrix[1][0] * product_matrix[0][1]

    #assert abs(det) - 1 <= 1e-3, "det(T) = " + str(det) + "\n\n" + str(product_matrix) + "\n\n" + str(type(product_matrix[0][0]))
    
    if abs(product_matrix[0][0] + product_matrix[1][1]) <= 2:       
        discriminant = math.sqrt(abs(trace ** 2 - 4 * det)) * 1j

        eigenvalue = (-trace + discriminant) / 2

        # Return the complex argument of the first eigenvalue
        # (where the eigenvalues should be complex conjugates)
        return sp.angle(eigenvalue)
    else:
        return sp.nan
        #raise ValueError()

def plot_fractal(spring_constant, mass_default, mass_ratio, generations, phase):
    plot_generations = []
    plot_angular_freqs = sp.arange(0, 2 * sp.pi, 0.001)

    for generation in generations:
        masses = lg.generate_masses_1D(generation, mass_ratio, mass_default)
        
        for afreq in plot_angular_freqs:
            phase_calculation(afreq, spring_constant, mass)
            

frequencies = sp.arange(0, 2 * sp.pi, 0.001)
phases = []

for frequency in frequencies:
    phases.append(phase_calculation(frequency, spring_constant, masses))
    

plt.figure("Transfer Matrix")

plt.scatter(phases, frequencies, marker=".")

plt.xlabel(r"Phase $\phi$ / $\rm{rad}$")
plt.ylabel(r"Normal Mode Frequency $\omega$ / $\rm{rad} \ \rm{s}^{-1}$")

plt.show()
