import scipy as sp
import scipy.linalg as spl

import matplotlib.pyplot as plt

import lattice.latticegeneration as lg
import conf.configloader as conf
import vectorise.pseudovectorised as pv

# Load parameters from configuration
spring_constant = conf.get_configuration_item("latticeGeneration", "springConstant")
mass_default = conf.get_configuration_item("latticeGeneration", "massDefault")
mass_ratio = conf.get_configuration_item("latticeGeneration", "massRatio")
generation = conf.get_configuration_item("latticeGeneration", "generationCount")

# Generate mass crystal
masses = sp.array(lg.generate_masses_1D(generation, mass_ratio, mass_default))

def matrix_multiplication(spring_constant, masses, angular_freq):
    """Evaluate the product of all of the tranfer matrices"""

    # Evaluate all of the a_i coefficients
    a = 2 - ((masses * angular_freq * angular_freq) / spring_constant)

    # Find all of the top rows of the matrices
    top_rows = sp.zeros(len(a) * 2)
    top_rows[::2] = a
    top_rows[1::2] = 1

    # Find all of the bottom rows of the matrices
    bottom_rows = sp.zeros(len(a) * 2)
    bottom_rows[::2] = -1

    # Split the large array into len(a) matrices
    matrices = sp.split(
        sp.matrix(
            sp.vstack((top_rows, bottom_rows))
        ), len(a), axis=1
    )

    # Start with the identity matrix
    current = sp.matrix(sp.identity(2, dtype=float))

    # Can't find a way to vectorise out this loop. Maybe I should
    # restructure the entire function.
    for matrix in matrices:
        current *= matrix

    return current

@pv.pseudo_vectorised("f", (1, 2))
def phase_calculation(angular_freq, spring_constant, masses):
    """Calculates the phase for a given frequency"""

    # Determine the product of all of the transfer matrices
    product_matrix = matrix_multiplication(spring_constant, masses, angular_freq)

    determinant = spl.det(product_matrix)

    #assert sp.fabs(determinant - 1.0) < 1e5, "det(T) = " + str(determinant) + " which is invalid (f = " + str(angular_freq) + ")"
    
    if sp.fabs(sp.trace(product_matrix)) <= 2:
        eigenvalues, eigenvectors = spl.eig(product_matrix)

        # Return the complex argument of the first eigenvalue
        # (where the eigenvalues should be complex conjugates)
        return sp.angle(eigenvalues[0])
    else:
        return sp.nan
        #raise ValueError()

frequencies = sp.arange(0, sp.pi, 0.01)
phases = phase_calculation(frequencies, spring_constant, masses)

plt.figure("Transfer Matrix")

plt.scatter(phases, frequencies, marker=".")

plt.xlabel(r"Phase $\phi$ / $\rm{rad}$")
plt.ylabel(r"Normal Mode Frequency $\omega$ / $\rm{rad} \ \rm{s}^{-1}$")

plt.show()
