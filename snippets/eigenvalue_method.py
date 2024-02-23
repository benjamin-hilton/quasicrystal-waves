import scipy as sp
import scipy.linalg as spl

import lattice.latticegeneration as lg

import matplotlib.pyplot as plt

def create_spring_matrix(N, k):
    m = sp.zeros(shape=(N, N), dtype="complex")

    sp.fill_diagonal(m, -2 * k)

    for i in xrange(N - 1):
        m[i, i + 1] = k

    for i in xrange(1, N):
        m[i, i - 1] = k

    return sp.matrix(m)

def create_mass_matrix(m1, m2, generation):
    diagonal = lg.generate_masses_1D(generation, m2 / m1, m1)

    N = len(diagonal)

    mass_matrix = sp.zeros(shape=(N, N))

    sp.fill_diagonal(mass_matrix, diagonal)

    return N, sp.matrix(mass_matrix)

def find_diagonal_matrix(N, mass_matrix, spring_matrix):
    m = mass_matrix.I * spring_matrix

    eigenvalues, eigenvectors = spl.eig(m)

    eigenvalue_matrix = sp.zeros(shape=(N, N), dtype="complex")
    sp.fill_diagonal(eigenvalue_matrix, eigenvalues)
    eigenvalue_matrix = sp.matrix(eigenvalue_matrix, dtype="complex").real

    eigenvector_matrix = sp.matrix(sp.hstack(eigenvectors).reshape(N, N), dtype="complex").real

    return eigenvalue_matrix, eigenvector_matrix    

def plot_generation_graphs():
    for generations in xrange(0, 10, 2):

        freqs = []
        ks = []

        for k in sp.arange(-20, 20 + 0.1, 0.1):

            N, mass_matrix = create_mass_matrix(1, 4, generations)
            spring_matrix = create_spring_matrix(N, k)

            D, P = find_diagonal_matrix(N, mass_matrix, spring_matrix)

            frequencies = sp.absolute(sp.sqrt(-spl.eig(D)[0]))

            for f in frequencies:
                ks.append(k)
                freqs.append(f)

        plt.figure(str(generations) + "th Generation")

        plt.xlabel(r"Spring Constant $k$ / $\rm{N} \ \rm{m}^{-1}$")
        plt.ylabel(r"Normal Mode Frequency $\omega$ / $\rm{rad} \ \rm{s}^{-1}$")

        plt.scatter(ks, freqs, marker=".")

def plot_fractal():
    k = 10
    m = 1

    for mr in sp.arange(1, 15, 1):

        gs = []
        ws = []

        for generations in xrange(0, 15, 2):
            N, mass_matrix = create_mass_matrix(m, m * mr, generations)
            spring_matrix = create_spring_matrix(N, k)

            D, P = find_diagonal_matrix(N, mass_matrix, spring_matrix)

            frequencies = sp.absolute(sp.sqrt(-spl.eig(D)[0]))

            for f in frequencies:
                gs.append(generations)
                ws.append(f)

        print "Plotting fractal"

        plt.figure("Fractal (mass ratio: " + str(mr) + ")")
        plt.xlabel(r"Allowed Frequencies $\omega$ / $\rm{rad} \ \rm{s}^{-1}$")
        plt.ylabel("Generation")
        plt.scatter(ws, gs)
        
plot_generation_graphs()
#plot_fractal()

plt.show()
