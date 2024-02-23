import scipy as sp
import scipy.linalg as spl

import matplotlib.pyplot as plt

import lattice.latticegeneration as lg
import conf.configloader as conf

import math
import cmath

import os

plt.rcParams['figure.max_open_warning'] = 50

#Load parameters from configuration
spring_constant = conf.get_configuration_item("latticeGeneration", "springConstant")
mass_default = conf.get_configuration_item("latticeGeneration", "massDefault")
mass_ratio = conf.get_configuration_item("latticeGeneration", "massRatio")
generation = conf.get_configuration_item("latticeGeneration", "generationCount")
generation_function = str(conf.get_configuration_item("latticeGeneration", "generationFunction"))
method_name = str(conf.get_configuration_item("latticeGeneration", "methodName"))
graph_size = conf.get_configuration_item("latticeGeneration", "size").values()


def create_spring_matrix(num_masses, spring_constant, phi):
    """Create the matrix representing the equations of motion"""
    
    m = sp.zeros(shape=(num_masses, num_masses), dtype="complex")

    # Fill the leading diagonal
    sp.fill_diagonal(m, -2 * spring_constant)

    # Fill the upper diagonal
    for i in xrange(num_masses - 1):
        m[i, i + 1] = spring_constant

    # Fill the lower diagonal
    for i in xrange(1, num_masses):
        m[i, i - 1] = spring_constant

    # Fill the corner elements with a complex phase
    m[-1, 0] = spring_constant * cmath.rect(1, -phi)
    m[0, -1] = spring_constant * cmath.rect(1, phi)

    return sp.matrix(m)

def create_mass_matrix(mass_ratio, mass_default, generation):
    """Create a matrix which has the masses on the leading diagonal"""

    diagonal = lg.generate_by_name(generation_function, generation, mass_ratio, mass_default)

    #print diagonal
    
    num_masses = len(diagonal)

    mass_matrix = sp.zeros(shape=(num_masses, num_masses))

    sp.fill_diagonal(mass_matrix, diagonal)

    return num_masses, sp.matrix(mass_matrix)    
    
def frequency_calculation(num_masses, spring_constant, phase, mass_matrix):
    spring_matrix = create_spring_matrix(num_masses, spring_constant, phase)
    
    eigenvalues = -spl.eig(mass_matrix.I * spring_matrix)[0]
    
    return sp.all(eigenvalues.imag <= 1e-5), eigenvalues.real

phases = sp.arange(0, sp.pi, 0.005)

def plot_omega_phi(spring_constant, mass_default, mass_ratio, generation, phases):
    num_masses, mass_matrix = create_mass_matrix(mass_ratio, mass_default, generation)

    plot_phases = []
    plot_angular_freqs = []
    
    for phase in phases:
        real, angular_freqs = frequency_calculation(num_masses, spring_constant, phase, mass_matrix)

        if real:
            for angular_freq in angular_freqs:
                plot_phases.append(phase)
                plot_angular_freqs.append(angular_freq)


    fig = plt.figure("Omega vs. Phi - Eigenvalue (" + str(generation) + " generations, " + method_name + ")", figsize=graph_size)

    plt.axis([0 - sp.pi / 10, 2.1 * sp.pi, -0.5, math.ceil(max(plot_angular_freqs))])

    plt.scatter(plot_phases, plot_angular_freqs, marker=".")

    plt.xlabel(r"Phase $\phi$ / $\rm{rad}$")
    plt.ylabel(r"Normal Mode Frequency $^2$ $\omega^2$ / $\rm{rad}^2 \ \rm{s}^{-2}$")

    fig.savefig(os.path.dirname(__file__) + os.sep + os.pardir + os.sep + os.pardir + os.sep + os.pardir + os.sep + "images" + os.sep + "Omega vs. Phi - Eigenvalue (" + str(generation) + " generations, " + method_name + ").png")
    

generations = sp.arange(1, generation)

def plot_multiple_generations(spring_constant, mass_default, mass_ratio, generations, phase, plot_func):
    for generation in generations:
        print "GENERATION:", generation
        plot_func(spring_constant, mass_default, mass_ratio, generation, phase)

mass_ratios = sp.arange(0.001, 10.01, 0.01)
frequencies = sp.arange(0.001, sp.pi * 2, 0.01)

def plot_omega_mass(spring_constant, mass_default, mass_ratios, generation, phase):
    plot_mass_ratios = []
    plot_angular_freqs = []
    
    for mass_ratio in mass_ratios:
        num_masses, mass_matrix = create_mass_matrix(mass_ratio, mass_default, generation)
        real, angular_freqs = frequency_calculation(num_masses, spring_constant, phase, mass_matrix)

        if real:
            for angular_freq in angular_freqs:
                plot_mass_ratios.append(mass_ratio)
                plot_angular_freqs.append(angular_freq) 

    fig = plt.figure("Omega vs. Mass - Eigenvalue (" + str(generation) + " generations, " + method_name + ")", figsize=graph_size)

    plt.axis([-2, 12, -0.5, 8])

    plt.scatter(plot_mass_ratios, plot_angular_freqs, marker=".")

    plt.xlabel(r"Mass ratio")
    plt.ylabel(r"Normal Mode Frequency $^2$ $\omega^2$ / $\rm{rad}^2 \ \rm{s}^{-2}$")

    fig.savefig(os.path.dirname(__file__) + os.sep + os.pardir + os.sep + os.pardir + os.sep + os.pardir + os.sep + "images" + os.sep + "Omega vs. Mass - Eigenvalue (" + str(generation) + " generations, " + method_name + ").png")
    
def plot_fractal(spring_constant, mass_default, mass_ratio, generations, phase):
    plot_generations = []
    plot_angular_freqs = []
    
    for generation in generations:
        print "GENERATION:", generation
        num_masses, mass_matrix = create_mass_matrix(mass_ratio, mass_default, generation)
        real, angular_freqs = frequency_calculation(num_masses, spring_constant, phase, mass_matrix)
        
        if real:
            plot_angular_freqs += list(angular_freqs)
            plot_generations += list(sp.ones(angular_freqs.size) * generation)
    
    fig = plt.figure("Eigenvalue Fractal" + "(" + method_name + ")", figsize=graph_size)
   
    plt.scatter(plot_angular_freqs, plot_generations, marker = "|", s = 150)
    
    plt.ylabel(r"Generation")
    plt.xlabel(r"Normal Mode Frequency $^2$ $\omega^2$ / $\rm{rad}^2 \ \rm{s}^{-2}$")
    
    fig.savefig(os.path.dirname(__file__) + os.sep + os.pardir + os.sep + os.pardir + os.sep + os.pardir + os.sep + "images" + os.sep + "Eigenvalue Fractal (" + method_name + ").png")
     
def plot_eigenvector_fractal(spring_constant, mass_default, mass_ratio, generations, phase):
    plot_generations = []
    plot_components = []

    for generation in generations:
        num_masses, mass_matrix = create_mass_matrix(mass_ratio, mass_default, generation)
        
        spring_matrix = create_spring_matrix(num_masses, spring_constant, phase)
        eigenvalues, eigenvectors = spl.eig(mass_matrix.I * spring_matrix)

        if sp.all(eigenvalues.imag <= 1e-5):
            print "GENERATION:", generation
            print eigenvalues[sp.argmax(-eigenvalues.real)]
            for x in eigenvectors[sp.argmax(eigenvalues.real)]:
                plot_generations.append(generation)
                plot_components.append(sp.absolute(x))

    fig = plt.figure("Eigenvector Fractal " + "(" + method_name + ")", figsize=graph_size)
   
    plt.scatter(plot_components, plot_generations, marker = "|", s = 150)
    
    plt.ylabel(r"Generation")
    plt.xlabel(r"Eigenvector components $rm{m}$")
    
    fig.savefig(os.path.dirname(__file__) + os.sep + os.pardir + os.sep + os.pardir + os.sep + os.pardir + os.sep + "images" + os.sep + "Eigenvector Fractal (" + method_name + ").png")

    
##plot_multiple_generations(spring_constant, mass_default, mass_ratio, generations, frequencies, plot_omega_phi)
##
##print "Finished omega vs. phi plot"
##
##plot_multiple_generations(spring_constant, mass_default, mass_ratios, generations, sp.pi / 2, plot_omega_mass)
##
##print "Finished omega vs. mass plot"

plot_fractal(spring_constant, mass_default, mass_ratio, sp.arange(1,101), sp.pi / 2)

plt.show()
