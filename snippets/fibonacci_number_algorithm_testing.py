import scipy as sp
import math

from timeit import timeit

def fib_(n=100000):
    first = 1
    second = 1

    for i in xrange(n):
        first, second = second, first + second

    return first

def fib_0(n=100000):
    """
    This is very cryptic and I don't understand it. Fast though.
    """
    
    v1, v2, v3 = 1, 1, 0    # initialise a matrix [[1,1],[1,0]]
    for rec in bin(n)[3:]:  # perform fast exponentiation of the matrix (quickly raise it to the nth power)
        calc = v2*v2
        v1, v2, v3 = v1*v1+calc, (v1+v3)*v2, calc+v3*v3
        if rec=='1':    v1, v2, v3 = v1+v2, v1, v2
    return v2

def fib_1(n=100000):
    fib_matrix = sp.matrix([[1, 1], [1, 0]], dtype=long)

    return (fib_matrix ** n)[1, 0]

def fib_2(n=100000):
    sqrt_5 = math.sqrt(5)
    first_term = ((1 + sqrt_5) / 2) ** n
    second_term = ((1 - sqrt_5) / 2) ** n

    return (first_term - second_term) / sqrt_5

def test(n):

    t  = timeit(fib_,  number=n) / float(n)
    t0 = timeit(fib_0, number=n) / float(n)
    t1 = timeit(fib_1, number=n) / float(n)

    try:
        t2 = timeit(fib_2, number=n) / float(n)
    except:
        print "Bad things are happening"

        print "AVERAGE TIME:"
        print "    ITERATIVE", t
        print "    SCIPY_MAT", t1
        print "    CRYPTIC_M", t0

    print "SPEEDS COMPARED TO DIAGONALISING:"
    print "    ITERATIVE", t / t2
    print "    SCIPY_MAT", t1 / t2
    print "    CRYPTIC_M", t0 / t2

from time import time

start = time()
for i in xrange(1, 1475):
    fib_2(i)

print time() - start, "s"
