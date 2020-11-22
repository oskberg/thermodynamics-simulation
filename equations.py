""" Module containing funcions used throughout the simulation

@ Oskar Hoegberg 13/02/2020

"""

import numpy as np

def boltzmann_distributution(x,a,c):
    return c * x * np.exp(-x**2/a)

def line(x,a,b):
    return a*x+b

def vdw_line(T,m):
    return m * T