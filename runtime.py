""" Simple module for running the simulation

Use this file to generate ball separation plots and 
boltzmann distributions.

@ Oskar Hoegberg 13/02/2020

"""

import numpy as np

from simulation import Simulation

# Creates a simulation
sim = Simulation(
    container_radius=20, balls=64, ball_radius=0.5, 
    ball_mass=1, temperature=300)

# Runs the simulation for 10_000 frames without animation
# Once completed prints some stats and plots separations
# and a Boltzmann distribution

#Set separations to `False` for much better performance
sim.run(
    10_000,
    animate=True, 
    stats=True,
    boltzmann=True, 
    separations=False,
    )
