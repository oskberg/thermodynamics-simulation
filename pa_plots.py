"""Module for plotting PA Diagrams of the simulation

This file plots pressure as a funtion of container area 
for different ball radii in order to investigate the effect
this has on the b coefficient in van der Waals equation.
                                                                      
It may take a long time to generate the same plots as are 
used in the lab report...

@ Oskar Hoegberg 13/02/2020

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

import equations as eq
from ball import Ball
from simulation import Simulation

# The values 64 and 300 are specific to cofiguration 64 balls
# and temperature 300. These are hardcoded due to the curve_fit
# function. They MUST be changed when changing the simulation.
def isotherm_vdw(A,b):
    return 64*300/(A-64*b)

def pa_plot(
    ball_r,
    temp = 300,
    num_in_interval = 7,
    ball_number = 64,
    repeats = 5,
    frames = 2_000,
):
    """Plots pressure against temperature

    Creates simulations with different temperature but keeps 
    everything else constant. It then records the pressure at 
    each temperature and plots it in a graph. For it to be able 
    to run multiple times with different ball radii, plt.plot() 
    mmust be called AFTER the function. Error calulations are 
    done by running the simulation ``repeats`` number of times.

    The default values are set up for a simulation but the ball 
    radius needs to be specified.

    Parameters
    ----------
    ball_r : float
        The radius of the balls
    temp : float
        Temperature to use in the simulation
    num_in_interval : int
        Number of data points to measure
    ball_number : int
        Number of balls
    repeats : int
        Number of times each simulation runs. Higher number means 
        more accurate simulation.
    frames : int
        Number of frames each simulation is run
    
    Returns
    -------
    b : float
        Returns the b coeffiecient in the van der Waals equation
    """

    # Create a dictionary to hold all statistics needed for a plot
    stats = {
        'A':[],
        'P':[],
        'P_std':[],
        'P_pred':[],
    }

    # 14 works well for 64 balls
    # replace 14 with the smallest radius required for the 
    # container at ball_r = 1
    container_r_interval = np.linspace(
        14 * ball_r, 14 * ball_r + 10, num_in_interval)

    # Collect pressure data for each temperature given in temp_interval
    for c_r in container_r_interval:

        # Create simulation with each tempereature desired, 
        # keeping all other variables constant 
        PT_simulation = Simulation(
            container_radius=c_r, 
            balls=ball_number, 
            ball_radius=ball_r,
            ball_mass=1, 
            temperature=temp,
            )

        # Array to hold the pressure data from each experiment
        pressure_arr = np.zeros(repeats)


        for i in range(repeats):
            
            # Run for at least 2_000 frames for reliable pressure 
            # results when using many balls
            new_sim_stats = PT_simulation.run(
                frames, 
                stats=True, 
            )

            # Add the last pressure reading to the data array
            pressure_arr[i] = new_sim_stats['p_tot'][-1]
        
        # Add temperature, pressure, and standard deviation of 
        # pressure to dictionary
        stats['A'].append(c_r ** 2 * np.pi)
        stats['P'].append(np.average(pressure_arr))
        stats['P_std'].append(np.std(pressure_arr))

    # Extract values from dictionary
    A = stats['A']
    P = stats['P']
    P_std = stats['P_std']

    # Create label for specific radius
    label_text = 'Radius ' + str(ball_r)
    plt.errorbar(A, P, fmt='none',yerr=P_std, capsize=2, c='b')

    # Fit a curve using the van der Waals equation to the pressure data
    opt1, cov1 = curve_fit(isotherm_vdw, A, P)

    # Plot the van der Waals line
    x = np.linspace(min(A), max(A), 100)
    plt.plot(
        x,isotherm_vdw(x, opt1[0]), label=label_text, 
        linestyle='--')
    
    # Calculate b coefficient 
    b = opt1[0]

    # Plot ideal gas law reference
    plt.plot(
        x,isotherm_vdw(x, 0), linestyle='-', 
        linewidth=0.8,c='#4169E1')
    
    return b

# Container for b values
b_values = []

# Run the simulation for different ball radii to see the effect 
# it has on the b coefficient. The values are the ones used in 
# the lab report. For only one line, use only one r value. Also
# the values in isotherm_vdw need to be changed.
rs = [0.5,1.0,1.3]
for r in rs:
    b_val = pa_plot(
        ball_r=r,
        temp = 300,
        num_in_interval = 5,
        ball_number = 64,
        repeats = 5,
        frames = 2_000,
        )

    b_values.append(b_val)
    
for i, r in enumerate(rs):
    print('Ball radius ', r, '\t-\t', b_values[i])

# Plot the final result

plt.title('PA diagram')
plt.xlabel('Area')
plt.ylabel('Pressure ')
plt.legend()
plt.show()

