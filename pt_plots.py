"""Module for plotting PT Diagrams of the simulation

This file plots pressure as a funtion of temperature for 
different ball radii in order to investigate the effect 
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


def pt_plot(
    ball_r,
    temp_interval = np.linspace(100, 500, 7),
    container_r = 20.0,
    ball_number = 64,
    repeats = 5,
    frames = 2_000,
):
    """Plots pressure against temperature

    Creates simulations with different temperature but keeps 
    everything else constant. It then records the pressure at 
    each temperature and plots it in a graph. For it to be 
    able to run multiple times with different ball radii, 
    plt.plot() mmust be called AFTER the function. Error 
    calulations are done by running the simulation ``repeats`` 
    number of times.

    The default values are set up for a simulation but the ball 
    radius needs to be specified.

    Parameters
    ----------
    ball_r : float
        The radius of the balls
    temp_interval : numpy array
        Array of temperatures to use in the simulation
    container_r : float
        Radius of container
    ball_number : int
        Number of balls
    repeats : int
        Number of times each simulation runs. Higher number 
        means more accurate simulation.
    frames : int
        Number of frames each simulation is run
    
    Returns
    -------
    b : float
        Returns the b coeffiecient in the van der Waals equation
    """

    # Create a dictionary to hold all statistics needed for a plot
    stats = {
        'T':[],
        'P':[],
        'P_std':[],
        'P_pred':[],
    }

    # Collect pressure data for each temperature given in 
    # temp_interval
    for temp in temp_interval:

        # Create simulation with each tempereature desired, 
        # keeping all other variables constant 
        PT_simulation = Simulation(
            container_radius=container_r, 
            balls=ball_number, 
            ball_radius=ball_r,
            ball_mass=1, 
            temperature=temp,
            )

        # Array to hold the pressure data from each experiment
        pressure_arr = np.zeros(repeats)


        for i in range(repeats):
            
            # Run for at least 2_000 frames for reliable 
            # pressure results when using many balls
            new_sim_stats = PT_simulation.run(
                frames, 
                stats=True, 
            )

            # Add the last pressure reading to the data array
            pressure_arr[i] = new_sim_stats['p_tot'][-1]
                                                                      
        # Add temperature, pressure, and standard deviation of 
        # pressure to dictionary
        stats['T'].append(temp)
        stats['P'].append(np.average(pressure_arr))
        stats['P_std'].append(np.std(pressure_arr))

    # Extract values from dictionary
    T = stats['T']
    P = stats['P']
    P_std = stats['P_std']

    # Create label for specific radius
    label_text = 'Radius ' + str(ball_r)
    plt.errorbar(
        T, P, fmt='none',yerr=P_std, 
        capsize=2, c='black')

    # Fit a curve using the van der Waals equation to the 
    # pressure data
    opt1, cov1 = curve_fit(eq.vdw_line, T, P)

    # Plot the van der Waals line
    x = np.linspace(min(T), max(T), 100)
    plt.plot(
        x,eq.vdw_line(x, opt1[0]), 
        label=label_text, linestyle='--')
    
    # Calculate b coefficient 
    container_area = container_r ** 2 * np.pi
    b = container_area / ball_number - 1/opt1[0]

    # Plot ideal gas law reference
    plt.plot(
        x,eq.vdw_line(x, ball_number / container_area),
        linestyle='-', linewidth=0.8,c='#4169E1')
    
    return b

# Container for b values
b_values = []

# Run the simulation for different ball radii to see the effect 
# it has on the b coefficient. The values are the ones used in 
# the lab report. For only one line, use only one r value.
rs = [0.5,1.0,1.3]
for r in rs:
    b_val = pt_plot(
        ball_r=r,
        temp_interval = np.linspace(100, 500, 7),
        container_r = 20.0,
        ball_number = 64,
        repeats = 5,
        frames = 2_000,
        )
    
    b_values.append(b_val)


for i, r in enumerate(rs):
    print('Ball radius ', r, '\t-\t', b_values[i])


# Plot the final result

plt.title('PT diagram')
plt.xlabel('Temperature (K$\cdot k_b$)')
plt.ylabel('Pressure (Pa/$k_b$)')
plt.legend()
plt.show()

