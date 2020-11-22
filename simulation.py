""" Module for the Simulation class

@ Oskar Hoegberg 13/02/2020

"""

import csv
import math
import random as rd

import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
import scipy.optimize as spo

import equations as eq
from ball import Ball


class Simulation:
    """Class containing the functionality to run the simulation

    Simulation is run with a container with an arbitrary 
    number of balls in it. Run simulation by calling ``run``. 
    Different statistics can be turned on and off using 
    different keywords when calling ``run``.

    """
    
    def __init__(
        self, container_radius=10, balls=0, 
        ball_mass=1.0, ball_radius=1.0, temperature=300
    ):
        """ Initialise a simulation

        Parameters
        ----------
        container_radius : float
            The radius of the container
        balls : int
            Number of balls in the simulation
        ball_mass : float
            Mass of a single ball
        ball_radius : float
            Radius of a single ball
        temperature : float
            The desired temperature of the system
        
        """

        self._time_passed = 0.0
        self._container = Ball(
            radius=container_radius, is_container=True)
        self._number_of_balls = balls
        self._balls = []
        self._dt_list = np.full((balls + 1, balls + 1), np.nan)
        self._stat_dict = {
            'r_sep': [],
            'r_mag': [],
            'v_tot': [],
            'p_tot': [],
            'T_tot': [],
            'ke_tot': [],
            'ke_avg': [],
            'v_dist': [],
            'pt_tot': [],
            'time_passed': [],
        }

        # For conistent results a seed may be set for the exact 
        # same initial condition
        # rd.seed(100)
        
        # kb is set to 1 for units to be nicer to work with
        kb = 1.0

        # Position of balls evenly spaced
        ball_positions = np.linspace(
            -1/np.sqrt(2)*(container_radius - ball_radius - 0.1),
             1/np.sqrt(2)*(container_radius - ball_radius - 0.1), 
             math.ceil(np.sqrt(balls))
             )

        # Speed required for simulation to have desired temperature
        ball_speed = np.sqrt(2 * kb * temperature / ball_mass)
        max_v = ball_speed * 2

        # Create balls
        # Set ball position sin a grid with side length
        # sqrt(number of balls)
        for x in range(math.ceil(np.sqrt(balls))):
            for y in range(math.ceil(np.sqrt(balls))):
                # If desired number of balls is reached, stop creating balls
                if len(self._balls) >= balls: 
                    break

                # Random x velocity from -max_v to +max_v
                ball_v_x = (rd.random() - 0.5) * max_v

                # y veloctiy is set limited for temperature to be reached
                rand_neg = (-1) ** rd.randint(0,10)
                ball_v_y =  rand_neg *  np.sqrt(ball_speed ** 2 - ball_v_x ** 2)

                ball_v = np.array([ball_v_x, ball_v_y])

                ball_pos = np.array([ball_positions[x], ball_positions[y]])
                ball = Ball(
                    mass=ball_mass, radius=ball_radius, 
                    velocity=ball_v, position=ball_pos)
                
                # Add the ball to the ball array
                self._balls.append(ball)
            
            if len(self._balls) >= balls: 
                break

        # add container LAST to the list of balls
        self._balls.append(self._container)

        # Calculate all collision times for dt_list once
        for i, ball_self in enumerate(self._balls):
            for j, ball_other in enumerate(self._balls):
                if j > i:
                    t_return = ball_self.time_to_collision(ball_other)
                    self._dt_list[i][j] = t_return



    def next_collision(self, stats=False):
        """ Funciton which takes the simulation through a collision

        Looks for the smalles dt in in dt_list. Then move all balls
        and perform collisions. Last it recaclulates some time in
        dt_list to keep it correct.

        Parameters
        ---------
        stats : bool
            True if statistics are wanted
        
        """
        
        # Smallest non nan value in dt list is next collision time
        smallest_dt = np.nanmin(self._dt_list)

        # The indices in the 2d array represent what balls are colliding
        smallest_dt_index = np.where(self._dt_list == smallest_dt)
        smallest_dt_i = smallest_dt_index[0][0]
        smallest_dt_j = smallest_dt_index[1][0]
        
        # Get the two balls in the collision
        ball_i = self._balls[smallest_dt_i]
        ball_j = self._balls[smallest_dt_j]

        # Move all the balls by dt
        for ball in self._balls:
            ball.move(smallest_dt)

        # Update the ellapsed time of the simulation
        self._time_passed += smallest_dt

        # Keep all dt so times between any collisions may be calculated 
        self._stat_dict['time_passed'].append(smallest_dt)

        # Time remaining for collisions is decreased by the time the balls were moved by
        self._dt_list -= smallest_dt

        # Collide the two balls which will collide first and measure impulse
        impulse_on_container = ball_i.collide(ball_j)

        # Add contribution of the impulse to the aggregate pressure
        container_circumference  = self._container.get_radius()*2*np.pi
        self._stat_dict['pt_tot'].append(
            np.linalg.norm(impulse_on_container/container_circumference))

        # The two balls involved in the collisions require recalculated collision times
        # Calculate these and update 
        for l, ball_other in enumerate(self._balls):
            dt_i = ball_i.time_to_collision(ball_other)
            dt_j = ball_j.time_to_collision(ball_other)

            self._dt_list[smallest_dt_i][l] = dt_i
            self._dt_list[smallest_dt_j][l] = dt_j


    def run(
        self, num_frames, animate=False, stats=False, 
        boltzmann=False, separations=False):
        """Takes the simulation through a set number of frames

        Run function takes the simulation through a specified 
        number of frames by calling ``next_collision`` in a 
        for loop. It also calulates statistics every frame and
        adds it to the statistics dictionary. Setting parameters 
        to `True` adds additional funcionality. 

        Parameters
        ---------
        num_frames : int
            The number of frames to run the simulation for
        animate : bool
            Set to `True` to plot the position of the balls 
            every frame.
        stats : bool
            Set to `True` to print pressure and temperature 
            statistics at the end.
        boltzmann : bool
            Set to `True` to plot the speed distribution at 
            the end of the simulation
        separations : bool
            SLOW, USE ONLY WHEN NEEDED! Set to `True` to plot 
            the ball separation and distance from origin. 

        Returns
        -------
        stat_dict : dictionary of lists
            Returns all the statistics collected over 
            the simulation

        """

        # Animate the simulation by updating patches of 
        # all balls every frame
        if animate:
            # Add one to the limits of the drawing to have 
            # some space around container
            radius = self._container.get_radius() + 1
            f = pl.figure(figsize=(5,5))
            ax = pl.axes(xlim=(-radius, radius), ylim=(-radius, radius))
            ax.add_artist(self._container.get_patch())

            # Loop over the ball array and update patches 
            for ball in self._balls:
                ax.add_patch(ball.get_patch())
        
        # Run the simulation for num_frames number of frames. 
        # Every frame, collect statistics about the system 
        # Add statistics to the statistics dictionary
        for frame in range(num_frames):
            self.next_collision(stats=stats)
            
            # Get total kinetic enegy of the system
            ke = self.kinetic_energy()
            self._stat_dict['ke_tot'].append(ke)
            
            # Get average kinetic energy by dividing by number 
            # of balls
            ke_avg = ke/self._number_of_balls
            self._stat_dict['ke_avg'].append(ke_avg)
            
            # Add the average pressure over an interval specified 
            # in pressure function
            pressure = self.pressure()[0]
            self._stat_dict['p_tot'].append(pressure)
            
            # Get current temperature
            T = self.temperature()
            self._stat_dict['T_tot'].append(T)

            # Add all current speeds to statistics
            # Use `extend` since ``v_dist`` is a list
            v_dist = self.velocity_distribution()
            self._stat_dict['v_dist'].extend(v_dist)

            # Separations calculation is very intensive so should 
            # only be run if specifically desired. Separations 
            # means the distance between balls or distance to origin
            if separations:
                r_mag, r_sep = self.separations()
                self._stat_dict['r_mag'].extend(r_mag)
                self._stat_dict['r_sep'].extend(r_sep)

            # Print out some statistics every 100 frames so we know 
            # the simulation is progressing as expected
            if frame % 100 == 0 or frame == num_frames:
                print('Collision number', frame)
                print('P\t-\t', pressure)
                print('T\t-\t', T,'\n')

            if animate:
                pl.pause(0.001)
        
        if animate:
            pl.show()

        # For Maxwell-Boltzmann plot, set boltzman to 
        # true in run call
        if boltzmann:
            # kb = 1
            # Get variables needed to plot the graph
            v_dist = np.array(self._stat_dict['v_dist'])
            v_avg = np.average(v_dist)
            print('\nAverage of velocities\t-\t', v_avg)
            print('Variance of velocities\t-\t', np.var(v_dist))
            
            # Plot histogram of speeds
            no_bins = 25
            n, bins, patches = plt.hist(
                v_dist, bins=no_bins, 
                label='Speed Data', 
                color='#4169E1')

            # Bins are the edges of the bins but we need the centres
            bin_centres = [0.5 * (bins[i] + bins[i+1]) for i in range(len(bins) - 1)]

            # Curve fit the histogram with Maxwell-Boltzmann curve
            opt_hist, cov_hist = spo.curve_fit(
                eq.boltzmann_distributution, bin_centres, n)
            
            # Plot fitted curve
            v = np.linspace(0, max(v_dist), 100)
            plt.plot(v, eq.boltzmann_distributution(
                v, opt_hist[0], opt_hist[1]), 
                label='Maxwell-Boltzmann Distribution', 
                c='#FF4500')
            
            # Set plot parameters and show plot
            plt.title('Speed Distribution')
            plt.xlabel('Speed')
            plt.ylabel('Number of Balls')
            plt.legend()
            plt.show()

        # If stats are enabled, stats are calculated
        if stats:
            
            # The actual pressure is the pressure measured last. 
            # This is assumed to be the most correct value.
            p_true, p_std = self.pressure()
            T_tot = np.array(self._stat_dict['T_tot'])

            # Calculate the predicted pressure accoring to the 
            # inverse square law. Also calulate the predicted 
            # ratio between pressure and temperature. 
            area = self._container.get_radius()**2*np.pi
            p_predicted = 1/(area) * self._number_of_balls * np.average(T_tot)
            ratio_predicted = 1/(area) * self._number_of_balls * 1
            
            # Print out the stats
            print(
                '\n\nN -', self._number_of_balls,'; T -', np.average(T_tot),
                '; L -', self._container.get_radius()*2*np.pi)
            print('Predicted Preassure\t-\t', p_predicted)
            print('Actual Preassure\t-\t', p_true)
            print('Predicted ration p/T\t-\t', ratio_predicted)
            print('Actual ration p/T\t-\t', p_true/np.average(T_tot))

            # If separations is enabled, creates plots showing the 
            # separation between balls and balls. Also balls and container.
            if separations:
                
                # r_sep is separation between balls. r_mag is distance 
                # form origin for balls
                r_sep = np.array(self._stat_dict['r_sep'])
                r_mag = np.array(self._stat_dict['r_mag'])
                
                # Check if there are balls outside of container or if 
                # balls are overlapping
                ball_radius = self._balls[0].get_radius()

                # Allow for a little bit of overlapp due to floating 
                # point error. This does not mean balls are missbehaving
                tolerance = ball_radius * 0.01

                # Check how many balls are outside container by looking 
                # at the number of values in r_mag which are larger 
                # than the container
                balls_outside = len(np.where(
                    r_mag > self._container.get_radius() - ball_radius + tolerance)[0] )
                
                # Check how many balls are stuck toghether by similar method
                balls_too_close = len(np.where(r_sep < ball_radius * 2 - tolerance)[0])

                print(
                    '\n\n\tBalls Outside Container:', 
                    balls_outside)
                print(
                    '\tBalls Stuck Together:', 
                    balls_too_close)

                bins=100
                
                plt.subplot(1,2,1)
                plt.title('Ball Distance from Centre')
                plt.hist(
                    r_mag, bins=bins, 
                    label='Data', 
                    color='#4169E1')
                plt.xlabel('Distance')
                plt.ylabel('Number of balls')
                plt.legend()
                
                plt.subplot(1,2,2)
                plt.title('Distance Between Balls')
                plt.hist(
                    r_sep, 
                    bins=bins, 
                    label='Data', 
                    color='#4169E1')
                plt.xlabel('Distance')
                plt.legend()
                
                plt.show()

        return self._stat_dict


    def kinetic_energy(self):
        """Sums the kinetic energy of all the balls"""
        ke = sum([0.5 * ball.get_mass() * np.linalg.norm(ball.get_velocity()) ** 2 for ball in self._balls[:-1]])
        return ke

    def temperature(self):
        """Calculates the temperature from the kinetic energy"""
        kb = 1
        T = self.kinetic_energy() * 2/2 *1/kb /(self._number_of_balls)
        return T

    def pressure(self):
        """"Calcualtes the average pressure an interval

        In the beginning the pressure is irregular and not reliable. 
        Taking the average over the entire simulation therefore
        introduces some error. It is better then to calulate the
        average over the last x number of frames. Through testing,
         1000 frames turned out to be a good number producing 
        relatively smooth pressure.
        
        Returns
        -------
        pressure : float
            The average pressure on the container the last 
            1000 frames. 
        """
        
        # Interval determined by testing
        # 1000 frames is good for up to 100 balls
        interval = 1000
        # Check that there are more pressure readings than 
        # the interval requires
        if len(self._stat_dict['time_passed']) > interval:
            # Get the time of the interval
            time = sum(self._stat_dict['time_passed'][-interval:])
            P_time = sum(self._stat_dict['pt_tot'][-interval:])
            std = np.std(self._stat_dict['pt_tot'][-interval:])
            # Return the average pressure over the time
            return P_time/time, std/time
        else:
            # If there are not enough pressure readings, return 0 
            return 0, 0

    def velocity_distribution(self):
        """Returns a list of the speeds for all balls."""
        return [np.linalg.norm(ball.get_velocity()) for ball in self._balls[:-1]]

    def separations(self):
        """Calculates the separation between balls and 
        their separation from the origin.

        This is a very slow funtion so should only be called 
        if needed. The iteration implicitly creates a 2d list 
        of all balls and calulates the between ball index i 
        and ball index j. To avoid double counting, only half 
        of this list is evaluated and i == j is exluded as this 
        is distance to self.

        Example array: 
        # x x x     
        # # x x
        # # # x
        # # # #
        => Only the x positions need to be evaluated

        Returns
        -------
        r_mag : list
            A list containing the separations from the origin 
            for every ball
        r_sep : list
            A list containing the separations between all 
            the balls

        """
        
        # Exlude container (last in self._balls) 
        # from calculations
        balls_only = self._balls[:-1].copy()
        
        r_mag = []
        r_sep = []
        
        # Loop over the balls to measure distance from origin
        for i, ball1 in enumerate(balls_only):
            r_mag.append(np.linalg.norm(ball1.get_position()))

            # Loop again to measure teh balls distance to 
            # other balls
            for ball2 in balls_only[i+1:]:
                separation = np.linalg.norm(ball1.get_position() 
                    - ball2.get_position())

                r_sep.append(separation)
        
        return r_mag, r_sep
