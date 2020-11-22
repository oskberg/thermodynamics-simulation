"""Simulation Classes for Project B: Thermodynamics Snookered.

This module contains all the classes needed for the project

@ Oskar Hoegberg 13/02/2020

"""

import numpy as np
import pylab as pl

class Ball:
    """Creates a two dimensional ball

    Balls are solid circular objects that move in the simulation.

    Attributes
    ----------
    is_container : bool
        Defines if the ball is a container or a ball
    
    """

    def __init__(
        self, mass=0.0, radius=0.0, velocity=[0.0, 0.0], 
        position=[0.0, 0.0], is_container=False):
        """ __init__ method of a ball

        Parameters
        ----------
        mass : float
            The mass of the ball
        radius : float
            The radius of the ball
        velocity : list
            The two dimensional velocity of the ball. Given as a list which 
            is converted into a numpy array
        position : list
            Position of the ball in two dimensional cartesian coordinates. 
            Given as list and get_s converted into numpy array
        is_container : bool
            Defines weather the ball is a container or not. `True` if ball 
            is a container, otherwise `False`

        """
        self._mass = mass
        self._radius = radius
        self._velocity = np.array(velocity)
        self._position = np.array(position)
        self.is_container = is_container

        # If self is a container create patch with only a border so it looks 
        # like a container.
        # Otherwise, create ball with fill colour.
        if is_container: 
            self._patch = pl.Circle(self.get_position(), self.get_radius(), 
            edgecolor='blue', facecolor='None')
        elif not is_container: 
            self._patch = pl.Circle(self.get_position(), 
            self.get_radius(), facecolor='red')

    def get_patch(self):
        return self._patch

    def set_mass(self, mass):
        self._mass = mass
        return self

    def get_mass(self):
        return self._mass

    def set_radius(self, radius):
        self._radius = radius
        return self

    def get_radius(self):
        return self._radius

    def set_velocity(self, velocity):

        self._velocity = np.array(velocity)
        return self

    def get_velocity(self):
        return self._velocity

    def set_position(self, position):
        self._position = np.array(position)
        self._patch.center = position
        return self

    def get_position(self):
        return self._position

    def move(self, dt):
        """Moves ball in time along velocity vector. Updates position of ball.

        Takes a time `dt` and calculates the distance it moves given its velocity. 
        This distance vector is added to the current position. 

        Parameters
        ----------
        dt : float
            Time to move the ball
        """
        self.set_position(self.get_position() + self.get_velocity() * dt) 

    def time_to_collision(self, other):
        """Calculates the time to collision between self and another ball. 

        Uses difference in postition and velocity between self and another ball to 
        calculate in what time the two will collide.

        Parameters
        ---------
        other : Ball object
            The other ball in the collision.

        Returns
        -------
        t_return : float
            The time until next collision.
                - Positive time if collision in the future
                - np.nan otherwise
        
        """

        # Define parameters of the equation
        v = self._velocity - other._velocity
        r = self._position - other._position
        
        # Self cannot collide with itself
        # Container shouldn't be considered
        # We do not want the container as the first object in collision function
        if (self.is_container) or (self is other):
            return np.nan

        # Ball - Container collision
        # There sould always be a positive and a negative solution
        # Return largest time
        elif other.is_container:
            R = abs(self._radius - other._radius)

            if not np.linalg.norm(r) > R or 1: 
                
                t1 = (-v.dot(r) 
                    + np.sqrt( v.dot(r) ** 2 - v.dot(v) * (r.dot(r) - R ** 2) )) / v.dot(v)
                t2 = (-v.dot(r) - np.sqrt( v.dot(r) ** 2 
                    - v.dot(v) * (r.dot(r) - R ** 2) )) / v.dot(v)

                return t1 if t1 > t2 else t2
        
        # Ball - Ball collision
        elif not other.is_container:
            # Check if moving towards each other
            alignment = r.dot(v)

            # Moving away so will not collide
            if alignment > 0.0:
                return np.nan
            
            # Approaching so return smallest value
            elif alignment < 0.0:
                R = self._radius + other._radius

                t1 = (-v.dot(r) 
                    + np.sqrt( v.dot(r) ** 2 - v.dot(v) * (r.dot(r) - R ** 2) )) / v.dot(v)
                t2 = (-v.dot(r) 
                    - np.sqrt( v.dot(r) ** 2 - v.dot(v) * (r.dot(r) - R ** 2) )) / v.dot(v)

                return t1 if t1 < t2 else t2

            if alignment == 0:
                print(self is other)


    def collide(self, other):
        """Updates velocities of two balls `self` and `other`. 

        Uses 2D collision mechanics to update the velocities of two different 
        ball objects. Uses current position, velocity and mass so has be brought
        together using the move function before colliding to get correct velocities.

        Parameters
        ---------
        other : Ball object
            The other ball in the collision.
                                                                                  
        Returns
        -------
        impulse_on_container : float
            Magnitude of impulse onto the container due to ball. 0 if collision between 
            two balls.
        
        """
        # Define vaiables in equation

        r1 = self.get_position()
        r2 = other.get_position()
        
        m1 = self.get_mass()
        m2 = other.get_mass()
        
        u1 = self.get_velocity()
        u2 = other.get_velocity()

        # Ball - Container collision 
        # Should not change velocity of container
        if other.is_container:
            v1 = u1 - 2 * (u1 - u2).dot(r1-r2) * (r1-r2) / (np.linalg.norm(r1-r2)**2)
            
            # Impulse is difference in momentum
            impulse_on_container = (v1 - u1) * m1

            self.set_velocity(v1)

            return impulse_on_container

        # Ball - Ball collision
        else:
            v1 = u1 - 2 * m2 / (m1 + m2) * (u1 - u2).dot(r1-r2) * (r1-r2) / (np.linalg.norm(r1-r2)**2)
            v2 = u2 - 2 * m1 / (m1 + m2) * (u2 - u1).dot(r2-r1) * (r2-r1) / (np.linalg.norm(r2-r1)**2)

            self.set_velocity(v1)
            other.set_velocity(v2)

            # Impulse on container is 0
            return 0

