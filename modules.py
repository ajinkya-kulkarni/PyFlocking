# This file contains the functions used in the main program.

import numpy as np
import os
import pandas
from scipy.spatial import KDTree

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm

########################################################################################

class ActiveParticle(object):
	'''
	Class definition for an active particle.

	Attributes:
	x (float): the x-coordinate of the particle.
	y (float): the y-coordinate of the particle.
	vx (float): the x-velocity of the particle.
	vy (float): the y-velocity of the particle.
	theta (float): the direction of the particle.

	Methods:
	__init__(self, x, y, vx, vy, theta): the constructor for the ActiveParticle class.
	'''

	# Define the class attributes
	x = 0.
	y = 0.
	vx = 0.
	vy = 0.
	theta = 0.

	def __init__(self, x, y, vx, vy, theta):
		'''
		Initialize an ActiveParticle object.

		Args:
		x (float): the x-coordinate of the particle.
		y (float): the y-coordinate of the particle.
		vx (float): the x-velocity of the particle.
		vy (float): the y-velocity of the particle.
		theta (float): the direction of the particle.

		Returns:
		None.
		'''

		# Set the attributes of the particle
		self.x = x
		self.y = y
		self.vx = vx
		self.vy = vy
		self.theta = theta

########################################################################################

def make_particle(x, y, vx, vy, theta):
	'''
	Create a new ActiveParticle object.

	Args:
	x (float): the x-coordinate of the particle.
	y (float): the y-coordinate of the particle.
	vx (float): the x-velocity of the particle.
	vy (float): the y-velocity of the particle.
	theta (float): the direction of the particle.

	Returns:
	actparticle (ActiveParticle): a new ActiveParticle object.
	'''

	# Create a new ActiveParticle object using the specified parameters
	actparticle = ActiveParticle(x, y, vx, vy, theta)

	# Return the new particle object
	return actparticle

########################################################################################

def write_data(particles, filename, N):
	'''
	Write particle data to a file.

	Args:
	particles (list): a list of Particle objects to write to the file.
	filename (str): the name of the file to write to.
	N (int): the number of particles to write to the file.

	Returns:
	None.
	'''

	# Create an empty numpy array to hold the particle data
	datarray = np.zeros((N, 5))

	# Loop through each particle in the list and store its attributes in the numpy array
	for i in range(N):
		p = particles[i]
		datarray[i,0] = p.x
		datarray[i,1] = p.y
		datarray[i,2] = p.vx
		datarray[i,3] = p.vy
		datarray[i,4] = p.theta

	# Convert the numpy array to a pandas dataframe
	df = pandas.DataFrame(datarray)

	# Set the column names for the dataframe
	df.columns = (['x', 'y', 'vx', 'vy', 'theta'])

	# Write the dataframe to the specified file
	df.to_csv(filename)

	# Return None
	return None

########################################################################################

def get_neighbors(ap, particles, r):
	'''
	Returns a list of particles that are within a distance r of a given target particle.

	Parameters
	----------
	ap : object
		A particle object that represents the target particle.
		Must have attributes `x` and `y`.

	particles : list of objects
		A list of particle objects that represents all particles in the system.
		Each object must have attributes `x` and `y`.

	r : float
		The distance threshold for identifying neighboring particles.

	Returns
	-------
	list of objects
		A list of particle objects that are within a distance `r` of the target particle `ap`.
		Excludes the target particle `ap` itself.	
	'''
	neighbors = []

	for p in particles:
		if (p != ap):
			dist = np.sqrt((ap.x - p.x)**2 + (ap.y - p.y)**2)

			if (dist <r):
				neighbors.append(p)

	return (neighbors)

########################################################################################

def get_average(neighbors):
	'''
	return the average direction vector of the neigbhors
	input: list of neighbors
	'''

	n_neighbors = len(neighbors)

	avg_theta = 0.

	if (n_neighbors > 0):
		for single_neighbor in neighbors:
			avg_theta += single_neighbor.theta
		avg_theta = avg_theta / n_neighbors

	return (avg_theta)

########################################################################################

def time_update(particles, rad_influence, eta, Pspeed, deltat, bcond, Xmin, Xmax, Ymin, Ymax):
	'''
	Function to update particle positions in time.

	Args:
	particles (list): a list of Particle objects with x, y, vx, vy, and theta attributes.
	rad_influence (float): the radius of influence for finding neighbors.
	eta (float): the amount of noise to be added to the average velocity.
	Pspeed (float): the speed of the particles.
	deltat (float): the time step for the update.
	bcond (str): the boundary condition for the particle movement ('periodic' or 'reflective').
	Xmin (float): the minimum value for the x coordinate.
	Xmax (float): the maximum value for the x coordinate.
	Ymin (float): the minimum value for the y coordinate.
	Ymax (float): the maximum value for the y coordinate.

	Returns:
	particles (list): the updated list of Particle objects.
	'''

	# Loop through each particle in the list
	for single_particle in particles:
		
		# Find neighbors of the particle within the radius of influence
		all_neighbors = get_neighbors(single_particle, particles, rad_influence)
		
		# Calculate the average velocity of the neighbors
		avg_theta = get_average(all_neighbors)
		
		# Add random noise to the average velocity
		randnoise = np.random.uniform(-np.pi, np.pi, size=1)
		
		# Set the new particle direction (theta)
		single_particle.theta = avg_theta + randnoise * eta
		
		# Calculate the new velocity components (vx and vy) using the speed (Pspeed) and direction (theta)
		single_particle.vx = Pspeed * np.cos(single_particle.theta)
		single_particle.vy = Pspeed * np.sin(single_particle.theta)

		# Update the particle position using the velocity components and the time step (deltat)
		single_particle.x = single_particle.x + single_particle.vx * deltat
		single_particle.y = single_particle.y + single_particle.vy * deltat

		# Set boundary conditions using the boundaryconditions function
		single_particle = boundaryconditions(single_particle, bcond, Xmin, Xmax, Ymin, Ymax)

	# Return the updated list of particles
	return particles


########################################################################################

def boundaryconditions(particle, btype, Xmin, Xmax, Ymin, Ymax):
	"""
	Set boundary conditions for the particle.
	
	Parameters
	----------
	particle : ActiveParticle
		The particle to apply the boundary conditions to.
	btype : str
		The boundary type. Can be 'periodic' or 'confined'.
	Xmin : float
		The minimum value of the x-axis.
	Xmax : float
		The maximum value of the x-axis.
	Ymin : float
		The minimum value of the y-axis.
	Ymax : float
		The maximum value of the y-axis.
	
	Returns
	-------
	particle : ActiveParticle
		The particle with updated position and velocity based on the boundary conditions.
	"""
	if btype == 'periodic':
		# Periodic boundary conditions
		
		if particle.x < Xmin:
			particle.x = Xmax + particle.x

		if particle.x > Xmax:
			particle.x = particle.x - Xmax

		if particle.y < Ymin:
			particle.y = Ymax + particle.y

		if particle.y > Ymax:
			particle.y = particle.y - Ymax

	elif btype == 'confined':
		# Confined boundary conditions
		
		if particle.x <= Xmin:
			particle.vx = -particle.vx
			particle.x = Xmin

		if particle.x >= Xmax:
			particle.vx = -particle.vx
			particle.x = Xmax

		if particle.y <= Ymin:
			particle.vy = -particle.vy
			particle.y = Ymin

		if particle.y >= Ymax:
			particle.vy = -particle.vy
			particle.y = Ymax

		particle.theta = np.arctan2(particle.vy, particle.vx)

	else:
		# Invalid boundary type
		print('Error: Invalid boundary type specified. Must be "periodic" or "confined".')
	
	return particle

########################################################################################

def plot_particles(df, Xmin, Xmax, Ymin, Ymax, titlename, imgfile, rad_influence):
	"""
	Plots the particles in the dataframe as black circles with a white interior,
	with the velocity vectors represented by front arrows of varying color.
	
	Parameters
	----------
	df : pandas DataFrame
		A dataframe containing the particle positions and velocities.
	Xmin : float
		The minimum value of the x-axis.
	Xmax : float
		The maximum value of the x-axis.
	Ymin : float
		The minimum value of the y-axis.
	Ymax : float
		The maximum value of the y-axis.
	titlename : str
		The title of the plot.
	imgfile : str
		The filename to save the plot as.
	
	Returns
	-------
	None
	"""
	
	# Calculate angles and magnitudes
	angles = np.arctan2(df['vx'], df['vy'])
	magnitudes = np.sqrt(df['vx']**2 + df['vy']**2)
	
	# Set up colormap and normalization
	norm = Normalize()
	colormap = cm.hsv
	norm.autoscale(angles)
	colors = colormap(norm(angles))
	
	# Plot circles and arrows
	fig, ax = plt.subplots(figsize = (5, 5), constrained_layout = True)
	
	for x, y, angle, magnitude, color in zip(df['x'], df['y'], angles, magnitudes, colors):
		circle_size = 0.02
		circle = plt.Circle((x, y), circle_size, facecolor='tab:orange', edgecolor=None, alpha = 0.8, lw=0.8, zorder=1)
		ax.add_artist(circle)

		# ax.arrow(x, y, 0.05 * magnitude/2 * np.sin(angle), 0.05 * magnitude/2 * np.cos(angle), color=color, width=0.002, head_width=0.01, length_includes_head=True, zorder=2)
		ax.arrow(x, y, 2*circle_size* np.sin(angle), 2*circle_size* np.cos(angle), color='black', alpha = 0.8, width=1e-4, head_width=1e-2, length_includes_head=True, zorder=2)

	# Set plot limits and title
	ax.set_xlim(Xmin, Xmax)
	ax.set_ylim(Ymin, Ymax)
	# ax.set_title(titlename, pad = 10)
	ax.set_aspect('equal')
	ax.set_xticks([])
	ax.set_yticks([])

	# Save plot
	fig.savefig(imgfile, dpi=100)

	plt.close()

########################################################################################