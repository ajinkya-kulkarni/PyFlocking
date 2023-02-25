import numpy as np
import os
import pandas
import imageio
import shutil
from tqdm import tqdm

import sys
# Don't generate the __pycache__ folder locally
sys.dont_write_bytecode = True 
# Print exception without the buit-in python warning
sys.tracebacklimit = 0 

########################################################################################

from modules import *
from parameters import *

os.system('clear')

########################################################################################

# Remove old animation file if present
animation_file = 'animation.gif'
if os.path.exists(animation_file):
	os.remove(animation_file)

# Create a temporary directory for storing png and csv files
temp_dir = 'tmp'
if os.path.exists(temp_dir):
	shutil.rmtree(temp_dir)
os.mkdir(temp_dir)

########################################################################################

# Generate random particle coordinations
x = np.random.uniform(Xmin, Xmax, size=N)
y = np.random.uniform(Ymin, Ymax, size=N)
theta = np.random.uniform(-np.pi, np.pi, size=N)

########################################################################################

particles = []
for i in range(N):
	single_particle = make_particle(0, 0, 0, 0, 0)
	single_particle.x = x[i]
	single_particle.y = y[i]
	single_particle.theta = theta[i]
	single_particle.vx = Pspeed * np.cos(single_particle.theta)
	single_particle.vy = Pspeed * np.sin(single_particle.theta)

	particles.append(single_particle)

########################################################################################

nsteps = int((Tend - Tstart)/deltat)

########################################################################################

print()

for i in tqdm(range(nsteps), desc='Simulating particles'):
	particles = time_update(particles, rad_influence, eta, Pspeed, deltat, bcond, Xmin, Xmax, Ymin, Ymax)
	filename = os.path.join(temp_dir, f'timestep_{i:0>{len(str(nsteps))}}.csv')
	with open(filename, 'w') as f:
		write_data(particles, f, N)

	imgfile = os.path.join(temp_dir, f'image_{i:0>{len(str(nsteps))}}.png')
	titlename = fr'$\eta = {eta}$'
	filename = os.path.join(temp_dir, f'timestep_{i:0>{len(str(nsteps))}}.csv')
	with open(filename, 'r') as f:
		df = pandas.read_csv(f)

	# Call plot_particles function
	fig = plot_particles(df, Xmin, Xmax, Ymin, Ymax, titlename, imgfile, rad_influence)

########################################################################################

# Get a list of all PNG files in the temporary directory
png_files = [os.path.join(temp_dir, f) for f in os.listdir(temp_dir) if f.endswith('.png')]

# Sort the files by filename
png_files.sort()

# Create a new list to store the images
images = []

# Loop through each file and add it to the list of images
for i in tqdm(range(len(png_files)), desc='Generating video'):
	filename = png_files[i]
	with open(filename, 'rb') as f:
		image = plt.imread(f)
		image = (image * 255).astype(np.uint8)
		images.append(image)

# Create the GIF animation from the images
imageio.mimsave('animation.gif', images, fps=10)

# Delete the temporary directory and its contents
shutil.rmtree(temp_dir)

print()

########################################################################################