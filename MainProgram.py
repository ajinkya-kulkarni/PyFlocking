import numpy
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
os.system('rm animation.gif')

########################################################################################

# Create a temporary directory for storing png and csv files
temp_dir = 'tmp'
if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)
os.mkdir(temp_dir)

########################################################################################

# Generate random particle coordinations
x = numpy.random.uniform(Xmin, Xmax, size=N)
y = numpy.random.uniform(Ymin, Ymax, size=N)
theta = numpy.random.uniform(-numpy.pi, numpy.pi, size=N)

########################################################################################

particles = []
for i in range(N):
	a = make_particle(0., 0., 0., 0., 0.)
	a.x = x[i]
	a.y = y[i]
	a.theta = theta[i]
	a.vx = Pspeed * np.cos(a.theta)
	a.vy = Pspeed * np.sin(a.theta)

	particles.append(a)

########################################################################################

nsteps = int((Tend - Tstart)/deltat)

########################################################################################

print()

for i in tqdm(range(nsteps), desc = 'Simulating particles'):
	particles = time_update(particles, rad_influence, eta, Pspeed, deltat, bcond, Xmin, Xmax, Ymin, Ymax)
	filename = os.path.join(temp_dir, 'timestep_' + str(i).zfill(5) + '.csv')
	write_data(particles, filename, N)

	filename = os.path.join(temp_dir, 'timestep_' + str(i).zfill(5) + '.csv')
	imgfile = os.path.join(temp_dir, 'image_' + str(i).zfill(5) + '.png')
	titlename = r'$\eta =$' + str(eta)

	df = pandas.read_csv(filename)

	# Call plot_particles function
	fig = plot_particles(df, Xmin, Xmax, Ymin, Ymax, titlename, imgfile)

########################################################################################

# Get a list of all PNG files in the temporary directory
png_files = [os.path.join(temp_dir, f) for f in os.listdir(temp_dir) if f.endswith('.png')]

# Sort the files by filename
png_files.sort()

# Create a new list to store the images
images = []

# Loop through each file and add it to the list of images
for i in tqdm(range(len(png_files)), desc = 'Generating video'):
	image = plt.imread(png_files[i])
	image = (image * 255).astype(numpy.uint8)
	images.append(image)

# Create the GIF animation from the images
imageio.mimsave('animation.gif', images, fps=10)

# Delete the temporary directory and its contents
shutil.rmtree(temp_dir)

print()

########################################################################################