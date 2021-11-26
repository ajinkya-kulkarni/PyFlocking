# Program for active matter
# Code author: Ajinkya Kulkarni <kulkajinkya@gmail.com>

import numpy, os, math, pandas
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

os.system('clear')

os.system('rm timestep_*.csv')

os.system('rm image_*.png')

### Initialize the parameters

# number of particles
N = 1000

# noise intensity
eta = 0.1

# neighbor radius
rad_influence = 1

# time step
deltat = 0.01

# Maximum time
Tstart = 0.
Tend = 5

# Speed of the particle
Pspeed = 1.

#Boundaries of the domain
Xmin = 0.
Xmax = 1.
Ymin = 0.
Ymax = 1.

#Boundary condition - either 'periodic' or 'confined'
bcond = 'periodic'

#define angles
pi = numpy.pi
twopi = 2. * numpy.pi

### Required functions

#define the class called particles
class ActiveParticle(object):
    x = 0.
    y = 0.
    vx = 0.
    vy = 0.
    theta = 0.

    # The class "constructor" - It's actually an initializer
    def __init__(self, x, y, vx, vy, theta):
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.theta = theta


def make_particle(x, y, vx, vy, theta):
    actparticle = ActiveParticle(x, y, vx, vy, theta)
    return actparticle


# Function to write data to filename
def write_data(particles, filename):
    import pandas
    '''
    function to write the data to a file name
    '''
    datarray = numpy.zeros((N, 5))
    for i in range(N):
        p = particles[i]
        datarray[i,0] = p.x
        datarray[i,1] = p.y
        datarray[i,2] = p.vx
        datarray[i,3] = p.vy
        datarray[i,4] = p.theta

    df = pandas.DataFrame(datarray)
    df.columns = (['x', 'y', 'vx', 'vy', 'theta'])
    df.to_csv(filename)
    return()

# Functions to get neighbors
def get_neighbors(ap, particles, r):
    '''
    identify list of neighbours for a given particle
    '''
    neighbors = []

    for p in particles:
        if (p != ap):
            dist = numpy.sqrt((ap.x - p.x)**2 + (ap.y - p.y)**2)

            if (dist <r):
                neighbors.append(p)

    return (neighbors)


# Function to get average theta
def get_average(neighbors):
    '''
    return the average direction vector of the neigbhors
    input: list of neighbors
    '''

    n_neighbors = len(neighbors)

    avg_theta = 0.

    if (n_neighbors > 0):
        for nei in neighbors:
            avg_theta += nei.theta
        avg_theta = avg_theta / n_neighbors

    return (avg_theta)


def time_update(particles):
    '''
    function to update particle positions in time
    input: list of particles
    '''
    for ap in particles:
        #find neighbours for the particles
        ngbh = get_neighbors(ap, particles, rad_influence)
        #compute average velocity from that of the neighbour
        avg_theta = get_average(ngbh)
        #add noise to this velocity
        randnoise = numpy.random.uniform(-pi, pi, size=1)
        #set the particle theta
        ap.theta = avg_theta + randnoise * eta
        #set the new velocity components
        ap.vx = Pspeed * math.cos(ap.theta)
        ap.vy = Pspeed * math.sin(ap.theta)

        ap.x = ap.x + ap.vx * deltat
        ap.y = ap.y + ap.vy * deltat

        #set boundary conditions
        ap = boundaryconditions(ap, bcond)

    return(particles)


def boundaryconditions(particle, btype):
    '''
    set boundary conditions for the particle
    here, periodic boundary conditions are set
    '''
    if (btype=='periodic'):

        if (particle.x < Xmin):
            particle.x = Xmax + particle.x

        if (particle.x > Xmax):
            particle.x = particle.x - Xmax

        if (particle.y < Ymin):
            particle.y = Ymax + particle.y

        if (particle.y > Ymax):
            particle.y = particle.y - Ymax

    elif (btype=='confined'):

        if (particle.x <= Xmin):
            particle.vx = - particle.vx
            particle.x = Xmin

        if (particle.x >= Xmax):
            particle.vx = - particle.vx
            particle.x = Xmax

        if (particle.y <= Ymin):
            particle.vy = - particle.vy
            particle.y = Ymin

        if (particle.y >= Ymax):
            particle.vy = - particle.vy
            particle.y = Ymax

        particle.theta = numpy.arctan2(particle.vy, particle.vx)

    else:
        print ('specify the boundary conditions fully, either periodic or confined')

    return (particle)

# Generate random particle coordinations
x = numpy.random.uniform(Xmin, Xmax, size=N)
y = numpy.random.uniform(Ymin, Ymax, size=N)
theta = numpy.random.uniform(-pi, pi, size=N)

particles = []
for i in range(N):
    a = make_particle(0., 0., 0., 0., 0.)
    a.x = x[i]
    a.y = y[i]
    a.theta = theta[i]
    a.vx = Pspeed * math.cos(a.theta)
    a.vy = Pspeed * math.sin(a.theta)

    particles.append(a)

nsteps = int((Tend - Tstart)/deltat)

for i in range(nsteps):
    print ("timestep is", i*deltat, "out of", Tend)
    particles = time_update(particles)
    filename = 'timestep_' + str(i).zfill(4) + '.csv'
    write_data(particles, filename)


#### Post-process

for i in range(nsteps):
    filename = 'timestep_' + str(i).zfill(4) + '.csv'
    imgfile = 'image_' + str(i).zfill(4) + '.png'
    titlename = r'$\eta =$' + str(eta)
    print (filename, imgfile)

    df = pandas.read_csv(filename)

    colors1 = numpy.arctan2(df['vx'], df['vy'])

    colors2 = numpy.sqrt(df['vx']*df['vx'] + df['vy']*df['vy'])

    norm = Normalize()
    colormap = cm.jet
    norm.autoscale(colors1)
    color=colormap(norm(colors1))

    plt.xlim(Xmin, Xmax)
    plt.ylim(Ymin, Ymax)

    plt.close()
    plt.quiver(df['x'], df['y'], df['vx'], df['vy'],color=colormap(norm(colors1)))

    clb=plt.colorbar()
    clb.ax.set_title(r'$\theta$')

    # plt.quiver(df['x'], df['y'], df['vx'], df['vy'])
    plt.title(titlename)
    plt.savefig(imgfile, dpi = 300)

#Stitch video

os.system('rm video.mp4')
os.system('rm timestep_*.csv')

cmd = 'cat *.png | ffmpeg -f image2pipe -i - video.mp4'
print ('stitch video - command is', cmd)
os.system(cmd)

os.system('rm image_*.png')
