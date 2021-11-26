
# Program for active matter - V1
# Code author: Ajinkya Kulkarni <kulkajinkya@gmail.com>

import numpy, os, math, pandas
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize


### Required functions

# In[ ]:

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


# In[ ]:

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


# In[ ]:

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


# In[ ]:

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


# In[ ]:

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
        ap = boundaryconditions(ap)

    return(particles)


# In[ ]:

def boundaryconditions(particle):
    '''
    set boundary conditions for the particle
    here, periodic boundary conditions are set
    '''

    if (particle.x < Xmin):
        particle.x = Xmax + particle.x


    if (particle.x > Xmax):
        particle.x = particle.x - Xmax

    if (particle.y < Ymin):
        particle.y = Ymax + particle.y

    if (particle.y > Ymax):
        particle.y = particle.y - Ymax

    return (particle)


### Initialize the parameters

# In[ ]:

# number of particles
N = 100

# noise intensity
eta = 0.1

# neighbor radius
rad_influence = 0.02

# time step
deltat = 0.01

# Maximum time
Tstart = 0.
Tend = 1.

# Speed of the particle
Pspeed = 1.

#Boundaries of the domain
Xmin = 0.
Xmax = 1.
Ymin = 0.
Ymax = 1.

#define angles
pi = numpy.pi
twopi = 2. * numpy.pi


# In[ ]:

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


# In[ ]:

write_data(particles, 'check.csv')


### Check data

# In[ ]:

df = pandas.read_csv('check.csv')


# In[ ]:

#plt.scatter(df['x'], df['y'])
#plt.scatter(df['vx'], df['vy'])
plt.quiver(df['x'], df['y'], df['vx'], df['vy'])
#plt.hist(df['theta'])
plt.show()


### Get neighbours of each particle

# In[ ]:

nsteps = int((Tend - Tstart)/deltat)
print nsteps


# In[ ]:

for i in range(nsteps):
    print "timestep is", i*deltat
    particles = time_update(particles)
    filename = 'timestep_' + str(i).zfill(2) + '.csv'
    write_data(particles, filename)



#### Post-process

# In[ ]:

for i in range(100):
    filename = 'timestep_' + str(i).zfill(2) + '.csv'
    imgfile = 'image_' + str(i).zfill(2) + '.png'
    titlename = 'eta = 0.1, timestep =' + str(i * deltat).zfill(3) + 'sec'
    print filename, imgfile

    df = pandas.read_csv(filename)

    colors = numpy.arctan2(df['vx'], df['vy'])
    norm = Normalize()
    colormap = cm.jet
    norm.autoscale(colors)
    color=colormap(norm(colors))

    plt.close()
    plt.quiver(df['x'], df['y'], df['vx'], df['vy'],color=colormap(norm(colors)))
    plt.title(titlename)
    plt.savefig(imgfile)


# In[3]:

#Stitch video
cmd = 'rm video.mp4'
os.system(cmd)

cmd = 'avconv -r 10 -i image_%4d.png -b:v 1000k video.mp4'
print 'stitch video - command is', cmd
os.system(cmd)
