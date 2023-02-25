# This file contains the simulation parameters used in the main program.

# The boundaries of the simulation domain of length L.
# This parameter defines the spatial extent of the simulation domain. The domain is a rectangle with lower left corner at (Xmin, Xmin + L) and upper right corner at (Xmax, Xmax + L).
L = 1
Xmin = 0
Xmax = Xmin + L
Ymin = 0
Ymax = Ymin + L

# The boundary condition used in the simulation, which can be either 'periodic' or 'confined'.
# This parameter controls the behavior of particles near the boundary of the simulation domain. If bcond is 'periodic', particles that move out of the domain on one side will re-enter the domain from the opposite side. If bcond is 'confined', particles that move out of the domain are reflected back into the domain.
bcond = 'periodic'

# The number of particles in the simulation.
N = 200

# The noise intensity, which determines how much random motion is added to each particle at each time step.
# This parameter controls the degree of stochasticity in the particle dynamics. Higher values of eta lead to more chaotic behavior.
# Range from 0 to 1
eta = 0.8

# The radius of influence of each particle, which determines which other particles it can interact with.
# This parameter controls the spatial range of interactions between particles. Larger values of rad_influence lead to longer-ranged interactions.
rad_influence = 0.2*L

# The time step used in the simulation.
# This parameter controls the granularity of the simulation time evolution. Smaller values of deltat lead to more accurate simulations, but also increase the computational cost.
deltat = 1e-2

# The starting and ending times of the simulation, respectively.
# This parameter controls the duration of the simulation. The simulation will run from Tstart to Tend.
Tstart = 0
Tend = 2

# The speed of each particle.
# This parameter controls the average magnitude of the particle velocities. Higher values of Pspeed lead to faster particle motion.
Pspeed = 1