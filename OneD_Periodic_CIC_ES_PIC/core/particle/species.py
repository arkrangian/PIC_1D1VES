import numpy as np
from scipy import constants

class Species:

    def __init__(self, info: dict):
        # Initialization
        self.name = info['name']                                        # Particle name
        self.mass = info['mass']                                        # times of me (2me, 3me -> 2, 3)
        self.charge = info['charge']                                    # times of elementary charge (3e, -e -> 3, -1)
        self.temperature = info['temperature']                          # ev
        self.Sptcl_number_density = info['Sptcl_number_density']        # Density of superparticle (weight is not included)
        self.Sptcl_weight = info['Sptcl_weight']                        # Weight of superparticle
        self.position_init = info['position_init']                      # Position init ('random', 'uniform'), 'uniform' means uniformly distributed
        self.velocity_init = info['velocity_init']                      # Velocity init ('normal', 'injection') todo!!! (injection is hardcoded)
        self.motion_fix = info['motion_fix']                            # If True, don't push particle
        self.boundary_condition = info['boundary_condition']            # ('periodic'), 'reflect' will be added(Todo)
        self.num_of_particle = None                                     # Insert by particle generator

        # Simultion Inform
        self.label: np.ndarray = None                                   # Lagrangian label for single particle tracking
        self.location: np.ndarray = None                                # Location of particles,                            SI
        self.velocity: np.ndarray = None                                # Velocity of particles,                            SI
        self.Efield_at_particle: np.ndarray = None                      # ElectricField at each particle,                   SI

    def pushParticle(self, timeStep: np.float64, func_boundary_check: callable):
        # If motion_fix, don't move
        if (self.motion_fix):
            return
        # Push
        ratio_factor = (constants.elementary_charge/constants.electron_mass)
        u_velocity = (self.charge/self.mass) * self.Efield_at_particle * timeStep * ratio_factor        # du/dt = (charge/mass) * E
        gamma = np.sqrt(1+(u_velocity/constants.speed_of_light)**2)                                     # r(gamma) = sqrt(1 + (u/c)^2)

        self.velocity += u_velocity/gamma                                                               # v = u / r(gamma)
        self.location += self.velocity * timeStep                                                       # dx/dt = v

        # Sort
        self.location, self.velocity, self.label = self.__argsort(self.location, self.velocity, self.label)

        # Boundary Check(boundary_check function should get location_array)
        func_boundary_check(self.location)

        # ReSort
        self.location, self.velocity, self.label = self.__argsort(self.location, self.velocity, self.label)

    def singleParticleObserv(self, index):
        rag_index = np.where(self.label == index)[0]
        print("loc: ", self.location[rag_index])
        print("vel: ", self.velocity[rag_index])

    def get_totalKE(self) -> np.float64:
        return 0.5 * self.mass * constants.electron_mass * np.sum(self.velocity * self.velocity)
    
    def get_velDistribution(self, vel_width) -> tuple[np.ndarray, np.ndarray]:
        hist, edges = np.histogram(self.velocity, bins=np.arange(0, max(self.velocity)+vel_width, vel_width), density=True)

        # Filter Nonzero
        non_zero_bins = hist > 0
        hist = hist[non_zero_bins]
        edges = edges[:-1][non_zero_bins]

        return (hist, edges)

    def __argsort(self, r: np.ndarray, v: np.ndarray, l: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        idx = np.argsort(r)
        r = r[idx]
        v = v[idx]
        l = l[idx]
        return (r,v,l)