import numpy as np
from .species import Species
from scipy import constants
from ..world.boundary import Boundary

class Particle_Set:

    def __init__(self, infos: tuple[dict]):
        self.particles: list[Species] = []
        for info in infos:
            self.particles.append(Species(info))

    def push_particles(self, timestep: np.float64, func_boundary_check: callable):
        for particle in self.particles:
            particle.pushParticle(timeStep=timestep, func_boundary_check=func_boundary_check)

    def generate_particles(self, boundary: Boundary):
        for particle in self.particles:
            particle.num_of_particle = int(particle.Sptcl_number_density * boundary.system_volume)
            print("num of particles(",particle.name," ) :",particle.num_of_particle)
            self.__generateLabel(particle, boundary)
            self.__generateLocation(particle, boundary)
            self.__generateVelocity(particle, boundary)

    def get_totalKE(self) -> np.float64:
        totalKE = 0
        for particle in self.particles:
            totalKE += particle.get_totalKE()
        return totalKE
    
    def get_particle_velDist(self, speciesName: str, vel_width: np.float64 ) -> tuple[np.ndarray, np.ndarray]:
        for particle in self.particles:
            if(particle.name == speciesName):
                return particle.get_velDistribution(vel_width)
        return tuple[np.array([]), np.array([])]

    
    def __generateLabel(self, particle: Species, boundary: Boundary):
        particle.label = np.arange(0, particle.num_of_particle)

    def __generateLocation(self, particle: Species, boundary: Boundary):
        if particle.position_init == 'uniform':
            particle.location = np.array([(boundary.xMax/particle.num_of_particle)*i for i in range(particle.num_of_particle)])
        elif particle.position_init == 'random':
            particle.location = np.sort(boundary.xMax * np.random.rand(particle.num_of_particle))

    # Non-relativistic
    def __generateVelocity(self, particle: Species, boundary: Boundary):
        stdev = np.sqrt((constants.elementary_charge * particle.temperature) / (particle.mass * constants.electron_mass))
        particle.velocity = np.random.normal(0, stdev, particle.num_of_particle)
    
