import numpy as np
from ..particle.species import Species
from ..world.boundary import Boundary
from scipy import constants

# 1-dimensional Periodic Cloud-In-Cell(1st order) Interpolator
class Interpolater_CIC:

    """
    chargeDensityGrid
    """
    @staticmethod
    def chargeDensityInterpolate(particles: list[Species], boundary: Boundary) -> None:
        chargeDensityGrid = np.zeros(shape=boundary.num_of_cells, dtype=np.float64)

        for particle in particles:
            # Charge Check
            if particle.charge == 0 :
                continue
            
            # Variable Setting
            count = np.zeros(boundary.num_of_cells)         # Particle Number in Each Cell
            sRDC = np.zeros(boundary.num_of_cells)          # Sum of Relative Distance(Particle_Location - Left_Grid) in Cell
            relDist = np.zeros(particle.num_of_particle)    # Each Distances BTW Left Grid and Particle

            # Calc Particle Number in Each Cell
            ptlGridIndex = (particle.location/boundary.cell_length).astype(int)
            unique, counts = np.unique(ptlGridIndex, return_counts=True)
            cellCountInform = dict(zip(unique,counts))

            # Get count, sRDC
            pivot = 0
            for index, num in cellCountInform.items():
                count[index] = num
                sRDC[index] = np.sum(particle.location[pivot:pivot+num]) - index * num * boundary.cell_length
                relDist[pivot:pivot+num] = particle.location[pivot:pivot+num] - index * boundary.cell_length
                pivot += num

            # Produce GridCharge
            tmpChargeDensityGrid = np.zeros(shape=boundary.num_of_cells, dtype=np.float64)
            tmpChargeDensityGrid[1:] = sRDC[0:-1] + (boundary.cell_length * count[1:] - sRDC[1:])
            tmpChargeDensityGrid[0] = sRDC[-1] + (boundary.cell_length * count[0] - sRDC[0])

            # Devide by GridLength
            tmpChargeDensityGrid = (tmpChargeDensityGrid * particle.charge * particle.Sptcl_weight)/(boundary.cell_length * boundary.cell_volume) # ParticleCharge는 normalized 된 값 (1e -> 1)

            # Add Result to Final Result
            chargeDensityGrid += tmpChargeDensityGrid

        boundary.cell_charge_density = chargeDensityGrid * constants.elementary_charge

    @staticmethod
    def interpolateFieldsAtParticles(particles: list[Species], boundary: Boundary) -> None:
        
        for particle in particles:
            count = np.zeros(boundary.num_of_cells)
            relDist = np.zeros(particle.num_of_particle)

            ptlGridIndex = (particle.location/boundary.cell_length).astype(int)
            unique, counts = np.unique(ptlGridIndex, return_counts=True)
            tmpCellCountInform = dict(zip(unique, counts))

            pivot = 0
            for index, num in tmpCellCountInform.items():
                count[index] = num
                relDist[pivot:pivot+num] = particle.location[pivot:pivot+num] - index * boundary.cell_length
                pivot += num

            electricFieldL = np.zeros(shape=particle.num_of_particle, dtype=np.float64)
            electricFieldR = np.zeros(shape=particle.num_of_particle, dtype=np.float64)

            pivot = 0
            for index, num in tmpCellCountInform.items():
                electricFieldL[pivot:pivot+num] = boundary.cell_Efield[index]
                if index+1 == boundary.num_of_cells:
                    electricFieldR[pivot:pivot+num] = boundary.cell_Efield[0]
                else:
                    electricFieldR[pivot:pivot+num] = boundary.cell_Efield[index+1]
                pivot += num

            particle.Efield_at_particle = ((boundary.cell_length - relDist)*electricFieldL + relDist*electricFieldR) / boundary.cell_length

