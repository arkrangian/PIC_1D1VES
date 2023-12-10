import numpy as np
from tqdm import tqdm
import csv
import os
import json

from .core.world.boundary import Boundary
from .core.particle.particle_set import Particle_Set
from .core.particle.species import Species
from .core.interpolation_lib.oneD_CIC_interpolator import Interpolater_CIC
from .effective_checker.effective_checker import Effective_checker

class PIC_1D1V_ES:

    def __init__(self, 
                 info_simulation: dict,
                 info_boundary: dict, 
                 info_particles: tuple[dict]
                 ) -> None:

        # Save simulation info
        self.profile_info = [info_simulation, info_boundary, info_particles]
        
        # Simulation Condition
        self.timestep: np.float64 = info_simulation['timestep']                               # Timestep
        self.number_of_timesteps: int = info_simulation['number_of_timesteps']                # Total iterations of simulation
        self.simulTime: np.float64 = self.timestep * self.number_of_timesteps                 # Total simulation time
        self.debye_length: np.float64 = info_simulation['debye_length']                       # DebyeLength of plasma

        # For Result save
        self.every = info_simulation['every']
        self.folder_name = info_simulation['folder_name']
        self.round = info_simulation['round']

        # Create simulation boundary
        self.boundary: Boundary = Boundary(boundary_info=info_boundary)

        # Create particles
        self.particles: Particle_Set = Particle_Set(infos=info_particles)

    def make_fileIO(self):
        os.makedirs(self.folder_name, exist_ok=True)

        path_profile = os.path.join(self.folder_name, 'profile.json')
        with open(path_profile, 'w') as json_file:
            json.dump(self.profile_info, json_file, indent=2)
        
        path_TotalE = os.path.join(self.folder_name, 'total_E.csv')
        self.f_TotalE = open(path_TotalE,'w',newline='')
        self.wr_Total_E = csv.writer(self.f_TotalE)

        path_Total_KE = os.path.join(self.folder_name, 'KE.csv')
        self.f_KE = open(path_Total_KE,'w',newline='')
        self.wr_KE = csv.writer(self.f_KE)

        path_Total_EPE = os.path.join(self.folder_name, 'EPE.csv')
        self.f_EPE = open(path_Total_EPE,'w',newline='')
        self.wr_EPE = csv.writer(self.f_EPE)

        path_cellE = os.path.join(self.folder_name, 'cell_E.csv')
        self.f_cellE = open(path_cellE,'w',newline='')
        self.wr_cellE = csv.writer(self.f_cellE)

        self.f_ion0VelDist, self.wr_ion0VelDist = self.__createFileIO(folder_name=self.folder_name, file_name= 'ion0_velDist.csv')
        self.f_electronVelDist, self.wr_electronVelDist = self.__createFileIO(folder_name=self.folder_name, file_name= 'electron_velDist.csv')

    def __createFileIO(self, folder_name: str, file_name: str):
        path = os.path.join(folder_name, file_name)
        fileIO = open(path, 'w', newline='')
        csvIO = csv.writer(fileIO)
        return fileIO, csvIO

    def effectiveCheck(self):
        isEffective, effective_inform = Effective_checker.isEffective(
            dx = self.boundary.cell_length, 
            dt = self.timestep, 
            debye_length= self.debye_length
        )
        
        print(effective_inform)

        if(isEffective == False):
            neglect = input("keep going? (type 'y' else terminate) : ")
            if(neglect == 'y'):
                None
            else:
                exit()

    def initSet(self):
        # Generate Particles
        self.particles.generate_particles(boundary = self.boundary)    

        # Initial Charge Density Interpolate (Particle -> Grid)
        Interpolater_CIC.chargeDensityInterpolate(
            particles=self.particles.particles,
            boundary=self.boundary
            )
        
        # Initial Calc of V,E
        self.boundary.calc_field()

        # InterPolate grid(E) to Particles
        Interpolater_CIC.interpolateFieldsAtParticles(
            particles=self.particles.particles,
            boundary=self.boundary
        )

    def simulStart(self) -> None:
        # Effective Check
        self.effectiveCheck()
        # Make FileIO
        self.make_fileIO()
        # InitSet
        self.initSet()

        for i in tqdm(range(self.number_of_timesteps), desc="simulation..", mininterval=1):

            # 1. Move Particles
            self.particles.push_particles(
                timestep=self.timestep,
                func_boundary_check=self.boundary.periodic_boundary_check
            )

            # 2. Interpolate particles(charge density) to grids
            Interpolater_CIC.chargeDensityInterpolate(
                particles=self.particles.particles, 
                boundary=self.boundary
                )

            # 3. Calculate V,E on the boundary
            self.boundary.calc_field()

            # 4. Interpolate grids(E) to particles
            Interpolater_CIC.interpolateFieldsAtParticles(
                particles=self.particles.particles, 
                boundary=self.boundary)

            # 5. Diagnotics & Save
            self.results_saver(i)
        
        self.simulationEnd()

    def results_saver(self, i):
        if(i%self.every == 0):
            # Kinetic, Potential, Total Energy
            KE = self.particles.get_totalKE()
            PE = self.boundary.get_totalEME()
            TotalE = KE + PE
            self.wr_Total_E.writerow([i, i*self.timestep, TotalE])
            self.wr_KE.writerow([i, i*self.timestep, KE])
            self.wr_EPE.writerow([i, i*self.timestep, PE])

            # Cell Electric Field
            self.wr_cellE.writerow(list(self.boundary.cell_Efield))

            # Ion(Maxwellian) distribution
            ion0_vel_width = 1000
            ion0_hist, ion0_edge = self.particles.get_particle_velDist('ion', vel_width=ion0_vel_width)
            self.wr_ion0VelDist.writerow([i, i*self.timestep] + list(ion0_edge))
            self.wr_ion0VelDist.writerow([i, i*self.timestep] + list(ion0_hist))

            elec_vel_width = 50000
            elec_hist, elec_edge = self.particles.get_particle_velDist('electron', vel_width=elec_vel_width)
            self.wr_electronVelDist.writerow([i, i*self.timestep] + list(elec_edge))
            self.wr_electronVelDist.writerow([i, i*self.timestep] + list(elec_hist))
        
    
    def simulationEnd(self) -> None:
        self.f_TotalE.close()
        self.f_KE.close()
        self.f_EPE.close()
        self.f_cellE.close()
        self.f_ion0VelDist.close()
        self.f_electronVelDist.close()