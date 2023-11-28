import numpy as np
from scipy import constants
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import factorized

from .maxwell_solver.oneD_periodic_solver import OneD_Periodic_Solver

"""
Boundary.py
=========================================
Geometry                : 1d Cartesian
EM_Boundary_Condition   : Periodic
=========================================
"""
class Boundary:

    def __init__(self, boundary_info : dict):
        # Init
        num_of_cells = boundary_info['num_of_cells']
        cell_length = boundary_info['cell_length']
        cell_area = boundary_info['cell_area']

        # System Parameters
        self.xMin: np.float64 = 0
        self.xMax: np.float64 = num_of_cells * cell_length
        self.cell_length: np.float64 = cell_length
        self.num_of_cells = num_of_cells       
        self.cell_volume = cell_length * cell_area
        self.system_volume = (self.xMax - self.xMin) * cell_area                                

        # Cell Quantity
        self.cell_charge_density : np.ndarray = np.zeros(shape=num_of_cells, dtype=np.float64)      # SI, C/(m^3)
        self.cell_Efield : np.ndarray = np.zeros(shape=num_of_cells, dtype=np.float64)              # SI, N/C
        self.cell_potential : np.ndarray = np.zeros(shape=num_of_cells, dtype=np.float64)           # SI, N*m/C

        # Solvers
        self.maxwell_solver = OneD_Periodic_Solver(
            cell_length = self.cell_length, 
            num_of_cells = self.num_of_cells
        )
        
    """
    Use external solver "OneD_Periodic_Solver"
    """
    def calc_field(self):
        self.maxwell_solver.poissonSolver(
            # input
            i_charge_density = self.cell_charge_density,
            # output
            o_Efield = self.cell_Efield,
            o_potential = self.cell_potential
        )

    def get_totalEME(self) -> np.float64:
        # 0.5 * ep0 * (E^2)
        return 0.5 * constants.epsilon_0 * np.sum(self.cell_Efield*self.cell_Efield) * self.cell_volume
    
    def get_totalEME2(self) -> np.float64:
        return 0.5 * np.sum(self.cell_charge_density * self.cell_potential) * self.cell_volume

    """
    Modify 'location_array' with Periodic boundary condition.
    Function directly get location_array's address and modify it. (numpy array)
    ex)
    if x ~ (0, 10)
    particle location -1 goes to 9
    particle location 21 goes to 1
    """
    def periodic_boundary_check(self, location_array: np.ndarray) :
        # left to right
        while True:
            leftFinish = None
            for index, r in enumerate(location_array):
                if r >= self.xMin:
                    leftFinish = index
                    break
            if leftFinish == 0:
                break
            location_array[0:leftFinish] += self.xMax

        # right to left
        while True:
            rightStart = None
            for index, r in reversed(list(enumerate(location_array))):
                if r < self.xMax:
                    rightStart = index+1
                    break
            if rightStart == len(location_array):
                break
            location_array[rightStart:] -= self.xMax

    """
    Modify 'location_array' with Reflective boundary condition.
    Function directly get location_array's address and modify it. (numpy array!)
    ex)
    if x ~ (0, 10)
    particle location -1 goes to 1
    particle location 21 goes to ??? (Todo)
    """
    def reflective_boundary_check(self, location_array: np.ndarray) :
        # Todo!!
        None
