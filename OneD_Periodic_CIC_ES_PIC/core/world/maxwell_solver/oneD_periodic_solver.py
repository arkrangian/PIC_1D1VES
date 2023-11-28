import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import factorized
from scipy import constants

"""
OneD_Peridoic_solver.py
=========================================
Geometry                : 1d Cartesian
EM_Boundary_Condition   : Periodic

Solve only V and E (not B)
=========================================
"""
class OneD_Periodic_Solver:

    def __init__(self, cell_length, num_of_cells):
        self.cell_length = cell_length
        self.num_of_cell = num_of_cells
        
        self.periodic_potential_solver = self.__potentialSolverMaker(num_of_cells)

    def __potentialSolverMaker(self, num_of_cell):
        n = num_of_cell - 1
        data = [[2]*n, [-1]*(n-1), [-1]*(n-1)]
        offsets = [0,1,-1]

        raplacianDiag = diags(data, offsets, shape=(n, n), format='csc', dtype=np.float64)
        return factorized(raplacianDiag)
    
    """
    input) 
        grid_charge_density array ( [rho_0, rho_1, ... , rho_GN-1] )
    output) 
        grid_potential array ( [V_0, V_1, ... , V_GN-1] )
        grid_Efield array ( [E_0, E_1, ... , E_GN-1] )
    """
    def poissonSolver(self, i_charge_density: np.ndarray, o_Efield: np.ndarray, o_potential: np.ndarray) :
        # Calculate Potential
        o_potential[1:] = self.periodic_potential_solver(i_charge_density[1:]) * (self.cell_length**2)
        o_potential[0] = 0

        # Calculate Efield
        o_Efield[1:-1] = -(o_potential[2:] - o_potential[:-2]) / (2 * self.cell_length)
        o_Efield[0] = -(o_potential[1] - o_potential[-1]) / (2 * self.cell_length)
        o_Efield[-1] = -(o_potential[0] - o_potential[-2]) / (2 * self.cell_length)

        # Devide by epsion-zero
        o_potential /=  constants.epsilon_0
        o_Efield /= constants.epsilon_0