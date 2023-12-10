import numpy as np
from scipy import constants
from OneD_Periodic_CIC_ES_PIC.oneD_peridoic_CIC_PIC import PIC_1D1V_ES

simulation = {
    # Simulation Condition
    'number_of_timesteps'   : 10000,
    'timestep'              : 1e-8,

    # Plasma Condition
    'debye_length'          : 0.0235,

    # Result Save
    'every'                 : 50,
    'folder_name'            : 'SLPIC_base2',
    'round'                 : 5
}

boundary = {
    # Boundary Condition
    'num_of_cells'       : 128,
    'cell_length'        : 1e-2,
    'cell_area'          : 5e-3,
}

particles = (
    {
        'name'                  : 'ion',
        'mass'                  : constants.proton_mass/constants.electron_mass,
        'charge'                : 1,
        'temperature'           : 0.1,         # ev
        'Sptcl_number_density'  : 1e7,
        'Sptcl_weight'          : 1000,
        'position_init'         : 'random', #'random' or 'uniform'
        'motion_fix'            : False,    # True or False
        'boundary_condition'    : 'periodic'    # or 'reflect'(not implemented todo!)
    },
    {
        'name'          : 'electron',
        'mass'          : 1,
        'charge'        : -1,
        'temperature'   : 1,         # ev
        'Sptcl_number_density': 1e7,
        'Sptcl_weight'  : 1000,
        'position_init' : 'random', # or 'uniform'
        'motion_fix'    : False,    # or True
        'boundary_condition'    : 'periodic'
    }
)

simul = PIC_1D1V_ES(
    info_simulation = simulation,
    info_boundary = boundary,
    info_particles = particles
    )

simul.simulStart()
