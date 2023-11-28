import numpy as np
from scipy import constants
from OneD_Periodic_CIC_ES_PIC.oneD_peridoic_CIC_PIC import PIC_1D1V_ES

simulation = {
    # Simulation Condition
    'number_of_timesteps'   : 1000,
    'timestep'              : 1e-6,

    # Plasma Condition
    'debye_length'          : 1e-3,

    # Result Save
    'every'                 : 1,
    'folder_name'            : 'dispersion_1',
    'round'                 : 5
}

boundary = {
    # Boundary Condition
    'num_of_cells'       : 1024,
    'cell_length'        : 1e-04,
    'cell_area'          : 1,
}

particles = (
    {
        'name'                  : 'ion',
        'mass'                  : constants.proton_mass/constants.electron_mass,
        'charge'                : 1,
        'temperature'           : 0,         # ev
        'Sptcl_number_density'  : 1e6,
        'Sptcl_weight'          : 100,
        'position_init'         : 'uniform', #'random' or 'uniform'
        'motion_fix'            : True,    # or True
        'boundary_condition'    : 'periodic'    # or 'reflect'
    },
    {
        'name'          : 'electron',
        'mass'          : 1,
        'charge'        : -1,
        'temperature'   : 1e-7,         # ev
        'Sptcl_number_density': 1e6,
        'Sptcl_weight'  : 100,
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
