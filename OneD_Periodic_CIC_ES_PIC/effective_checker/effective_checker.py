import scipy.constants as constant

class Effective_checker:

    @staticmethod
    def isEffective(dx, dt, debye_length) -> tuple[bool, dict]:
        debye_check, dx_debye_ratio = Effective_checker.__debyeCheck(dx, debye_length)
        cfl_check, mindt = Effective_checker.__cflCheck(dx, dt)

        return (debye_check and cfl_check), {'debye_ratio': dx_debye_ratio, 'mindt': mindt}
    
    """
    Check
    dx << debyeLength
    """
    @staticmethod
    def __debyeCheck(dx, debye_length) -> tuple[bool, float]:
        if(dx < debye_length*0.5):
            debyeCheck = True
        else:
            debyeCheck = False
        
        dxRatio = debye_length/dx

        return (debyeCheck, dxRatio)

    """
    Check
    c * dt < t_cfl
    if 1d -> t_cfl = dx
    """
    @staticmethod
    def __cflCheck(dx, dt) -> tuple[bool, float]:
        if dt < (dx/constant.speed_of_light):
            cflCheck = True
        else:
            cflCheck = False

        mindt = (dx/constant.speed_of_light)

        return (cflCheck, mindt)