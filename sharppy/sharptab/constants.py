''' Frequently used constants '''

consts = ["RMISSD", "ROCP", "ZEROCNK", "G", "TOL"]
funcs = ['MS2KTS', 'KTS2MS', 'M2FT', 'FT2M', 'QC']
__all__ = consts + funcs


# Constants
RMISSD = -9999.0            # Missing Flag
ROCP = 0.28571426           # R divided by Cp
ZEROCNK = 273.15            # Zero Celsius in Kelvins
G = 9.80665                 # Earth's Gravity
TOL = 1e-8                  # Floating Point Tolerance



# Functions
def MS2KTS(val):
    '''
    Convert meters per second to knots

    Inputs
    ------
        val     (float)         Speed (m/s)

    Returns
    -------
        Val converted to knots (float)
    '''
    return val * 1.94384449


def KTS2MS(val):
    '''
    Convert knots to meters per second

    Inputs
    ------
        val     (float)         Speed (kts)

    Returns
    -------
        Val converted to meters per second (float)
    '''
    return val / 1.94384449


def M2FT(val):
    '''
    Convert meters to feet

    Inputs
    ------
        val     (float)         Distance (m)

    Returns
    -------
        Val converted to feet (float)
    '''
    return val * 3.2808399


def FT2M(val):
    '''
    Convert feet to meters

    Inputs
    ------
        val     (float)         Distance (ft)

    Returns
    -------
        Val converted to meters (float)
    '''
    return val / 3.2808399


def QC(val):
    '''
    Check to see if given data is missing

    Inputs
    ------
        val     (float)         Value to check

    Returns
    -------
        1 if value is good
        0 if value is bad
    '''
    if val < -998.0 or val > 2.0e5: return 0
    return 1
