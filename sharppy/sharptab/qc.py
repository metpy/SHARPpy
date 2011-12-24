''' Quality Control Routines '''

__all__ = ['qc']



def qc(val):
    if val < -998.0 or val > 2.0e5: return 0
    return 1
