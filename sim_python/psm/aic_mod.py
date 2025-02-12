import numpy as np
from psm.integ_mod import *
#________________________________________________
def aic_g(array, delta, param, sigma):
    """
    Calcula el AIC para el modelo Gompertz.
    """
    array_log = np.log(array)
    array_1 = -(param * array * array_log) / ((sigma**2) * (array**2))
    array_2 = ((param**2) * (array**2)) / ((sigma**2) * (array**2))
    
    i_1 = ito_integrate(array_1, array)
    i_2 = integrate(array_2, delta)
    
    aic = -2 * (i_1 - i_2 / 2)
    return aic
#________________________________________________
def aic_v(array, delta, param, sigma, l_inf=9999999999.0):
    array_1 = param / ((sigma**2) * (l_inf - array))
    array_2 = (param**2) / (sigma**2)
    #
    i_1 = ito_integrate(array_1, array)
    i_2 = delta * array_2
    #
    aic = -2 * (i_1 - i_2 / 2)
    return aic
#________________________________________________
def aic_l(array,delta, param, sigma):
    """
    Calcula el valor del AIC (Akaike Information Criterion) para el modelo Logístico.

    Parameters
    ----------
    delta : float
        Incremento del proceso de Wiener
    r : float
        Valor de r
    sigma : float
        Valor de sigma
    npoints : int
        Tamaño de la simulación
    path : numpy.ndarray
        Array de valores donde se guarda la simulación

    Returns
    -------
    aic : float
        Akaike Information Criterion
    """
    array_1 = -(param * array * (1.0 - array)) / (sigma**2 * (array**2))
    array_2 = (param**2 * ((array * (1.0 - array))**2)) / (sigma**2 * (array**2))
    i_1 = ito_integrate(array_1, array)
    i_2 = integrate(array_2, delta)
    aic = -2 * (i_1 - i_2 / 2)
    return aic
#________________________________________________
def aic_(type_dif, array, delta, param, sigma, l_inf=9999999999.0):
    if type_dif in {'g'}:
        y = aic_g(array, delta, param, sigma)
    elif type_dif in {'l'}:
        y = aic_l(array, delta, param, sigma)
    elif type_dif in {'v'}:
        y = aic_v(array, delta, param, sigma, l_inf)
    else:
        y = None
    return y
