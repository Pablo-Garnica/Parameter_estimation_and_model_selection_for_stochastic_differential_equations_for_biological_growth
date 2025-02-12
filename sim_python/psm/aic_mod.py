import numpy as np
from psm.integ_mod import *
#________________________________________________
def aic_g(array, delta, param, sigma):
    """
    Calcula el AIC para el modelo Gompertz.
    Parameters:
        array (numpy.ndarray): Array de simulacion Gompertz
        delta (float): Incremento del proceso de Wiener.
        param (float): Valor del parametro (beta)
        sigma (float): Valor sigma
    Returns:
        np.float valor aic
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
    """
    Calcula el AIC para el modelo Von Bert.
    Parameters:
        array (numpy.ndarray): Array de simulacion Von Bert
        delta (float): Incremento del proceso de Wiener.
        param (float): Valor del parametro (kappa)
        sigma (float): Valor sigma
        l_inf (float): Limite superior
    Returns:
        np.float valor aic
    """
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
    Calcula el AIC para el modelo Logistic.
    Parameters:
        array (numpy.ndarray): Array de simulacion Logistic
        delta (float): Incremento del proceso de Wiener.
        param (float): Valor del parametro (r)
        sigma (float): Valor sigma
    Returns:
        np.float valor aic
    """
    array_1 = -(param * array * (1.0 - array)) / (sigma**2 * (array**2))
    array_2 = (param**2 * ((array * (1.0 - array))**2)) / (sigma**2 * (array**2))
    i_1 = ito_integrate(array_1, array)
    i_2 = integrate(array_2, delta)
    aic = -2 * (i_1 - i_2 / 2)
    return aic
#________________________________________________
def aic_(type_dif, array, delta, param, sigma, l_inf=9999999999.0):
    """
    Calcula el AIC segun el modelo seleccionado
    Parameters:
        type_dif ('g'|'l'|'v'): Tipo de modelo
            g:Gompertz
            l:Logistic
            v:Von Bert
        array (numpy.ndarray): Array de simulacion Gompertz
        delta (float): Incremento del proceso de Wiener.
        param (float): Valor del parametro segun el modelo
            g:beta
            l:r
            v:kappa
        sigma (float): Valor sigma
        l_inf (float): Limite superior
    Returns:
        np.float valor aic
    """
    if type_dif in {'g'}:
        y = aic_g(array, delta, param, sigma)
    elif type_dif in {'l'}:
        y = aic_l(array, delta, param, sigma)
    elif type_dif in {'v'}:
        y = aic_v(array, delta, param, sigma, l_inf)
    else:
        y = None
    return y
