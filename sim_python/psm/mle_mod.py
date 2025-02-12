import numpy as np
from psm.integ_mod import *
#________________________________________________
def mle_g(array, delta):
    """
    Calcula el MLE (Maximum Likelihood Estimator) para el modelo Gompertz.

    Parameters:
        npoints (int): Número de observaciones.
        path (numpy.ndarray): Observaciones de la SDE (stochastic differential equation).
        delta (float): Incremento del proceso de Wiener.

    Returns:
        float: Estimado b.
    """
    l_array = np.log(array)
    l_a = len(array)
    num = (
        (l_a - 1) * 
        np.sum(l_array[:-1] * l_array[1:]) - 
        np.sum(l_array[1:]) * np.sum(l_array[:-1])
    )
    dem = (
        (l_a - 1) * 
        np.sum(l_array[:-1]**2) - 
        np.sum(l_array[:-1])**2
    )
    j = num / dem
    bhat = -np.log(j) / delta
    return bhat
#________________________________________________
def mle_l(array, delta):
    """
    Calcula el MLE (Maximum Likelihood Estimator) para el modelo Logístico.

    Parameters:
        npoints (int): Número de observaciones.
        path (numpy.ndarray): Observaciones de la SDE (stochastic differential equation).
        delta (float): Incremento del proceso de Wiener.

    Returns:
        float: Estimador r.
    """
    array_0 = (1.0 - array)
    array_1 = array_0/ array
    i_1 = ito_integrate(array_1, array)
    i_2 = integrate(array_0**2, delta)
    rhat = i_1 / i_2
    return rhat
#________________________________________________
def mle_v(array, delta, l_inf=9999999999.0):
    """
    Calcula el MLE (Maximum Likelihood Estimator) para el modelo Von Bertalanffy.

    Parameters:
        npoints (int): Número de observaciones.
        path (numpy.ndarray): Observaciones de la SDE (stochastic differential equation).
        delta (float): Incremento del proceso de Wiener.
        l_inf (float): Límite superior.

    Returns:
        float: Estimador kappa.
    """
    l_a = len(array)
    time = delta * (l_a - 1)
    array_int = 1.0 / (l_inf - array)
    i_1 = ito_integrate(array_int, array)
    kappahat = i_1 / time
    return kappahat
#________________________________________________
def mle_(type_dif, array, delta, l_inf=9999999999.0):
    if type_dif in {'g'}:
        y = mle_g(array, delta)
    elif type_dif in {'l'}:
        y = mle_l(array, delta)
    elif type_dif in {'v'}:
        y = mle_v(array, delta, l_inf)
    else:
        y = None
    return y
#________________________________________________
def mle_times(array, delta,size ,l_inf=9999999999.0):
    """
    Calcula el MLE (Maximum Likelihood Estimator).

    Parameters:
        path (numpy.ndarray): Observaciones.
        delta (float): Incremento del proceso de Wiener.
        linf (float): Límite superior.
        size (int): Número de segmentos para calcular MLE.
        npoints (int): Número total de observaciones.

    Returns:
        tuple: Arrays con los valores de sigmahat y kappahat.
    """
    sigmahat = np.array([])
    kappahat = np.array([])
    l_a = len(array)
    geo_brown = array - l_inf
    size_jump = (l_a - 1) // (size - 1)
    array_1 = np.array([size_jump * (i + 1) for i in range(size)])
    #
    for i, n in enumerate(array_1):
        sigmahat[i], kappahat[i] = mle_gb(geo_brown[:n*size_jump], delta)
    return sigmahat, kappahat
#________________________________________________
def mle_gb(array, delta):
    """
    Calcula el MLE (Maximum Likelihood Estimator) para el modelo Geo Brown.

    Parameters:
        path (numpy.ndarray): Observaciones transformadas.
        delta (float): Incremento del proceso de Wiener.
        npoints (int): Número de observaciones.

    Returns:
        tuple: sigmahat y kappahat.
    """
    n_points = len(array)
    # dif_array = np.diff(array)
    log_array_2 = np.log(array[1:])
    log_array_1 = np.log(array[:-1])
    
    a = np.sum(log_array_2)
    b = np.sum(log_array_1)
    c = np.sum(log_array_2**2)
    d = np.sum(log_array_2 * log_array_1)
    e = np.sum(log_array_1**2)
    
    sigmahat = np.sqrt(
        (-a**2 + 2*a*b - b**2 + c*n_points - 2*d*n_points + e*n_points) /
        ((n_points**2) * delta)
    )
    kappahat = -((a - b) / (delta * (n_points - 1)) + (sigmahat**2) / 2)
    return sigmahat, kappahat
