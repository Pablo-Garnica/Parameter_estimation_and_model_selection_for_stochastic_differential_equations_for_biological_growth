import numpy as np
from psm.integ_mod import integrate
#________________________________________________
def qua_var_g(array, delta):
    """
    Calcula la variación cuadrática para el modelo Gompertz

    Parameters
    ----------
    npoints : int
        Número de observaciones
    path : numpy.ndarray
        Observaciones de la SDE (ecuación diferencial estocástica)
    delta : float
        Incremento del proceso de Wiener

    Returns
    -------
    sigmahat : float
        Estimador sigma
    """
    # Inicialización de las variables
    l_a = len(array)
    array = np.log(array)  # Logaritmo natural de los valores en path
    num = (
        (l_a - 1) * 
        np.sum(array[:-1] * array[1:]) - 
        np.sum(array[1:]) * np.sum(array[:-1])
    )
    dem = (
        (l_a - 1) * 
        np.sum(array[:-1] ** 2) - 
        np.sum(array[:-1]) ** 2
    )
    j = num / dem
    k = (
        (
            j * np.sum(array[:-1]) - 
            np.sum(array[1:])
        ) / 
        (l_a - 1)
    )
    h = (
        np.sum(
            (array[1:] - j * array[:-1] + k) ** 2
        ) / 
        (l_a - 1)
    )
    bhat = -np.log(j) / delta
    sigmahat = (
        np.sqrt(
            (h * 2.0 * bhat) / 
            (1.0 - np.exp(-2.0 * bhat * delta))
        )
    )
    return sigmahat
#________________________________________________
def qua_var_l(array, delta):
    """
    Calcula la variación cuadrática para el modelo Lofistic

    Parameters
    ----------
    npoints : int
        Número de observaciones
    path : numpy.ndarray
        Observaciones de la SDE (ecuación diferencial estocástica)
    delta : float
        Incremento del proceso de Wiener

    Returns
    -------
    sigmahat : float
        Estimador sigma
    """
    # Calcular num y dem mediante un bucle
    differences = np.diff(array)
    num = np.sum(differences**2)  # Numerador: suma de las diferencias al cuadrado
    # Calcula la suma de los cuadrados de los valores consecutivos
    dem = np.sum(array[1:]**2 + array[:-1]**2)  # Denominador
    # Calcula el estimador sigma
    sigmahat = np.sqrt(2.0 * num / (dem * delta))
    return sigmahat
#________________________________________________
def qua_var_v(array, delta, l_inf=9999999999.0):
    """
    Calcula la variación cuadrática para el modelo Lofistic con un límite superior

    Parameters
    ----------
    npoints : int
        Número de observaciones
    path : numpy.ndarray
        Observaciones de la SDE (ecuación diferencial estocástica)
    delta : float
        Incremento del proceso de Wiener
    l_inf : float
        Límite superior

    Returns
    -------
    sigmahat : float
        Estimador sigma
    """
    array_inf = (l_inf - array) ** 2  # Se calcula el cuadrado de la diferencia con el límite superior
    # Llamada a la función Integrate para calcular 'dem'
    dem = integrate(array_inf, delta)
    # Calcular num sumando las diferencias al cuadrado entre puntos consecutivos de path
    num = np.sum((array[1:]-array[:-1])**2)
    # Calcular sigmahat
    sigmahat = np.sqrt(num / dem)
    return sigmahat
#________________________________________________
def qua_var(type_dif, array, delta, l_inf=9999999999.0):
    if type_dif in {'g'}:
        y = qua_var_g(array, delta)
    elif type_dif in {'l'}:
        y = qua_var_l(array, delta)
    elif type_dif in {'v'}:
        y = qua_var_v(array, delta, l_inf)
    else:
        y = None
    return y