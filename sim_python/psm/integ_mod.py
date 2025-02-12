import numpy as np
#________________________________________________
def integrate(array, delta):
    """
    Calcula la integral de Riemann utilizando la regla del trapecio.

    Parameters
    ----------
    npoints : int
        Tamaño del arreglo de valores sobre el cual se calculará la integral
    path : numpy.ndarray
        Array de valores (imagen de la función) sobre los cuales se calculará la integral
    delta : float
        Distancia entre las particiones (se asumió partición uniforme)

    Returns
    -------
    integral : float
        Resultado de la integral de Riemann
    """
    # Calcula los términos promedios entre puntos consecutivos
    array_avg = (array[:-1] + array[1:]) / 2.0  # Promedio de cada subintervalo
    # Calcula la integral como suma de las áreas de los trapecios
    integ = delta * np.sum(array_avg)
    return integ
#________________________________________________
def ito_integrate(array_x,array_w):
    """
    Calcula la integral de Itô utilizando la fórmula de diferencias finitas.

    Parameters
    ----------
    npoints : int
        Tamaño del arreglo de valores sobre el cual se calculará la integral
    path : numpy.ndarray
        Array de valores (dominio de la función) sobre los cuales se calculará la integral
    W : numpy.ndarray
        Array de valores (imagen de la función) sobre los cuales se calculará la integral

    Returns
    -------
    integral : float
        Resultado de la integral de Itô
    """
    # Diferencias sucesivas de W
    delta_w = np.diff(array_w)
    # Promedio entre puntos consecutivos de path
    averages_path = (array_x[:-1] + array_x[1:]) / 2.0
    # Producto término a término y suma para calcular la integral
    integral = np.sum(delta_w * averages_path)
    return integral
