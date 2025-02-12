import numpy as np
#________________________________________________
def integrate(array, delta):
    """
    Calcula la integral de Riemann utilizando la regla del trapecio.
    Parameters:
        array (np.array): Array de trayectoria
        delta (float): Incremento del proceso de Wiener
    Returns:
        float: Resultado de la integral de Riemann
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
        array_x (np.array): Array de trayectoria
        array_w (np.array): Array de trayectoria aplicando la funcion 
            definida por el modelo
    Returns
    float Resultado de la integral de Itô
    """
    # Diferencias sucesivas de W
    delta_w = np.diff(array_w)
    # Promedio entre puntos consecutivos de path
    averages_path = (array_x[:-1] + array_x[1:]) / 2.0
    # Producto término a término y suma para calcular la integral
    integral = np.sum(delta_w * averages_path)
    return integral