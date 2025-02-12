import numpy as np
from psm.sim_mod import *
from psm.aic_mod import *
from psm.mle_mod import *
from psm.qua_mod import *
#________________________________________________
def choose_data(array, n_obs):
    """
    Selecciona de forma homogenea con tamaño n_obs
    Parameters:
        array (np.array): Array de trayectoria Logistic
        n_obs (float): Numero de observaciones 
            de la muestra final
    Returns:
        np.array con tamaño n_obs
    """
    l_a = len(array)
    size_jump = (l_a - 1) // (n_obs - 1)
    labels = [size_jump * (i - 1) for i in range(1, n_obs + 1)]
    array_1 = [array[label] for label in labels]
    return np.array(array_1)
#________________________________________________
def em_path(type_model,param,sigma,delta,n_iter_choose,start,end,n_obs,l_inf=9999999999.0):
    """
    Simulación de trayectoria segun el modelo
    de n_obs * n_iter_choose
    Parameters:
        type_model ('g'|'l'|'v'): Tipo de modelo
            g:Gompertz
            l:Logistic
            v:Von Bert
        param (float): Valor del parametro segun el modelo
            g:beta
            l:r
            v:kappa
        sigma (float): Valor sigma
        delta (float): Incremento del proceso de Wiener
        start (float): Valor inicial del movimiento Browniano
        n_iter_choose (int): Numero de iteraciones
            entre punto y punto
        n_obs (float): Numero de observaciones 
            de la muestra final
        l_inf (float): Limite superior
    Returns:
        np.array simulacion segun el modelo completa 
    """
    array_1 = sim(type_model,param,sigma,delta,start,n_iter_choose,l_inf)
    array_2 = sim(type_model,-param,sigma,delta,end,n_iter_choose,l_inf)
    #Ordenar en reversa
    array_2 = array_2[::-1]
    #Caso 1
    if np.max(array_1)<=np.max(array_2):
        for i in range(0, n_obs):
            if array_1[i] <= array_2[i]:
                arr_for = np.append(array_1[:i], array_2[i:])
                arr_for = np.append(np.array(start),arr_for)
                return arr_for
    #Caso 2
    elif np.max(array_1)>=np.max(array_2):
        for i in range(0, n_obs):
            if array_2[i] <= array_1[i]:
                arr_for = np.append(array_1[:i], array_2[i:])
                arr_for = np.append(np.array(start),arr_for)
                return arr_for
    #Caso 3
    else:
        return array_1
#________________________________________________
def em_array(type_model,array,param,sigma,delta,n_iter_choose,n_obs,l_inf):
    """
    Simulación de trayectoria segun el modelo para acompletar 
    de n_obs * n_iter_choose
    Parameters:
        type_model ('g'|'l'|'v'): Tipo de modelo
            g:Gompertz
            l:Logistic
            v:Von Bert
        array(np.array): array que se busca acompletar
        param (float): Valor del parametro segun el modelo
            g:beta
            l:r
            v:kappa
        sigma (float): Valor sigma
        delta (float): Incremento del proceso de Wiener
        start (float): Valor inicial del movimiento Browniano
        n_iter_choose (int): Numero de iteraciones para acompletar 
            entre punto y punto
        n_obs (float): Numero de observaciones 
            de la muestra final
        l_inf (float): Limite superior
    Returns:
        np.array simulacion segun el modelo completa 
        de tamaño n_obs * n_iter_choose
    """
    array_final = np.array([])
    for n in range(0,len(array)-1):
        start_for = array[n]
        end_for = array[n+1]
        #
        arr_for = em_path(type_model,param,sigma,delta,n_iter_choose,start_for,end_for,n_obs,l_inf=9999999999.0)
        if n == 0:
            array_final = arr_for
        else:
            array_final = np.vstack((array_final, arr_for))
    return array_final