import numpy as np
from scipy.stats import truncnorm
#________________________________________________
def difussion_param(x,sigma,type_dif,l_inf=9999999999.0):
    """
    Calcula el parametro diffusion
    Parameters:
        x (float): Valor en el que se evalua la diffusion
        sigma (float): Valor sigma
        type_dif ('g'|'l'|'v'): Tipo de modelo
            g:Gompertz
            l:Logistic
            v:Von Bert
        l_inf (float): Limite superior
    Returns:
        np.float parametro diffusion
    """
    if type_dif in {'l','g'}:
        y=sigma*x
    elif type_dif in {'v'}:
        y=sigma*(l_inf-x)
    else:
        print('error')
    return y
#________________________________________________
def drift_param(x,param,type_dif,l_inf=9999999999.0):
    """
    Calcula el parametro drift
    Parameters:
        x (float): Valor en el que se evalua la diffusion
        sigma (float): Valor sigma
        type_dif ('g'|'l'|'v'): Tipo de modelo
            g:Gompertz
            l:Logistic
            v:Von Bert
        l_inf (float): Limite superior
    Returns:
        np.float parametro drift
    """
    if type_dif in {'g'}:
        y=-param*x*np.log(x)
    elif type_dif in {'l'}:
        y=param*x*(1-x)
    elif type_dif in {'v'}:
        y=param*(l_inf-x)
    else:
        y = None
    return y
#________________________________________________
def milstein_step(start,delta,drift,diff,sigma):
    """
    Aplica el metodo Milstein para resolver ecuaciones estocasticas
    Parameters:
        start (float): Valor inicial del movimiento Browniano
        delta (float): Incremento del proceso de Wiener
        drift (float): Valor drift
        diff (float): Valor driffucion
        sigma (float): Valor sigma
    Returns:
        np.float Simulación de resolucion de ecuaciones estocasticas
    """
    x = np.random.normal()
    #Paso browneano
    w = np.sqrt(delta) * x
    end = start + drift * delta + diff * w + (0.5) * diff * sigma * (w**2 - delta)
    return end
#________________________________________________
def sim(type_dif,param,sigma,delta,start,n_iter,l_inf=9999999999.0,flg_print=False):
    """
    Simulación de trayectoria segun el modelo
    Parameters:
        type_dif ('g'|'l'|'v'): Tipo de modelo
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
        n_iter (int): Numero de iteraciones
        l_inf (float): Limite superior
        flg_print (bool): Flg para ejecutar prints o no
    Returns:
        np.array simulacion segun el modelo
    """
    array = np.array([])
    for i in range(0,n_iter):
        if flg_print:
            iter_batch = np.ceil(n_iter/100)
            if i%iter_batch==0:
                txt = f'%{round((i/n_iter) * 100, 2)}'
                print(txt)
        if i==0:
            x = start
            array = np.append(array,start)
        else:
            drift = drift_param(x,param,type_dif,l_inf)
            diff = difussion_param(x,sigma,type_dif,l_inf)
            x = milstein_step(x,delta,drift,diff,sigma)
            array = np.append(array,x)
    if flg_print:
        print(f'%100.0')
    return array
#________________________________________________
def sim_reply(type_dif,param,sigma,delta,start,n_iter,n_reply,l_inf=9999999999.0,flg_print=False):
    """
    Simulación de trayectoria segun el modelo
    Parameters:
        type_dif ('g'|'l'|'v'): Tipo de modelo
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
        n_iter (int): Numero de iteraciones
        n_reply (int): Numero de trayectorias generadas
        l_inf (float): Limite superior
        flg_print (bool): Flg para ejecutar prints o no
    Returns:
        np.array simulaciones de tamaño n_iter * n_reply
    """
    for i in range(0,n_reply):
        if flg_print:
            iter_batch = np.ceil(n_reply/100)
            if i%iter_batch==0:
                txt = f'%{round((i/n_reply) * 100, 2)}'
                print(txt)
        if i == 0:
            array = sim(type_dif,param,sigma,delta,start,n_iter,l_inf)
        else:
            array = np.vstack((
                array, 
                sim(type_dif,param,sigma,delta,start,n_iter,l_inf)
            ))
    if flg_print:
        print(f'%100.0')
    return array
#________________________________________________
def sim_choose(array, n_steps, l_inf=9999999999.0):
    """
    Elige una trayectoria a partir de una serie de trayectorias.
    Parameters:
        array (numpy.ndarray): Matriz de iteraciones con los valores simulados (shape: niter x npoints).
        n_steps (int): Número de pasos para el algoritmo.
        l_inf (float, optional): Límite superior. Si no se especifica, se establece a un valor muy alto.
    Returns:
        numpy.ndarray: Matriz de trayectorias seleccionadas (shape: n_steps).
    """
    array_out = np.array([])
    try:
        l_a = len(array[0])
        array_out = np.append(array_out,array[0, 0])
    except:
        l_a = 1
        array_out = np.append(array_out,array[0])
    size_jump = (l_a - 1) // (n_steps - 1)
    ls_points = [size_jump * (i - 1) + 1 for i in range(1, n_steps + 1)]
    #
    for i in range(1, n_steps):
        array_out = np.append(array_out, choose_obs(array[:, ls_points[i]], array_out[i - 1], l_inf))
    return array_out
#________________________________________________
def choose_obs(array, obs, l_inf=9999999999.0):
    """
    Realiza una selección basada en observaciones previas y distribuciones truncadas.
    Parameters:
        array (numpy.ndarray): Muestra actual (1D array).
        obs (float): Observación previa.
        l_inf (float): Límite superior.
    Returns:
        float: Nueva observación seleccionada.
    """
    variance = np.var(array)
    stddev = np.sqrt(variance)
    sorted_array = np.sort(array)
    if variance < 0.0000001:
        variance = 0.0001
    # Generar una muestra de una distribución normal truncada
    lower_bound = obs - variance
    upper_bound = l_inf
    trunc_normal_sample = truncnorm.rvs(
        (lower_bound) / stddev,
        (upper_bound - obs) / stddev,
        loc=obs,
        scale=stddev,
        size=1
    )[0]
    # Seleccionar el valor más cercano que sea mayor o igual al trunc_normal_sample
    for value in sorted_array:
        if value >= trunc_normal_sample:
            return value
    # Si no se encuentra un valor mayor, devolver el último
    return sorted_array[-1]