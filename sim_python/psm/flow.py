from psm.aic_mod import *
from psm.mle_mod import *
from psm.qua_mod import *
from psm.sim_mod import *
from psm.em_mod import *
from psm.graph import *
from numpy.random import default_rng
#________________________________________________
msg_seed = "Fije una semilla aleatoria,(Si no es un entero la semilla será 31416) "
msg_l_inf = "Introduce el valor del limite superior ,Si no es un float el limite será 9999999999"
#Mensajes
msg_type = '''Introduce el modelo que quieres calcular:
    "g : Gompertz"
    "l : Logistic"
    "v : Von Bert"
'''
msg_param  = '''-------------------------------------------------------------------
Seleccionaste el modelo {type_dif}
Para este modelo se requieren los siguientes argumentos
Introduce el valor {model_name_param} en ,0.6
'''
msg_sigma = 'Introduce el valor sigma en ,0.001: \n'
msg_delta = 'Introduce el valor delta en ,0.001: \n'
msg_start = 'Introduce el valor inicial de la simulación en ,0.01: \n'
msg_n_iter = 'Introduce el numero de valores para la simulacion en ,10,001: \n'
msg_type_execute = '''-------------------------------------------------------------------
    Que tipo de ejecición desea aplicar
    n : Ejecución normal del modelo
    t : Ejecución del modelo para seleccionar trayectoria
    e : Ejecución ocupando el algoritmo em
'''
msg_n_ex = '''type_dif:{type_dif}
qua:{qua}
aic:{aic}
mle:{mle}
---------------------------------------------
'''
msg_e_ex ='''type_dif:{type_dif}
    param mean {param_mean}
    sigma mean {sigma_mean}
---------------------------------------------
'''
msg_n_reply = 'Introduce el numero de trayectorias,100: \n'
msg_n_step = 'Introduce el tamaño de la muestra,10,000: \n'
msg_n_obs = 'Introduce el numero de observaciones,10: \n'
msg_plt_type = '''El tipo de grafico que requieres (si no se elije alguno de estos valores no se crearan loa graficos)
    sim = Simulación
    sig = Comportamiento sigma
    par = Comportamiento del parametro
'''
#________________________________________________
#Catalogos
d_type_m = {
    'g':'Gompertz',
    'l':'Logistic',
    'v':'Von Bert'
}
d_param_m = {
    'g':'beta',
    'l':'r',
    'v':'kappa'
}
#________________________________________________
def input_general_param(flg_print=False):
    """
    Flujo de parametros generales
    Returns:
        diccionario con parametros generales
    """
    d_flow = {}
    d_exe = {}
    try:
        seed = int(input(msg_seed))
    except:
        seed = 31416
    rng = default_rng(seed=seed)
    type_dif = input(msg_type)
    if type_dif=='v':
        try:
            l_inf = float(input(msg_l_inf))
        except:
            l_inf = 9999999999.0
        d_exe = d_exe | {'l_inf':l_inf}
    else:
        d_exe = d_exe | {'l_inf':9999999999.0}
    d_flow = d_flow | {'type_dif':d_type_m[type_dif], 'model_name_param':d_param_m[type_dif]}
    param = float(input(msg_param.format(**d_flow)))
    sigma = float(input(msg_sigma))
    delta = float(input(msg_delta))
    start = float(input(msg_start))
    n_iter = int(input(msg_n_iter))
    #
    d_exe = d_exe | {
        'type_dif':type_dif,
        'param':param,
        'sigma':sigma,
        'delta':delta,
        'start':start,
        'n_iter':n_iter,
        'd_flow':d_flow,
        'flg_print':flg_print,
    }
    #
    if flg_print:
        print(f'type_dif:{type_dif}')
        print(f'param:{param}')
        print(f'sigma:{sigma}')
        print(f'delta:{delta}')
        print(f'start:{start}')
        print(f'n_iter:{n_iter}')
    return d_exe
#________________________________________________
def normal_exe(type_dif,param,sigma,delta,start,n_iter,l_inf,d_flow,flg_print=False):
    """
    Flujo de ejecución normal
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
    """
    print(type_dif,param,sigma,delta,start,n_iter,l_inf,d_flow)
    array_0 = sim(type_dif,param,sigma,delta,start,n_iter,l_inf)
    qua =  qua_var(type_dif,array_0,delta)
    aic = aic_(type_dif, array_0, delta, param, sigma)
    mle = mle_(type_dif, array_0, delta)

    d_flow = d_flow | {
        'qua':qua,
        'aic':aic,
        'mle':mle,
    }
    d_return = {
        'array_0':array_0,
        'qua':qua,
        'aic':aic,
        'mle':mle,
    }
    x = input(msg_n_ex.format(**d_flow))
    if flg_print:
        print(msg_n_ex.format(**d_flow))
    return d_return
#________________________________________________
def trajectories_exe(type_dif,param,sigma,delta,start,n_iter,l_inf,d_flow,flg_print):
    """
    Flujo de ejecución para selección de trayectorias
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
    """
    n_reply = int(input(msg_n_reply))
    n_step = int(input(msg_n_step))

    array_0 = sim_reply(type_dif,param,sigma,delta,start,n_iter,n_reply,l_inf)
    array_1 = sim_choose(array_0,n_step,l_inf)
    qua =  qua_var(type_dif,array_1,delta)
    aic = aic_(type_dif, array_1, delta, param, sigma)
    mle = mle_(type_dif, array_1, delta)

    d_flow = d_flow | {
        'qua':qua,
        'aic':aic,
        'mle':mle,
    }
    #
    d_retun = {
        'array_sim_reply':array_0,
        'array_sim_choose':array_1,
        'qua':qua,
        'aic':aic,
        'mle':mle,
    }
    x = input(msg_n_ex.format(**d_flow))
    if flg_print:
        print(msg_n_ex.format(**d_flow))
    return d_retun
#________________________________________________
def em_exe(type_dif,param,sigma,delta,start,n_iter,l_inf,d_flow,flg_print):
    """
    Flujo de ejecución información incompleta
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
    """
    n_reply = int(input(msg_n_reply))
    n_obs = int(input(msg_n_obs))
    #
    if flg_print:
        print(f'n_reply:{n_reply}')
        print(f'n_obs:{n_obs}')
    #
    n_iter_chose = n_iter//n_obs
#________________________________________________
def input_type_exe(d_exe,flg_print=False):
    """
    Flujo de ejecución segun la elección de tipo
    Parameters:
        d_exe: output de input_general_param
    Returns:
        diccionario con información del la ejecución
    """
    type_execute = input(msg_type_execute)
    if flg_print:
        print(f'type_execute:{type_execute}')
    if type_execute == 'n':
        return normal_exe(**d_exe)|{'type_execute':type_execute}
    elif type_execute == 't':
        return trajectories_exe(**d_exe)|{'type_execute':type_execute}
    elif type_execute == 'e':
        return em_exe(**d_exe)|{'type_execute':type_execute}
    else:
        raise RuntimeError('type_execute debe estar en {"n","t","e"}')
#_________________________________________________________________
def input_plt(d_exe,d_f):
    """
    Flujo para exportar dato y graficos del modelo
    Parameters:
        d_exe: output de input_general_param
        d_f: output de normal_exe
    """
    if d_f['type_execute'] == 'n':
        plt_type = input(msg_plt_type)
        #_________________________________________________________________
        if plt_type in {'sim','sig','par'}:
            zoom = False if plt_type=='sim' else True
            flg_export = True
            d_plt = get_d_info(plt_type,d_exe,d_f)
            d_export = json_extract(d_plt['array'],plt_type,d_exe,d_f,flg_export)
            plt_sim(**d_plt,zoom=zoom,plt_type=plt_type,flg_export=flg_export)
            return d_export