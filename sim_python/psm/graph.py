import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import json
import datetime
from psm.mle_mod import  mle_
from psm.qua_mod import  qua_var
#________________________________________________
def plt_sim(array,delta,title,y_title,plt_type =None, y_line_val=None,zoom=True,flg_export=False):
    df = pd.DataFrame(array)
    y = df[0]
    x = (df.index)*delta
    fig, ax = plt.subplots()
    # Crear la gráfica de dispersión
    ax.scatter(x,y, s=0.1)
    # Añadir etiquetas y título
    ax.set_title(title) 
    ax.set_xlabel('Time')
    if zoom:
        q1 = y.quantile(0.99)
        q3 = y.quantile(0.01)
        # Calcular el rango intercuartílico (IQR)
        iqr = q3 - q1
        # Definir los límites inferior y superior
        l_inf = q1 - 3 * iqr 
        l_sup = q3 + 3 * iqr
        # Aplicar zoom en el eje y
        ax.set_ylim(l_sup,l_inf)
    if y_line_val:
        ax.axhline(y_line_val, color='red',linewidth=0.8,label=f'Valor estimado:{y_line_val}')
        ax.legend()
    ax.set_ylabel(y_title)
    if flg_export:
        fig.savefig(f"{title} {plt_type}.png")
    return fig, ax
#________________________________________________
def get_d_sim(d_exe,d_f):
    d_plt = {
        'array' : pd.Series(d_f.get('array_0')),
        'delta' : d_exe.get('delta'),
        'title' : f"Simulation {d_exe['d_flow'].get('type_dif')}",
        'y_title' : '',
    }
    return d_plt
#________________________________________________
def get_d_sig(d_exe,d_f):
    ss = pd.Series(d_f.get('array_0'))
    type_dif = d_exe['type_dif']
    delta = d_exe.get('delta')
    l_inf = d_exe.get('l_inf')
    #
    qua_var_lamd = lambda x: qua_var(type_dif,np.array(x),delta,l_inf)
    ss = ss.expanding().apply(qua_var_lamd)
    #
    d_plt = {
        'array' : ss,
        'delta' : delta,
        'title' : f"Simulation {d_exe['d_flow'].get('type_dif')}",
        'y_title' : 'Sigma',
        'y_line_val' : d_exe.get('sigma'),
    }
    return d_plt
#________________________________________________
def get_d_par(d_exe,d_f):
    ss = pd.Series(d_f.get('array_0'))
    type_dif = d_exe['type_dif']
    delta = d_exe.get('delta')
    l_inf = d_exe.get('l_inf')
    #
    mle_lamd = lambda x: mle_(type_dif,np.array(x),delta,l_inf)
    ss = ss.expanding().apply(mle_lamd)
    #
    d_plt = {
        'array' : ss,
        'delta' : delta,
        'title' : f"Simulation {d_exe['d_flow'].get('type_dif')}",
        'y_title' : d_exe['d_flow'].get('model_name_param'),
        'y_line_val' : d_exe.get('param'),
    }
    return d_plt
#________________________________________________
def get_d_info(plt_type,d_exe,d_f):
    if plt_type == 'sim':
        return get_d_sim(d_exe,d_f)
    if plt_type == 'sig':
        return get_d_sig(d_exe,d_f)
    if plt_type == 'par':
        return get_d_par(d_exe,d_f)
    else:
        return {}
#________________________________________________
def json_extract(data,type_data,d_exe,d_f,flg_export=False):
    d_json = {
        'model':d_exe['d_flow'].get('type_dif'),
        'type_data':type_data,
        'param':d_exe.get('param'),
        'sigma':d_exe.get('sigma'),
        'delta':d_exe.get('delta'),
        'start':d_exe.get('start'),
        'n_iter':d_exe.get('n_iter'),
        'qua':d_f.get('qua'),
        'aic':d_f.get('aic'),
        'mle':d_f.get('mle'),
        'datetime': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'data': [None if np.isnan(x) else x for x in data]
    }
    name_file = f'{d_exe['d_flow'].get('type_dif')}_{type_data}.json'
    if flg_export:
        with open(name_file, "w", encoding="utf-8") as archivo:
            json.dump(d_json, archivo, indent=4, ensure_ascii=False,default=lambda x: None if x is np.nan else x)
        return d_json

#________________________________________________
def get_d_plt(type_dif,array,delta,type_graph,y_val=None,l_inf=9999999999.0):
    d_cat = {
        'g':'Gompertz',
        'l':'Logistic',
        'v':'Von Bert',
    }
    d_par = {
        'g':r"$\beta$",
        'l':r"$\rho$",
        'v':r"$\kappa$",
    }
    d_plt = {}
    ss = pd.Series(array)
    #
    if type_graph == 'mle':
        mle_lamd = lambda x: mle_(type_dif,np.array(x),delta,l_inf)
        ss = ss.expanding().apply(mle_lamd)
        d_plt = d_plt | {'y_title' : d_par[type_dif]}
        
    elif type_graph == 'qua':
        qua_var_lamd = lambda x: qua_var(type_dif,np.array(x),delta,l_inf)
        ss = ss.expanding().apply(qua_var_lamd)
        d_plt = d_plt | {'y_title' : r'$\sigma$'}
    else:
        d_plt = d_plt | {'y_title' : ''}
    d_plt = d_plt | {
        'array' : ss,
        'delta' : delta,
        'title' : f"Simulation {d_cat[type_dif]}",
        'y_line_val':y_val
    }
    return d_plt
#________________________________________________
def matrix_plt(axes_list, max_col=None):
    if not axes_list:
        raise ValueError("La lista de ejes no puede estar vacía")

    max_col = max_col or len(axes_list)  # Si max_col es None, usa el tamaño total de la lista
    l_row = int(np.ceil(len(axes_list) / max_col))  # Calcula el número de filas necesarias

    fig_final, axes = plt.subplots(l_row, max_col, figsize=(5 * max_col, 5 * l_row))
    axes = np.atleast_2d(axes)  # Asegurar que axes sea siempre una matriz 2D

    for ax_orig, ax_dest in zip(axes_list, axes.flat):  # Recorremos la matriz como array 1D
        # Copiar líneas
        for line in ax_orig.lines:
            ax_dest.plot(line.get_xdata(), line.get_ydata(), color=line.get_color(), linestyle=line.get_linestyle())

        # Copiar scatter
        for coll in ax_orig.collections:
            offsets = coll.get_offsets()
            sizes = coll.get_sizes()
            colors = coll.get_facecolor()
            ax_dest.scatter(offsets[:, 0], offsets[:, 1], s=sizes, color=colors)

        # Copiar el título del eje Y
        ax_dest.set_ylabel(ax_orig.get_ylabel())

        # Copiar la línea horizontal si existe
        for hline in ax_orig.get_lines():
            ax_dest.axhline(y=hline.get_ydata()[0], color=hline.get_color(), linestyle=hline.get_linestyle())

        # Copiar límites de los ejes
        ax_dest.set_xlim(ax_orig.get_xlim())
        ax_dest.set_ylim(ax_orig.get_ylim())

    # Ocultar los subplots vacíos si la cantidad de gráficos no llena toda la cuadrícula
    for ax in axes.flat[len(axes_list):]:
        ax.set_visible(False)

    return fig_final, axes
#________________________________________________
def interval_95(array):
    mean = np.mean(array)
    std_dev = np.std(array, ddof=1)  # Desviación estándar muestral
    std_err = std_dev / np.sqrt(len(array))  # Error estándar
    # Intervalo de confianza del 95%
    inf, sup = stats.norm.interval(0.95, loc=mean, scale=std_err)
    return inf, sup
#________________________________________________
def matrix_sel(df,type_exe,type_estimation):
    name_param = 'real_value' if type_estimation == 'mle' else 'sigma'
    ls_sel = [
        'type_model',
        name_param,
        f'{type_estimation}_{type_exe}',
        f'{type_estimation}_{type_exe}_inf',
        f'{type_estimation}_{type_exe}_sup'
    ]
    d_ren = {
        f'{type_estimation}_{type_exe}':'estimate',
        f'{type_estimation}_{type_exe}_inf':'interval_inf',
        f'{type_estimation}_{type_exe}_sup':'interval_sup',
    }
    return df[ls_sel].rename(columns = d_ren)