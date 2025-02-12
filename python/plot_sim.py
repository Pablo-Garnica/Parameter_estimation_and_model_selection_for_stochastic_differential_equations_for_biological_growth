#_______________________________________________________________________________
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
#_______________________________________________________________________________
types_m = {
    'g':'gompertz',
    'l':'logistic',
    'v':'von_bert',
}
types_title = {
    'g':'Gompertz',
    'l':'Logistic',
    'v':'Von Bert',
}
types_param = {
    'g':'beta',
    'l':'r',
    'v':'kappa',
}
#_______________________________________________________________________________
def read_data(model_type):
    #Leer txt
    with open(f'{types_m[model_type]}_param.txt', 'r') as file:
        model_param = file.read()
    with open(f'{types_m[model_type]}_sigma.txt', 'r') as file:
        model_sigma = file.read()
    #Diccionario
    data = {
        'Params':model_param.split(),
        'Sigma':model_sigma.split(),
    }
    #DataFrame
    df = pd.DataFrame(data)
    df = df.astype(float)
    df = df.reset_index()
    df['index'] = df['index'] + 1
    df = df.rename(columns={'index':'x'})
    return df
#_______________________________________________________________________________
def plot_data(x, y,flg_sigma=True, type_sim='g', zoom=False): 
    # Crear la figura y los ejes
    fig, ax = plt.subplots()
    # Crear la gráfica de dispersión
    ax.scatter(x, y, s=0.1)
    # Añadir etiquetas y título
    ax.set_title(types_title[type_sim]) 
    ax.set_xlabel('Número de iteraciones')
    if zoom:
        q1 = y.quantile(0.99)
        q3 = y.quantile(0.01)
        # Calcular el rango intercuartílico (IQR)
        iqr = q3 - q1
        # Definir los límites inferior y superior
        l_inf = q1 - 1.5 * iqr
        l_sup = q3 + 1.5 * iqr
        # Aplicar zoom en el eje y
        ax.set_ylim(l_inf, l_sup)
    #Agregar linea horizontal
    y_ = list(y)
    ax.axhline(y_[-1], color='red',linewidth=0.8,label=f'Valor estimado: {y_[-1]}')
    if flg_sigma:
        ax.set_ylabel('sigma')
    else:
        ax.set_ylabel(types_param[type_sim])
    #Solo 3 decimales
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    return fig
#_______________________________________________________________________________
try:
    df = read_data('g')
    #
    plt_sigma = plot_data(df['x'],df['Sigma'],flg_sigma=True,type_sim='g',zoom=True)
    plt_param = plot_data(df['x'],df['Params'],flg_sigma=False ,type_sim='g',zoom=True)
    #
    plt_sigma.savefig(f'{types_m["g"]}_plt_sigma.png')
    plt_param.savefig(f'{types_m["g"]}_plt_{types_param["g"]}.png')
except:
    pass
#_______________________________________________________________________________
try:
    df = read_data('l')
    #
    plt_sigma = plot_data(df['x'],df['Sigma'],flg_sigma=True,type_sim='l',zoom=True)
    plt_param = plot_data(df['x'],df['Params'],flg_sigma=False ,type_sim='l',zoom=True)
    #
    plt_sigma.savefig(f'{types_m['l']}_plt_sigma.png')
    plt_param.savefig(f'{types_m['l']}_plt_{types_param["l"]}.png')
except:
    pass
#_______________________________________________________________________________
try:
    df = read_data('v')
    #
    plt_sigma = plot_data(df['x'],df['Sigma'],flg_sigma=True,type_sim='v',zoom=True)
    plt_param = plot_data(df['x'],df['Params'],flg_sigma=False ,type_sim='v',zoom=True)
    #
    plt_sigma.savefig(f'{types_m["v"]}_plt_sigma.png')
    plt_param.savefig(f'{types_m["v"]}_plt_{types_param["v"]}.png')
except:
    pass