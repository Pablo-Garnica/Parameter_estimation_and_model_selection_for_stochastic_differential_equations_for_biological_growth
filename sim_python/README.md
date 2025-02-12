# Parameter estimation and model selection for stochastic differential equations for biological growth
## Creación entorno virtual
En la carpeta del proyecto ejecutar
```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

## Ejecución codigo
```bash
.venv\Scripts\activate
python main.py
```

## Output ejecución normal
Si se ejecuta la aplicacion en tipo n (normal), se puede exportar la información de la simulación en json y la grafica
se exporta en la misma ruta la cual se ejecuta el main.py
## Configuraciones
### Normal
- type_dif = 'v'
- param = 0.6
- sigma = 0.001
- delta = 0.001
- start = 0.01
- n_iter = 10001
- n_reply = 5
### Trayectorias
- type_dif = 'g'
- param = 0.6
- sigma = 0.01
- delta = 0.001
- start = 0.01
- n_iter = 10001
- n_reply = 100
- n_step = 10000 | 5000 | 1000
### EM
- type_dif = 'l'
- param = 0.6
- sigma = 0.1
- delta = 0.001
- start = 0.01
- n_iter = 101
- n_reply = 3
- n_obs = 10