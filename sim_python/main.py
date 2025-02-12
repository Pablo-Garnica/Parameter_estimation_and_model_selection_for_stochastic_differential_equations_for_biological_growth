'''
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
'''
from psm.flow import *
from psm.graph import *
#_________________________________________________________________
d_exe = input_general_param()
d_f = input_type_exe(d_exe)
d_export =  input_plt(d_exe,d_f)