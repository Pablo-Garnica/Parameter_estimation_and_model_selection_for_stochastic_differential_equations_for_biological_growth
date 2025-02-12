import numpy as np
from psm.sim_mod import *
from psm.aic_mod import *
from psm.mle_mod import *
from psm.qua_mod import *
#________________________________________________
def choose_data(array, n_obs):
    l_a = len(array)
    size_jump = (l_a - 1) // (n_obs - 1)
    labels = [size_jump * (i - 1) for i in range(1, n_obs + 1)]
    array_1 = [array[label] for label in labels]
    return np.array(array_1)
#________________________________________________
def em_path(type_model,param,sigma,delta,n_iter_choose,start,end,n_obs,l_inf=9999999999.0):
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
    array_final = np.array([])
    for n in range(0,len(array)-1):
        start_for = array[n]
        end_for = array[n+1]
        #
        arr_for = em_path(type_model,param,sigma,delta,n_iter_choose,start_for,end_for,n_obs,l_inf)
        if n == 0:
            array_final = arr_for
        else:
            array_final = np.vstack((array_final, arr_for))
    return array_final

# def em_array(array,type_dif,param_0,sigma_0,delta,n_obs):
#     n_iter = len(array)
#     n_iter_chose = (n_iter // n_obs)
#     array_final = np.array([])
#     print(n_iter_chose)
#     for n in range(0,len(array[:-1])):
#         start_f = array[n]
#         end_f = array[n+1]
#         print(f'{start_f} -> {end_f}')
#     # for x in array[:-1]:
#         arr_for_1 = sim(type_dif=type_dif,param=param_0,sigma=sigma_0,delta=delta,start=start_f,n_iter=n_iter_chose)
#         arr_for_2 = sim(type_dif=type_dif,param=-param_0,sigma=sigma_0,delta=delta,start=end_f,n_iter=n_iter_chose)
#         arr_for_2 = arr_for_2[::-1]
#         #Cruce entre trayectorias
#         arr_for = np.array([])
#         if np.max(arr_for_1)<=np.max(arr_for_2):
#             print(1,np.max(arr_for_1),np.max(arr_for_2))
#             for i in range(1, n_obs):
#                 # if max(arr_for_1[:i]) > arr_for_2[i]:
#                 if arr_for_1[i] <= arr_for_2[i]:
#                     arr_for = np.append(arr_for_1[:i], arr_for_2[i:])
#                     arr_for = np.append(np.array(start_f),arr_for)
#                     print(i,arr_for_1[i],arr_for_2[i])
#                     break
#             # elif max(arr_for_2[:i]) > arr_for_1[i]:
#         elif np.max(arr_for_1)>=np.max(arr_for_2):
#             print(2,np.max(arr_for_1),np.max(arr_for_2))
#             for i in range(1, n_obs):
#                 if arr_for_2[i] >= arr_for_1[i]:
#                     # print('-------------------- caso 2 --------------------')
#                     arr_for = np.append(arr_for_1[:i], arr_for_2[i:])
#                     arr_for = np.append(np.array(start_f),arr_for)
#                     print(i,arr_for_1[i],arr_for_2[i])
#                     # arr_for = np.concatenate(array_1[i:],arr_for_2[:i])
#                     break
#         else:
#             arr_for = arr_for_1
#             break
        # for i in range(1, n_iter_chose):
        #     # if max(arr_for_1[:i]) > arr_for_2[i]:
        #     if arr_for_1[i] >= arr_for_2[i]:
        #         arr_for = np.append(arr_for_1[:i], arr_for_2[i:])
        #         arr_for = np.append(np.array(start_f),arr_for)
        #         print('-------------------- caso 1 --------------------')
        #         break
        #     #
        #     # elif max(arr_for_2[:i]) > arr_for_1[i]:
        #     elif arr_for_2[i] >= arr_for_1[i]:
        #         print('-------------------- caso 2 --------------------')
        #         arr_for = np.append(arr_for_1[:i], arr_for_2[i:])
        #         arr_for = np.append(np.array(start_f),arr_for)
        #         # arr_for = np.concatenate(arr_for_1[i:],arr_for_2[:i])
        #         break
        # if len(arr_for) == 0:
        #     arr_for = arr_for_1
        #     print('-------------------- caso raro --------------------')
        #Array final
#         array_final = np.append(array_final,arr_for)
#     array_final = array_final[:n_iter]
#     return array_final

# def em_mc(type_dif, param_0, sigma_0,start,n_obs, delta, n_iter, n_reply, l_inf=9999999999.0):
#     arr_qua = np.array([])
#     arr_mle = np.array([])
#     n_iter_chose = n_iter // n_obs
#     print(f'n_iter_chose:{n_iter_chose}')
#     md = n_reply//10
#     for i in range(0,n_reply):
#         if i%md ==0 :
#             msg =f'Avance {np.round((i / n_reply)*100,1)}%'
#             print(msg)
#         array = sim(type_dif=type_dif,param=param_0,sigma=sigma_0,delta=delta,start=start,n_iter=n_iter)
#         array_ch = choose_data(array,delta=delta,nobs=n_obs)
#         array_final = em_array(array_ch,type_dif,param_0,sigma_0,delta,n_obs)
#         #
#         qua = qua_var(type_dif, array_final, delta, l_inf=l_inf)
#         arr_qua =  np.append(arr_qua,qua) 
#         mle = mle_(type_dif, array, delta, l_inf=l_inf)
#         arr_mle = np.append(arr_mle, mle)
#     d_r = {
#         'qua':arr_qua,
#         'mle':arr_mle,
#     }
#     return d_r
# #________________________________________________
# def diffusion_bridge(type_dif, param, sigma, delta, start, end,n_iter_chose, l_inf=9999999999.0):
#     #
#     while True:
#         #
#         array_1 = sim(type_dif,param,sigma,delta,start,n_iter_chose,l_inf)
#         array_2 = sim(type_dif,-param,sigma,delta,end,n_iter_chose,l_inf)
#         #
#         # return array_1
#         min_1 = np.min(array_1)
#         min_2 = np.min(array_2)
#         #
#         max_1 = np.max(array_1)
#         max_2 = np.max(array_2)
#         #
#         if (min_1>max_2) | (min_2>max_1):
#             break
# #     # #
#     for i in range(1, n_iter_chose):
#         #
#         if max(array_1[:i]) > array_2[i]:
#             array_3 = np.concatenate((array_1[:i], array_2[i:]))
#             return array_3
#         #
#         elif max(array_2[:i]) > array_1[i]:
#             array_3 = np.concatenate((array_2[:i], array_1[i:]))
#             return array_3 
#         else:
#             return array_1
# #________________________________________________
# def complete_bridge_fn(type_dif, param, sigma, n_iter_chose, delta_bridge, n_obs, array, l_inf=9999999999.0):
#     # complete_bridge = np.zeros((n_iter_chose + 1) * (n_obs - 1) + 1)
#     array_complete = np.array(array[0])
#     for i in range(1, n_obs):
#         start = array[i - 1]
#         end = array[i]
#         array_bridge = diffusion_bridge(type_dif, param, sigma, delta_bridge, start, end,n_iter_chose, l_inf)
#         array_complete = np.append(start,array_bridge)
#     array_complete = np.append(array_bridge,end)
#     return array_complete
# #________________________________________________
# def em_mc(
#         type_dif, 
#         param_0, 
#         sigma_0, 
#         array, 
#         n_obs, 
#         delta, 
#         n_iter_chose, 
#         n_reply, 
#         n_error=1, 
#         l_inf=9999999999.0
#     ):
#     array_ = np.log(array) if type_dif == "g" else array
#     #
#     param_vec = np.array([param_0])
#     sigma_vec = np.array([sigma_0])
#     #
#     delta_bridge = delta
#     #
    
#     for i in range(1, n_reply):
#         bhs = np.array([])
#         shs = np.array([])
#         for j in range(n_error):
#             array_brd_com = complete_bridge_fn(type_dif, param_vec[i-1], sigma_vec[i-1], n_iter_chose, delta_bridge, n_obs, array_, l_inf)
#             mle = mle_(type_dif,array_brd_com, delta_bridge, l_inf)
#             qua = qua_var(type_dif, array_brd_com, delta_bridge, l_inf)
#             bhs = np.append(bhs,mle)
#             shs = np.append(shs,qua)
#             if j ==0:
#                 array_conat = array_brd_com
#             else:
#                 array_conat = np.vstack(array_conat, array_brd_com)
#         sigma_vec = np.append(sigma_vec,np.mean(shs))
#         param_vec = np.append(param_vec,np.mean(bhs))
#     return param_vec, sigma_vec, array_conat