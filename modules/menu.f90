module menu
    use usefull_func
    use rand_seed
    use simulation
    use mle
    use qua
    use aic
    use rand_seed
    use em
    implicit none
    contains
    subroutine menu_random_seed()
        !-------------------------------------------------------------------
        !> \brief En la aplicación de consola da la opción de seleccionar 
        !> y fijar la semilla aleatorea
        !>
        !-------------------------------------------------------------------
        integer :: seed
        character :: option
        print *, "-------------------------------------------------------------------"
        print *, "-------------------------------------------------------------------"
        print *, "Esta de acuerdo con el la semilla aleatorea sea 31416"
        print *, "'y' para estar de acuerdo cualquier otro valor para cambiar"
        read *, option
        if (option.eq."y") then
            seed = 31416
        else
            print *, "Introduce el valor de la semilla aleatoria"
            read *, seed
        end if
        call random_seed_set_up(seed)
    end subroutine
    subroutine menu_linf(linf)
        !-------------------------------------------------------------------
        !> \brief Pregunta si esta bien el valor por de default, en caso de
        !> no estar de acuerdo cambiarlo
        !>
        !> \param[in] linf(real*8) Limite superior
        !> estocasticas
        !-------------------------------------------------------------------
        real*8, intent(out) :: linf
        character :: option
        print *, "Esta de acuerdo que el limite superior sea 999999999999999999999999999999.00"
        print *, "'y' para estar de acuerdo cualquier otro valor para cambiar"
        read *, option
        if (option.eq."y") then
            linf = 999999999999999999999999999999.00
        else
            print *, "Introduce el valor del limite superior"
            read *, linf
        end if
    end subroutine
    subroutine menu_model_type(type_model)
        !-------------------------------------------------------------------
        !> \brief Escoge el modelo que se quiere simular y evaluar
        !>
        !> \param[in] type_model(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !-------------------------------------------------------------------
        character, intent(out) :: type_model
        print *, "-------------------------------------------------------------------"
        print *, "-------------------------------------------------------------------"
        print *, "Introduce el modelo que quieres calcular:"
        print *, "g : Gompertz"
        print *, "l : Logistic"
        print *, "v : Von Bert"
        read *, type_model
    end subroutine
    subroutine dict_models(type_model,name_model,model_name_param,&
                            param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom)
        !-------------------------------------------------------------------
        !> \brief Diccionario de valores condicionales a partir de la
        !> selección del modelo
        !>
        !> \param[in] type_model(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[out] name_model(character*20) Nombre del modelo
        !> \param[out] model_name_param(character*20) Nombre del parametro 
        !> del modelo
        !> \param[out] param_recom(character*100) Intervalo o dominio para
        !> seleccionar el valor del parametro y un parametro recomendado
        !> \param[out] sigma_recom(character*100) Intervalo o dominio para
        !> seleccionar el valor sigma y un parametro recomendado
        !> \param[out] delta_recom(character*100) Intervalo o dominio para
        !> seleccionar el valor delta y un parametro recomendado
        !> \param[out] xstart_recom(character*100) Intervalo o dominio para
        !> seleccionar el valor de inicio del la simulació y un parametro 
        ! recomendado
        !> \param[out] npoints_recom(character*100) Numero recomendado de
        !> simulaciones
        !-------------------------------------------------------------------
        character, intent(in) :: type_model
        character*20, intent(out) :: name_model,model_name_param
        character*100, intent(out) :: param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom
        if (type_model.eq."g") then
            name_model = "Gompertz"
            model_name_param = "beta"
            sigma_recom = "Valor recomendado 0.1"
            param_recom = "Valor recomendado 0.6"
            delta_recom = "Valor recomendado 0.001"
            xstart_recom = "Valor recomendado 0.01"
            npoints_recom = "Valor recomendado 10,000"
        else if (type_model.eq."l") then
            name_model = "Logistic"
            model_name_param = "r"
            sigma_recom = "Valor recomendado 0.1"
            param_recom = "Valor recomendado 0.6"
            delta_recom = "Valor recomendado 0.001"
            xstart_recom = "Valor recomendado 0.01"
            npoints_recom = "Valor recomendado 10,000"
        else if (type_model.eq."v") then
            name_model = "Von Bert"
            model_name_param = "kappa"
            sigma_recom = "Valor recomendado 0.1"
            param_recom = "Valor recomendado 0.6"
            delta_recom = "Valor recomendado 0.001"
            xstart_recom = "Valor recomendado 0.01"
            npoints_recom = "Valor recomendado 10,000"
        else
            print *, "Error"
        end if
    end subroutine


    subroutine menu_type_execution(type_execution)
        !-------------------------------------------------------------------
        !> \brief Nos da el input si se tiene que ejecutar el algoritmo em
        ! 
        !> \param[in] type_execution(character) Si el parametro es "y" entonces
        !> se ejecuta el algoritmo em, en otro caso no se ejecuta
        !-------------------------------------------------------------------
        character, intent(out) :: type_execution
        print *, "-------------------------------------------------------------------"
        print *, "Que tipo de ejecición desea aplicar"
        print *, "'n' : Ejecución normal del modelo"
        print *, "'t' : Ejecución del modelo para seleccionar trayectoria"
        print *, "'e' : Ejecución ocupando el algoritmo em"
        read *, type_execution            
    end subroutine







    subroutine menu_input_param(type_model,type_execution,model_param,sigma,delta,xstart,npoints,&
            nobs,niter,nsteps,nmc,linf)
        !-------------------------------------------------------------------
        !> \brief La aplicación de consola hace que el usuario agregue los 
        !> inputs segun el modelo
        !>
        !> \param[in] type_model(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[out] param(real*8) Valor de ??? si type es:
        !>   "v" : param es en realidad el valor kappa
        !>   "g" : param es beta
        !>   "l" : param es en realidad el valor r
        !> \param[out] sigma(real*8) ??
        !> \param[out] delta(real*8) Incremento del proceso de Wiener
        !> \param[out] xstart(real*8) Valor inicial en el que se evalua el 
        !> modelo
        !> \param[out] npoints(integer) Tamaño de la simulación
        !> \param[out] linf[optional](real*8) Limite superior
        !-------------------------------------------------------------------
        character, intent(in) :: type_model,type_execution
        real*8, intent(out) :: model_param,sigma,delta,xstart
        integer, intent(out) :: npoints,nobs,niter,nsteps,nmc
        real*8, intent(out), optional :: linf
        character*20 :: name_model,model_name_param
        character*100 :: param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom
        call dict_models(type_model,name_model,model_name_param, &
        param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom)
        print *, "-------------------------------------------------------------------"
        print *, "Seleccionaste el modelo ", trim(name_model)
        print *, "Para este modelo se requieren los siguientes argumentos"
        print *, "Introduce el valor ",trim(model_name_param)," en " ,trim(param_recom)
        read *, model_param
        print *, "Introduce el valor sigma en ",trim(sigma_recom)
        read *, sigma
        print *, "Introduce el valor delta en ",trim(delta_recom)
        read *, delta
        print *, "Introduce el valor inicial de la simulación en ",trim(xstart_recom)
        read *, xstart
        print *, "Introduce el numero de valores para la simulacion en ",trim(npoints_recom)
        read *, npoints
        if (type_model.eq."v") then
            call menu_linf(linf)
        end if
        print *, "-------------------------------------------------------------------"
        print *, "Los parametros son :"
        print *, trim(model_name_param), " : ",model_param
        print *, "sigma : ",sigma
        print *, "delta : ",delta
        print *, "Valor inicial : ",xstart
        print *, "Numero de simulaciones : ",npoints
        if (type_model.eq."v") then
            print *, "Limite superior : ",linf
        end if
        if (type_execution.eq."t") then
            print *, "-------------------------------------------------------------------"
            print *, "Introduce el numero de iteraciones del modelo para seleccionar"
            print *, "la trayectoria" 
            read*,niter
            print *, "-------------------------------------------------------------------"
            print *, "Introduce el numero de pasos para ejecutar el modelo para "
            print *, "seleccionar la trayectoria" 
            read*,nsteps
        end if
        if (type_execution.eq."e") then
            print *, "-------------------------------------------------------------------"
            print *, "Introduce el numero de para realizar la muestra del y ejecutar "
            print *, "el modelo em"
            read*,nobs
            print *, "-------------------------------------------------------------------"
            print *, "Introduce el numero de iteraciones del modelo ejecutar el modelo"
            print *, "em" 
            read*,niter
            print *, "-------------------------------------------------------------------"
            print *, "Introduce el numero de pasos del modelo ejecutar el modelo"
            print *, "em" 
            read*,nsteps
            print *, "-------------------------------------------------------------------"
            print *, "Introduce el nmc del modelo ejecutar el modelo em"
            read*,nmc
        end if
        print *, "-------------------------------------------------------------------"
    end subroutine
    subroutine menu_result(type_model,type_execution,model_param,sigma,delta,xstart,npoints,nobs,&
            niter,nsteps,nmc,linf,paramhat,sigmahat,aic_param,path)
        !-------------------------------------------------------------------
        !> \brief Calcula y muestra los valores MLE (Maximum Likelihood 
        !> Estimator), AIC (Akaike information criterion), Quadratic 
        !> variation y paramhat que es ????
        !>
        !> \param[in] type_model(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[in] type_execution(character) Debe de estar en los siguientes
        !>   "n" : Ejecución normal del modelo
        !>   "t" : Ejecución del modelo para seleccionar trayectoria
        !>   "e" : Ejecución ocupando el alhoritmo em
        !> \param[in] model_param(real*8) ????
        !> \param[in] sigma(real*8) Valor de ?????
        !> \param[out] sigma(real*8) ????
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] xstart(real*8) Valor inicial en el que se evalua el 
        !> modelo
        !> \param[in] npoints(integer) Tamaño de la simulación
        !> \param[in] linf[optional](real*8) Limite superior
        !> \param[out] paramhat(real*8) ????
        !> \param[out] sigmahat(real*8) ????
        !> \param[out] aic_param(real*8) AIC (Akaike information criterion)
        !> \param[out] Path de simulación
        !-------------------------------------------------------------------
        character, intent(in) :: type_model,type_execution
        real*8, intent(in) :: model_param,sigma,delta,xstart
        integer, intent(in) :: npoints,nobs,niter,nsteps,nmc
        real*8, intent(in), optional :: linf
        real*8, intent(out) :: paramhat,sigmahat,aic_param
        real*8, intent(out) :: path(npoints)
        character*20 :: name_model,model_name_param
        character*100 :: param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom
        real*8 :: path_iter(niter,npoints),path_out(nsteps),path_out_e(nobs)
        real*8 :: delta_bridge
        real*8 :: sigma_vec(nmc),beta_vec(nmc)

        call dict_models(type_model,name_model,model_name_param, &
        param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom)


        if (type_execution.eq."n") then
            call SIM(type_model,model_param,sigma,delta,xstart,npoints,path,linf)
            call MLE_(type_model,npoints,path,delta,paramhat,linf)
            call Qua_Var(type_model,npoints,path,delta,sigmahat,linf)
            call AIC_(type_model,delta,paramhat,sigmahat,npoints,path,aic_param,linf)
            print *, "Quadratic variation : ",sigmahat
            print *, "MLE (Maximum Likelihood Estimator) : ",trim(model_name_param) ,"hat : ",paramhat
            print *, "AIC (Akaike information criterion) : ",aic_param
            print *, "-------------------------------------------------------------------"
            print *, "-------------------------------------------------------------------"

        else if (type_execution.eq."t") then
            ! Simulación de elección de trayectorias
            ! Necesita niter y nsteps
            call SIM_ITER(type_model,model_param,sigma,delta,xstart,npoints,niter,path_iter,linf)
            call SIM_CHOOSE(path_iter,npoints,niter,nsteps,path_out,linf)
            call MLE_(type_model,nsteps,path_out,delta,paramhat,linf)
            call Qua_Var(type_model,nsteps,path_out,delta,sigmahat,linf)
            call AIC_(type_model,delta,paramhat,sigmahat,nsteps,path_out,aic_param,linf)
            print *, "Quadratic variation : ",sigmahat
            print *, "MLE (Maximum Likelihood Estimator) : ",trim(model_name_param) ,"hat : ",paramhat
            print *, "AIC (Akaike information criterion) : ",aic_param
            print *, "-------------------------------------------------------------------"
            print *, "-------------------------------------------------------------------"  
        else if (type_execution.eq."e") then
            print*,'En construcción em'

            call SIM(type_model,model_param,sigma,delta,xstart,npoints,path,linf)
            
            call choose_data(delta,npoints,path,nobs,delta_bridge,path_out_e)
            print*,'choose_data',path_out_e

            call EM_MC(type_model,model_param,sigma,path_out_e,nobs,delta_bridge,niter,nmc,nsteps,sigma_vec,beta_vec,linf)

            print *, "Quadratic variation : ",sigma_vec
            print *, "MLE (Maximum Likelihood Estimator) : ",trim(model_name_param) ,"hat : ",beta_vec

            print *, "-------------------------------------------------------------------"
            print *, "-------------------------------------------------------------------"  
        else
            print *, "Error"
        end if
    end subroutine
    subroutine date_char(date)
        !-------------------------------------------------------------------
        !> \brief Calcula la fecha y hora en el formato YYYYMMDDHHMMSS
        !>
        !> \param[out] date(character*14) Fecha y hora en el formato 
        !> YYYYMMDDHHMMSS
        !-------------------------------------------------------------------
        character*14, intent(out) :: date
        integer :: date_array(8)
        character*4 :: year
        character*2 :: month, day, hour, minute, second
        ! Llamar a la subrutina DATE_AND_TIME para obtener la fecha y la hora actuales
        call date_and_time(values=date_array)
        ! Extraer la fecha y la hora del array date_array y convertir a cadenas
        write(year, '(I4)') date_array(1)
        write(month, '(I2.2)') date_array(2)
        write(day, '(I2.2)') date_array(3)
        write(hour, '(I2.2)') date_array(5)
        write(minute, '(I2.2)') date_array(6)
        write(second, '(I2.2)') date_array(7)
        date = trim(year) //  trim(month) // trim(day)// trim(hour)// trim(minute)// trim(second)
    end subroutine
    subroutine menu_export(type_model,path,npoints)
        !-------------------------------------------------------------------
        !> \brief Menu que pregunta si quiere exportar la simulación
        !> Si se exporta se hace con el nombre 
        !> [nombre_modelo][fecha_hora].txt
        !>
        !> \param[in] type_model(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[in] path(real*8) Matriz de la simulación que se quiere
        !> exportar
        !> \param[in] npoints(integer) Tamaño de la simulación
        !-------------------------------------------------------------------
        character, intent(in) :: type_model
        integer, intent(in) :: npoints
        real*8, intent(in) :: path(npoints)
        character :: option
        character*14 :: date
        character*50 :: name_file
        call date_char(date)
        if (type_model.eq."g") then
            name_file = "gompertz" //  trim(date) // ".txt"
        else if (type_model.eq."l") then
            name_file = "logistic" //  trim(date) // ".txt"
        else if (type_model.eq."v") then
            name_file = "von_bert" //  trim(date) // ".txt"
        else
            print *, "Error"
        end if
        print *, "Quiere exportar la simulacion"
        print *, "'y' para estar exportar, cualquier otro valor para NO exportar"
        read*,option
        if (option.eq."y") then
            open(1,file=name_file)
            write(1,*) path
            endfile(1)
            close(1)
            print*,"Se exporto exitosamente el archivo: ",trim(name_file)
            print *, "-------------------------------------------------------------------"
        else
            print *, "-------------------------------------------------------------------"
        end if
    end subroutine
    subroutine menu_export_params(type_model,path_sim,npoints,param_1,sigma_1,delta,linf)
        !-------------------------------------------------------------------
        !> \brief Menu que pregunta si quiere exportar la simulación
        !> Si se exporta se hace con el nombre 
        !> [nombre_modelo][fecha_hora].txt
        !>
        !> \param[in] type_model(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[in] path(real*8) Matriz de la simulación que se quiere
        !> exportar
        !> \param[in] npoints(integer) Tamaño de la simulación
        !-------------------------------------------------------------------
        character, intent(in) :: type_model
        integer, intent(in) :: npoints
        real*8, intent(in) :: path_sim(npoints),param_1,sigma_1,delta
        real*8 :: path_sigma(npoints),path_param(npoints)
        character :: option
        character*14 :: date
        character*50 :: name_file,name_file_exp
        integer :: i
        real*8, intent(in) , optional :: linf
        real*8 inf, time
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        print *, "Quiere exportar la simulacion"
        print *, "'y' para estar exportar, cualquier otro valor para NO exportar"
        read*,option
        if (option.eq."y") then

            path_sigma(1) = param_1
            path_param(1) = sigma_1
            do i=2,npoints
                call MLE_(type_model,i,path_sim(1:i),delta,path_param(i),inf)
                call Qua_Var(type_model,i,path_sim(1:i),delta,path_sigma(i),inf)
            end do

            call date_char(date)
            if (type_model.eq."g") then
                name_file = "gompertz_"
            else if (type_model.eq."l") then
                name_file = "logistic_"
            else if (type_model.eq."v") then
                name_file = "von_bert_"
            else
                print *, "Error"
            end if
            name_file_exp = trim(name_file) // "param.txt"
            open(1,file=name_file_exp)
            write(1,*) path_param
            endfile(1)
            close(1)
            print*,"Se exporto exitosamente el archivo: ",trim(name_file)
            print *, "-------------------------------------------------------------------"
            name_file_exp = trim(name_file) // "sigma.txt"
            open(2,file=name_file_exp)
            write(2,*) path_sigma
            endfile(2)
            close(2)

            name_file_exp = trim(name_file) // "sigma_input.txt"
            open(3,file=name_file_exp)
            write(3,*) sigma_1
            endfile(3)
            close(3)

            name_file_exp = trim(name_file) // "param_input.txt"
            open(4,file=name_file_exp)
            write(4,*) param_1
            endfile(4)
            close(4)
            time = npoints * delta
            name_file_exp = trim(name_file) // "time.txt"
            open(5,file=name_file_exp)
            write(5,*) time
            endfile(5)
            close(5)
        else
            print *, "-------------------------------------------------------------------"
        end if
    end subroutine
end module menu