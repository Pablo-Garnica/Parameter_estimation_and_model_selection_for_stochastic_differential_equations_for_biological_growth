module menu
    use rand_seed
    use simulation
    use mle
    use qua
    use aic
    use rand_seed
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
            param_recom = "(0,inf) valor recomendado 12.0"
            sigma_recom = "(0,inf) valor recomendado 12.0"
            delta_recom = "(0,inf) valor recomendado 12.0"
            xstart_recom = "(0,inf) valor recomendado 12.0"
            npoints_recom = "Natural mayor a 1 valor recomendado 100"
        else if (type_model.eq."l") then
            name_model = "Logistic"
            model_name_param = "r"
            param_recom = "[50,inf) valor recomendado 7.0"
            sigma_recom = "[50,inf) valor recomendado 7.0"
            delta_recom = "[50,inf) valor recomendado 7.0"
            xstart_recom = "[50,inf) valor recomendado 7.0"
            npoints_recom = "Natural mayor a 1 valor recomendado 100"
        else if (type_model.eq."v") then
            name_model = "Von Bert"
            model_name_param = "kappa"
            param_recom = "[10,inf) valor recomendado 0.203"
            sigma_recom = "[10,inf) valor recomendado 0.203"
            delta_recom = "[10,inf) valor recomendado 0.203"
            xstart_recom = "[10,inf) valor recomendado 0.203"
            npoints_recom = "Natural mayor a 1 valor recomendado 100"
        else
            print *, "Error"
        end if
    end subroutine
    subroutine menu_input_param(type_model,model_param,sigma,delta,xstart,npoints,linf)
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
        character, intent(in) :: type_model
        real*8, intent(out) :: model_param,sigma,delta,xstart
        integer, intent(out) :: npoints
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
        print *, "-------------------------------------------------------------------"
    end subroutine
    subroutine menu_result(type_model,model_param,sigma,delta,xstart,npoints,linf,paramhat,sigmahat,aic_param)
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
        !-------------------------------------------------------------------
        character, intent(in) :: type_model
        real*8, intent(in) :: model_param,sigma,delta,xstart
        integer, intent(in) :: npoints
        real*8, intent(in), optional :: linf
        real*8, intent(out) :: paramhat,sigmahat,aic_param
        real*8 :: path(npoints)
        character*20 :: name_model,model_name_param
        character*100 :: param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom
        call dict_models(type_model,name_model,model_name_param, &
        param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom)
        call SIM(type_model,model_param,sigma,delta,xstart,npoints,path,linf)
        call MLE_(type_model,npoints,path,delta,paramhat,linf)
        call Qua_Var(type_model,npoints,path,delta,sigmahat,linf)
        call AIC_(type_model,delta,paramhat,sigmahat,npoints,path,aic_param,linf)
        print *, "Quadratic variation : ",sigmahat
        print *, "MLE (Maximum Likelihood Estimator) : ",trim(model_name_param) ,"hat : ",paramhat
        print *, "AIC (Akaike information criterion) : ",aic_param
        print *, "-------------------------------------------------------------------"
        print *, "-------------------------------------------------------------------"
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
end module menu