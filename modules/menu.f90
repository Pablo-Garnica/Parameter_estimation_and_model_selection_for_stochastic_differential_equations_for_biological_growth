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
        character, intent(out) :: type_model
        print *, "-------------------------------------------------------------------"
        print *, "-------------------------------------------------------------------"
        print *, "Introduce el modelo que quieres calcular:"
        print *, "g : Gompertz"
        print *, "l : Logistic"
        print *, "v : Von Bert"
        read *, type_model
    end subroutine

    subroutine dict_models(type_model,name_model,model_name_param)
        character, intent(in) :: type_model
        character*20, intent(out) :: name_model,model_name_param
        if (type_model.eq."g") then
            name_model = "Gompertz"
            model_name_param = "beta"
        else if (type_model.eq."l") then
            name_model = "Logistic"
            model_name_param = "r"
        else if (type_model.eq."v") then
            name_model = "Von Bert"
            model_name_param = "kappa"
        else
            print *, "Error"
        end if
    end subroutine
    subroutine menu_input_param(type_model,model_param,sigma,delta,xstart,npoints,linf)
        character, intent(in) :: type_model
        real*8, intent(out) :: model_param,sigma,delta,xstart
        integer, intent(out) :: npoints
        real*8, intent(out), optional :: linf
        character*20 :: name_model,model_name_param
        call dict_models(type_model,name_model,model_name_param)
        print *, "-------------------------------------------------------------------"
        print *, "Seleccionaste el modelo ", trim(name_model)
        print *, "Para este modelo se requieren los siguientes argumentos"
        print *, "Introduce el valor ",trim(model_name_param)
        read *, model_param
        print *, "Introduce el valor sigma"
        read *, sigma
        print *, "Introduce el valor delta"
        read *, delta
        print *, "Introduce el valor inicial de la simulación"
        read *, xstart
        print *, "Introduce el numero de valores para la simulacion"
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
    subroutine menu_result(type_model,model_param,sigma,delta,xstart,npoints,linf,paramhat,sigmahat)
        character, intent(in) :: type_model
        real*8, intent(in) :: model_param,sigma,delta,xstart
        integer, intent(in) :: npoints
        real*8, intent(in), optional :: linf
        real*8, intent(out) :: paramhat,sigmahat
        real*8 :: path(npoints)
        character*20 :: name_model,model_name_param
        call dict_models(type_model,name_model,model_name_param)
        call SIM(type_model,model_param,sigma,delta,xstart,npoints,path,linf)
        call MLE_(type_model,npoints,path,delta,paramhat,linf)
        call Qua_Var(type_model,npoints,path,delta,sigmahat,linf)
        print *, "Quadratic variation : ",sigmahat
        print *, "MLE (Maximum Likelihood Estimator) : ",trim(model_name_param) ,"hat : ",paramhat
        print *, "-------------------------------------------------------------------"
        print *, "-------------------------------------------------------------------"
    end subroutine
end module menu