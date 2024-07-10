module aic
    use integral
    implicit none
contains
    subroutine AIC_GOM(delta,b,sigma,npoints,path,aic)
        !-------------------------------------------------------------------
        !> \brief Calcula el valor del AIC (Akaike information criterion)
        !> para el modelo Gompertz
        ! 
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] b(real*8) ???
        !> \param[in] sigma(real*8) Valor la ????
        !> \param[in] npoints(integer) Tamaño de la simulación
        !> \param[in] path(real*8) Matriz donde se guarda la simulación
        !> \param[out] aic(real*8) Akaike information criterion
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: delta,b,sigma
        integer, intent(in) :: npoints
        real*8, intent(out) :: aic
        real*8 :: path(npoints),I1,I2,pathint1(npoints),pathint2(npoints)
        pathint1(:)=-(b*path(:)*LOG(path(:)))/((sigma**2)*(path(:)**2))
        pathint2(:)=((b**2)*(path(:)**2))/((sigma**2)*(path(:)**2))
        call ItoIntegrate(npoints,pathint1,path,I1)
        call Integrate(npoints,pathint2,delta,I2)
        aic=-2*(I1-I2/2)
        return
    end subroutine
    subroutine AIC_VON(delta,kappa,sigma,npoints,path,aic,linf)
        !-------------------------------------------------------------------
        !> \brief Calcula el valor del AIC (Akaike information criterion)
        !> para el modelo Von Bert
        !> 
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] kappa(real*8) ????
        !> \param[in] sigma(real*8) Valor la ????
        !> \param[in] npoints(integer) Tamaño de la simulación
        !> \param[in] path(real*8) Matriz donde se guarda la simulación
        !> \param[out] aic(real*8) Akaike information criterion
        !> \param[in] linf(real*8) Limite superior
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: delta,kappa,sigma
        integer, intent(in) :: npoints
        real*8, intent(out) :: aic
        real*8, intent(in), optional  :: linf
        real*8 :: path(npoints),I1,I2,pathint1(npoints),pathint2(npoints),inf
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        pathint1(:)=kappa/((sigma**2)*(inf-path(:)))
        pathint2(:)=(kappa**2)/(sigma**2)
        call ItoIntegrate(npoints,pathint1,path,I1)
        call Integrate(npoints,pathint2,delta,I2)
        aic=-2*(I1-I2/2)
        return
    end subroutine
    subroutine AIC_LOG(delta,r,sigma,npoints,path,aic)
        !-------------------------------------------------------------------
        !> \brief Calcula el valor del AIC (Akaike information criterion)
        !> para el modelo Logistico
        ! 
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] r(real*8) ????
        !> \param[in] sigma(real*8) Valor la ????
        !> \param[in] npoints(integer) Tamaño de la simulación
        !> \param[in] path(real*8) Matriz donde se guarda la simulación
        !> \param[out] aic(real*8) Akaike information criterion
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: delta,r,sigma
        integer, intent(in) :: npoints
        real*8, intent(out) :: aic
        real*8 :: path(npoints),I1,I2,pathint1(npoints),pathint2(npoints)
        pathint1(:)=-(r*path(:)*(1.0-path(:)))/((sigma**2)*(path(:)**2))
        pathint2(:)=((r**2)*((path(:)*(1.0-path(:)))**2))/((sigma**2)*(path(:)**2))
        call ItoIntegrate(npoints,pathint1,path,I1)
        call Integrate(npoints,pathint2,delta,I2)
        aic=-2*(I1-I2/2)
    return
    end subroutine
    subroutine AIC_(type_model,delta,param,sigma,npoints,path,aic_param,linf)
        !-------------------------------------------------------------------
        !> \brief Calcula el valor del AIC (Akaike information criterion)
        ! 
        !> \param[in] type(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] kappa(real*8) ????
        !> \param[in] param(real*8) Valor la ????
        !> \param[in] npoints(integer) Tamaño de la simulación
        !> \param[in] path(real*8) Matriz donde se guarda la simulación
        !> \param[out] aic(real*8) Akaike information criterion
        !> \param[in] linf(real*8) Limite superior
        !-------------------------------------------------------------------
        character, intent(in) :: type_model
        real*8, intent(in) :: delta,param,sigma
        integer, intent(in) :: npoints
        real*8, intent(out) :: aic_param
        real*8, intent(in), optional  :: linf
        real*8 :: path(npoints),inf
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        if (type_model.eq."g") then
            call AIC_GOM(delta,param,sigma,npoints,path,aic_param)
        else if (type_model.eq."l") then
            call AIC_LOG(delta,param,sigma,npoints,path,aic_param)
        else if (type_model.eq."v") then
            call AIC_VON(delta,param,sigma,npoints,path,aic_param,linf)
        else
            print *, "Error"
        end if
    return
    end subroutine
end module aic