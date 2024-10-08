module qua
    use integral
    implicit none
contains
    subroutine Qua_Var_G(npoints,path,delta,sigmahat)
        !-------------------------------------------------------------------
        !> \brief Calcula de Quadratic variation para el modelo Gompertz
        ! 
        !> \param[in] npoints(integer) Numero de observaciones
        !> \param[in] path(real*8) Observaciones de la SDE (stochastic 
        !> differential equation)
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[out] sigmahat(real*8) Estimador sigma
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: npoints
        real*8, intent(in) :: delta
        real*8, intent(out) :: sigmahat
        real*8 :: bhat,j,k,dem,num,h,path(npoints),lpath(npoints)
        lpath(:)=LOG(path(:))
        num=((npoints-1)*SUM(lpath(1:(npoints-1))*lpath(2:npoints))-SUM(lpath(2:npoints))*SUM(lpath(1:(npoints-1))))
        dem=((npoints-1)*SUM(lpath(1:(npoints-1))**2.0)-SUM(lpath(1:(npoints-1)))**2.0)
        j=num/dem
        k=(j*SUM(lpath(1:(npoints-1)))-SUM(lpath(2:npoints)))/(npoints-1)
        h=SUM((lpath(2:npoints)-j*lpath(1:npoints-1)+k)**2.0)/(npoints-1)
        bhat=-LOG(j)/delta
        sigmahat=SQRT((h*2.0*bhat)/(1.0-EXP(-2.0*bhat*delta)))
        return
    end subroutine
    subroutine Qua_Var_L(npoints,path,delta,sigmahat)
        !-------------------------------------------------------------------
        !> \brief Calcula de Quadratic variation para el modelo Lofistic
        ! 
        !> \param[in] npoints(integer) Numero de observaciones
        !> \param[in] path(real*8) Observaciones de la SDE (stochastic 
        !> differential equation)
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[out] sigmahat(real*8) Estimador sigma
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: npoints
        real*8, intent(in)  :: delta
        real*8, intent(out)  :: sigmahat
        real*8 :: path(npoints),num,dem
        integer :: i
        num=0.0
        dem=0.0
        do i=2,npoints
            num=num+(path(i)-path(i-1))**2
            dem=dem+path(i)**2+path(i-1)**2
        enddo
        sigmahat=SQRT(2.0*num/(dem*delta))
        return
    end subroutine
    subroutine Qua_Var_VB(npoints,path,delta,linf,sigmahat)
        !-------------------------------------------------------------------
        !> \brief Calcula de Quadratic variation para el modelo Lofistic
        ! 
        !> \param[in] npoints(integer) Numero de observaciones
        !> \param[in] path(real*8) Observaciones de la SDE (stochastic 
        !> differential equation)
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] linf(real*8) Limite superior
        !> \param[out] sigmahat(real*8) Estimador sigma
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: npoints
        real*8, intent(in)  :: delta,linf
        real*8, intent(out)  :: sigmahat
        real*8 :: path(npoints),num,dem,pathint(npoints)
        integer :: i
        num=0.0
        pathint(:)=(linf-path(:))**2
        call Integrate(npoints,pathint,delta,dem)
        do i=2,npoints
          num=num+(path(i)-path(i-1))**2
        enddo
        sigmahat=SQRT(num/dem)     
        return
      end subroutine      
      subroutine Qua_Var(type_model,npoints,path,delta,sigmahat,linf)
        !-------------------------------------------------------------------
        !> \brief Calcula de Quadratic variation para segun el modelo 
        !> seleccionado
        ! 
        !> \param[in] type_model(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[in] npoints(integer) Numero de observaciones
        !> \param[in] path(real*8) Observaciones de la SDE (stochastic 
        !> differential equation)
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[out] sigmahat(real*8) Estimador sigma
        !> \param[in] linf(real*8) Limite superior
        !-------------------------------------------------------------------
        implicit none
        character, intent(in) :: type_model
        integer, intent(in) :: npoints
        real*8 :: path(npoints)
        real*8, intent(in) :: delta
        real*8, intent(out) :: sigmahat
        real*8, intent(in) , optional :: linf
        real*8 :: inf
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        if (type_model.eq."g") then
            call Qua_Var_G(npoints,path,delta,sigmahat)
        else if (type_model.eq."l") then
            call Qua_Var_L(npoints,path,delta,sigmahat)
        else if (type_model.eq."v") then
            call Qua_Var_VB(npoints,path,delta,linf,sigmahat)
        else
            print *, "Error"
        end if
        return
    end subroutine
end module qua