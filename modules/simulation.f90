module simulation
    implicit none
    contains
    function unif() result(y)
        !-------------------------------------------------------------------
        !> \brief Calcula un numero aleatorio con distibucion uniforme [0,1]
        !>
        !> \return y(real*8) numero aleatorio con distibucion uniforme [0,1]
        !-------------------------------------------------------------------
        implicit none
        real*8 :: y
        call random_number(y)
        return
    end function unif
    function boxmuller() result(y)
        !-------------------------------------------------------------------
        !> \brief Calcula un numero aleatorio con distibucion normal con 
        ! media 0 y desviación estándar 1 atraves del metodo boxmuller
        !>
        !> \return y(real*8) Numero aleatorio con distibucion normal con 
        ! media 0 y desviación estándar 1
        !-------------------------------------------------------------------
        integer iset
        real*8 fac,gset,rsq,v1,v2,y
        save iset,gset
        data iset/0/
        1 if (iset.eq.0) then
            v1=2.0*unif()-1.0
            v2=2.0*unif()-1.0
            rsq=v1**2+v2**2
        if (rsq.ge.1..or.rsq.eq.0) goto 1
            fac=sqrt(-2.*log(rsq)/rsq)
            gset=v1*fac
            y=v2*fac
            iset=1
        else
            y=gset
            iset=0
        endif
        return 
    end function boxmuller
    subroutine normalvar(x)
        !-------------------------------------------------------------------
        !> \brief Ejecuta la función boxmuller y guarda el resultado en el
        ! parametro x
        !>
        !> \param[out] x(real*8) Numero aleatorio con distibucion normal con 
        ! media 0 y desviación estándar 1
        !-------------------------------------------------------------------
        real*8, intent(out) :: x
        x=boxmuller()
        return
    end subroutine
    subroutine BrownianStep(delta,startx,endx)
        !-------------------------------------------------------------------
        !> \brief Simulación del movimiento Browniano en el tiempo delta e
        ! >incio en startx
        !>
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] startx(real*8) Valor inicial del movimiento Browniano
        !> \param[out] endx(real*8) Simulación en el tiempo delta
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: delta,startx
        real*8, intent(out) :: endx
        real*8 :: var
        call normalvar(var) 
        endx=startx+sqrt(delta)*var
        return
    end subroutine
    subroutine MilsteinStep(delta,startx,drix,difx,sigma,endx)
        !-------------------------------------------------------------------
        !> \brief Aplica el metodo Milstein para resolver ecuaciones 
        !> estocasticas
        !>
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] startx(real*8) Valor inicial del movimiento Browniano
        !> \param[in] drix(real*8) Valor de drift
        !> \param[in] difx(real*8) Valor de diffusion
        !> \param[in] sigma(real*8) Valor la ????
        !> \param[out] endx(real*8) Simulación de resolucion de ecuaciones 
        !> estocasticas
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: delta,startx,drix,difx,sigma
        real*8, intent(out) :: endx
        real*8 :: W,xx
        xx=0.0
        call BrownianStep(delta,xx,W)
        endx=startx+drix*delta+difx*W+(1.0/2.0)*difx*sigma*(W**2-delta)
        return
    end subroutine
    subroutine DiffusionParameter_L_G(sigma,x,y)
        !-------------------------------------------------------------------
        !> \brief Calcula el parametro diffusion para los modelos Gompertz y
        !> Logistic
        ! 
        !> \param[in] sigma(real*8) Valor de ???
        !> \param[in] x(real*8) Valor en el que se evalua la diffusion
        !> \param[out] y(real*8) Valor del parametro de diffusion
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: sigma,x
        real*8, intent(out) :: y
        y=sigma*x
        return
    end subroutine
    subroutine DiffusionParameter_V_B(sigma,linf,x,y)
        !-------------------------------------------------------------------
        !> \brief Calcula el parametro diffusion para el modelo Von Bert
        ! 
        !> \param[in] sigma(real*8) Valor de ???
        !> \param[in] x(real*8) Valor en el que se evalua la diffusion
        !> \param[in] linf(real*8) Limite superior
        !> \param[out] y(real*8) Valor del parametro de diffusion
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: sigma,x,linf
        real*8, intent(out) :: y
        y=sigma*(linf-x)
        return
      end subroutine
    subroutine DiffusionParam(type,sigma,x,y,linf)
        !-------------------------------------------------------------------
        !> \brief Calcula el parametro diffusion para el modelo Von Bert,
        !> Gompertz y Logistic
        ! 
        !> \param[in] type(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[in] sigma(real*8) Valor de ???
        !> \param[in] x(real*8) Valor en el que se evalua la diffusion
        !> \param[out] y(real*8) Valor del parametro de diffusion
        !> \param[in][optional] linf(real*8) Limite superior
        !-------------------------------------------------------------------
        implicit none    
        character, intent(in) :: type
        real*8, intent(in) :: sigma,x
        real*8, intent(in), optional :: linf
        real*8, intent(out) :: y
        real*8 :: inf
        !
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        !
        if ((type.eq."g").or.(type.eq."l")) then
            y=sigma*x
        else if (type.eq."v") then
            y=sigma*(inf-x)
        else
            print*,"Error"
        endif
        return
    end subroutine
    subroutine DriftParameter_G(beta,x,y)
        !-------------------------------------------------------------------
        !> \brief Calcula el parametro drift para el modelo Gompertz
        ! 
        !> \param[in] beta(real*8) Valor de ???
        !> \param[in] x(real*8) Valor en el que se evalua la drift
        !> \param[out] y(real*8) Valor del parametro de drift
        !-------------------------------------------------------------------
        implicit none 
        real*8, intent(in) :: beta,x
        real*8, intent(out) :: y
        y=-beta*x*log(x)
        return
    end subroutine
    subroutine DriftParameter_L(r,x,y)
        !-------------------------------------------------------------------
        !> \brief Calcula el parametro drift para el modelo Logistic
        ! 
        !> \param[in] r(real*8) Valor de ???
        !> \param[in] x(real*8) Valor en el que se evalua la drift
        !> \param[out] y(real*8) Valor del parametro de drift
        !-------------------------------------------------------------------
        implicit none 
        real*8, intent(in) :: r,x
        real*8, intent(out) :: y
        y=r*x*(1.0-x)
        return
      end subroutine
    subroutine DriftParameter_V(kappa,linf,x,y)
        !-------------------------------------------------------------------
        !> \brief Calcula el parametro drift para el modelo Von Bert
        ! 
        !> \param[in] kappa(real*8) Valor de ???
        !> \param[in] linf(real*8) Limite superior
        !> \param[in] x(real*8) Valor en el que se evalua la drift
        !> \param[out] y(real*8) Valor del parametro de drift
        !-------------------------------------------------------------------
        implicit none 
        real*8, intent(in) :: kappa, linf, x
        real*8, intent(out) :: y
        y=kappa*(linf-x)
        return
    end subroutine
    subroutine DriftParam(type,beta,x,y,linf)
        !-------------------------------------------------------------------
        !> \brief Calcula el parametro drift para el modelo Von Bert
        !> Gompertz y Logistic
        ! 
        !> \param[in] type(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[in] beta(real*8) Valor de ??? si type es:
        !>   "v" : beta es en realidad el valor kappa
        !>   "g" : beta es bet 
        !>   "l" : beta es en realidad el valor r
        !> \param[in] x(real*8) Valor en el que se evalua la drift
        !> \param[out] y(real*8) Valor del parametro de drift
        !> \param[in][optional] linf(real*8) Limite superior
        !-------------------------------------------------------------------
        implicit none 
        real*8 :: beta,x,y,inf
        real*8 , optional :: linf
        character :: type
        !
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        !
        if (type.eq."g") then
            y=-beta*x*log(x)
        else if (type.eq."l") then
            y=beta*x*(1.0-x)
        else if (type.eq."v") then
            y=beta*(inf-x)
        else
            print*,"Error"
        endif
        return
    end subroutine
    subroutine SIM(type,param,sigma,delta,x,npoints,path,linf)
        !-------------------------------------------------------------------
        !> \brief Calcula la resolución de las ecuaciones estocasticas de 
        !> los modelos
        ! 
        !> \param[in] type(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[in] param(real*8) Valor de ??? si type es:
        !>   "v" : param es en realidad el valor kappa
        !>   "g" : param es beta
        !>   "l" : param es en realidad el valor r
        !> \param[in] sigma(real*8) Valor la ????
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] x(real*8) Valor inicial en el que se evalua el modelo
        !> \param[in] npoints(integer) Tamaño de la simulación
        !> \param[out] path(real*8) Matriz de valores de la simulación
        !> \param[in][optional] linf(real*8) Limite superior
        !-------------------------------------------------------------------
        implicit none
        character, intent(in) :: type
        real*8, intent(in) :: param,sigma,delta,x
        integer, intent(in) :: npoints
        real*8, intent(out) :: path(npoints)
        real*8, intent(in) , optional :: linf
        real*8 y1,y2,inf
        integer :: i
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        !
        path(1)=x
        do i=2,npoints
            call DriftParam(type,param,path(i-1),y1,inf)
            call DiffusionParam(type,sigma,path(i-1),y2,inf)
            call MilsteinStep(delta,path(i-1),y1,y2,sigma,path(i))
        enddo
        ! print*,path
        return
    end subroutine
end module simulation
