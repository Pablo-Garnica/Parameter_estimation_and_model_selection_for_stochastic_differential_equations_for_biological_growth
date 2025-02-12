module simulation
    use usefull_func
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

    subroutine SIM_ITER(type,param,sigma,delta,x,npoints,niter,path,linf)
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
        !> \param[in] npoints(integer) Numero de ejecuciones del modelo
        !> \param[out] path(real*8) Matriz de valores de la simulación
        !> \param[in][optional] linf(real*8) Limite superior
        !-------------------------------------------------------------------
        implicit none
        character, intent(in) :: type
        real*8, intent(in) :: param,sigma,delta,x
        integer, intent(in) :: npoints,niter
        real*8, intent(out) :: path(niter,npoints)
        real*8, intent(in) , optional :: linf
        real*8 :: path_for(npoints),inf
        integer :: i
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        do i=1,niter
            call SIM(type,param,sigma,delta,x,npoints,path_for,inf)
            path(i,:)=path_for
        enddo
    end subroutine
    ! subroutine SIM_CHOOSE(npoints,niter,nsteps,path_out,linf)
    subroutine SIM_CHOOSE(path,npoints,niter,nsteps,path_out,linf)
        !-------------------------------------------------------------------
        !> \brief Elige una trayectoria a partir de una serie de trayectorias
        ! 
        !> \param[in] path(real*8) Matriz de iteraciones de matrices valores de 
        !> la simulación
        !> \param[in] npoints(integer) Tamaño de la simulación
        !> \param[in] nsteps(integer) Numero de pasos realizados para el algoritmo
        !> \param[in] niter(integer) Numero de iteraciones
        !> \param[in] path_out(real*8) Matriz de trayectorias seleccionadas
        !> \param[in][optional] linf(real*8) Limite superior
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: npoints, nsteps, niter
        real*8, intent(in) :: path(niter,npoints)
        real*8, intent(out) :: path_out(nsteps)
        real*8, intent(in) , optional :: linf
        ! real*8 :: path(niter,npoints)
        integer i,size_jump,puntos(nsteps),ini(nsteps)
        real*8 a,inf
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        size_jump=(npoints-1)/(nsteps-1)   
        do i=1,nsteps
            puntos(i)=(size_jump)*(i-1)+1
        enddo
        call random_number (a)
        ini(1)=FLOOR(SIZE(path(:,1))*a)+1
        path_out(1)= path(ini(1),1)
        do i=2,nsteps
            call CHOOSE_OBS(path(:,puntos(i)),niter,path_out(i-1),inf,path_out(i))
        enddo        
        end subroutine
    subroutine CHOOSE_OBS(sample,ss,obs,Infi,nobs)
        !-------------------------------------------------------------------
        !> \brief ??????
        ! 
        !> \param[in] sample(real*8) ??
        !> \param[in] ss(integer) Tamaño de la simulación
        !> \param[in] obs(integer) Numero de pasos observaciones
        !> \param[in] Infi(integer) Limite superior
        !> \param[in] nobs(real*8) ????
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: sample(ss),obs,Infi
        integer, intent(in) :: ss
        real*8, intent(out) :: nobs
        integer ban,i
        real*8 Mean,Variance,StdDev,x(1),saort(ss)
        call MeanVar(sample,ss,Mean,Variance,StdDev)
        call sort(sample,ss,saort)
        call trun_normal(1,obs,Variance,obs-Variance,Infi,x)
        ban=0
        i=1
        do while(ban<1)
            if(saort(i)>x(1))then
                nobs=saort(i)
                ban=1
            else if(i==ss)then
                nobs=saort(i)
                ban=1
            else
                ban=0
                i=i+1
            endif
        enddo 
    end subroutine choose_obs
    subroutine MeanVar(x,n, Mean, Variance, StdDev)
        !-------------------------------------------------------------------
        !> \brief Calcula estadisticos mean, variance y desviación estandar
        ! 
        !> \param[in] x(real*8) Vector de numeros reales
        !> \param[in] n(integer) Longitud de vector
        !> \param[out] Mean(integer) Promedio
        !> \param[out] Variance(integer) Varianza
        !> \param[out] StdDev(real*8) Desviación estandar
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: x(n)
        real*8, intent(out) :: Mean, Variance, StdDev
        real*8 :: Suma,SumSQR
        Suma=SUM(x)
        SumSQR = SUM(x*x)
        Mean = Suma/ n
        Variance = (SumSQR - Suma*Suma/n)/(n-1)
        StdDev   = SQRT(Variance)
    end subroutine MeanVar
    subroutine sort(x,n,xs)
        !-------------------------------------------------------------------
        !> \brief Ordena el vector input
        ! 
        !> \param[in] x(real*8) Vector de numeros reales
        !> \param[in] n(integer) Longitud de vector
        !> \param[out] xs(integer) Vector ordenado
        !-------------------------------------------------------------------
        integer, intent(in) :: n
        real*8, intent(in) :: x(n) 
        real*8, intent(out) :: xs(n) 
        real*8 :: temp 
        integer :: i,j 
        xs=x 
        do 1 i=1,n-1 
          do 2 j=i+1,n 
            if (xs(i).gt.xs(j)) then 
              temp = xs(i) 
              xs(i) = xs(j) 
              xs(j) = temp 
            end if 
          2 continue 
        1 continue 
    end subroutine sort
    subroutine trun_normal(n,mean,var,inf,sup,x)
        !-------------------------------------------------------------------
        !> \brief Función normal truncada
        ! 
        !> \param[in] n(integer) Longitud de vector
        !> \param[in] mean(real*8) Promedio
        !> \param[in] var(real*8) Varianza
        !> \param[in] inf(real*8) Limite superior
        !> \param[in] sup(real*8) ???
        !> \param[out] x(real*8) Normal truncada
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: mean,var,inf,sup
        real*8, intent(out) :: x(n)
        real*8 alpha,beta,cdinf,cdsup,xt,pp
        real*8 :: sd,u
        integer :: i
        sd=sqrt(var)
        alpha=(inf-mean)/sd
        beta=(sup -mean)/sd
        do i=1,n
            call random_number(u)
            call cdnormal(alpha,cdinf)
            call cdnormal(beta,cdsup)
            pp=cdinf+u*(cdsup-cdinf)
            if(pp==1.00) then
                pp=0.99999995
            endif  
            call normal_01_cdf_inv(pp,xt)
            x(i)=sd*xt+mean
        enddo
        return
    end subroutine trun_normal
    subroutine cdnormal(x,prob)
        !-------------------------------------------------------------------
        !> \brief ???????
        ! 
        !> \param[in] x(real*8) Valor al que se le aplicará la función
        !> \param[out] prob(real*8) Probabilidad
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: x
        real*8, intent(out) :: prob
        real*8 :: xe,err
        xe=x/sqrt(2.0D0)
        call error(xe,err)
        prob= 1.0D0/2.0D0*(1.0D0+err)
    return
    end  
    subroutine error(x, err)
        !-------------------------------------------------------------------
        !> \brief Funcion de error
        ! 
        !> \param[in] x(real*8) Valor al que se le aplicará la función
        !> \param[out] err(real*8) Error
        !-------------------------------------------------------------------
        IMPLICIT double precision (A-H,O-Z)
        double precision :: EPS, PI, X2, ER, R, C0
        integer :: K
        ! Asignación de constantes
        EPS = 1.0D-15
        PI = 3.141592653589793D0
        ! Calcular X^2
        X2 = X * X
        ! Primera condición para X < 3.5
        IF (DABS(X) .LT. 3.5D0) THEN
            ER = 1.0D0
            R = 1.0D0
            DO 10 K = 1, 50
                R = R * X2 / (K + 0.5D0)
                ER = ER + R
                IF (DABS(R) .LE. DABS(ER) * EPS) GO TO 15
            10 CONTINUE
        15  C0 = 2.0D0 / DSQRT(PI) * X * DEXP(-X2)
            ERR = C0 * ER
        ! Segunda condición para X >= 3.5
        ELSE
            ER = 1.0D0
            R = 1.0D0
            DO 20 K = 1, 12
                R = -R * (K - 0.5D0) / X2
                ER = ER + R
            20 CONTINUE
            C0 = DEXP(-X2) / (DABS(X) * DSQRT(PI))
            ERR = 1.0D0 - C0 * ER
        END IF
        ! Ajustar el signo si X < 0
        IF (X .LT. 0.0D0) THEN
            ERR = -ERR
        END IF
        RETURN
    end subroutine error
end module simulation
