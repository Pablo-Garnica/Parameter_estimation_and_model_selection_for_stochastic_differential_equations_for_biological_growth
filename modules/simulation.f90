module simulation
    implicit none
    contains
    !Buscar que todas las subrutinas sean tipo input output
    !INTENT(IN), INTENT(OUT) or INTENT(INOUT).
    !En funciones con result para el out
    function unif() result(y)
        implicit none
        real*8 :: y
        call random_number(y)
        return
    end function unif

    function boxmuller() result(y)
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
        real*8 x
        x=boxmuller()
        return
    end subroutine

    subroutine BrownianStep(delta,startx,endx)
        implicit none
        real*8 delta,startx,endx,var
        call normalvar(var) 
        endx=startx+sqrt(delta)*var
        return
    end subroutine

    subroutine MilsteinStep(delta,startx,drix,difx,sigma,endx)
        implicit none
        real*8 delta,startx,endx,drix,difx,der
        real*8 sigma,W,xx,brow
        xx=0.0
        call BrownianStep(delta,xx,brow)
        der=sigma
        W=brow
        endx=startx+drix*delta+difx*W+(1.0/2.0)*difx*der*(W**2-delta)
        return
    end subroutine
    !Gompertz y logistic
    subroutine DiffusionParameter_L_G(sigma,x,y)
        implicit none
        real*8 sigma,x,y
        y=sigma*x
        return
    end subroutine
    subroutine DiffusionParameter_V_B(sigma,linf,x,y)
        implicit none
        real*8 sigma,x,y,linf
        y=sigma*(linf-x)
        return
      end subroutine
    !!!!! Diffusion Gompertz, Logistic, Von Ber
    subroutine DiffusionParam(type,sigma,x,y,linf)
        implicit none    
        real*8 :: sigma,x,y,inf
        real*8 , optional :: linf
        character :: type
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
    !!!!!!Gompertz
    subroutine DriftParameter_G(beta,x,y)
        implicit none 
        real*8 beta,x,y
        y=-beta*x*log(x)
        return
    end subroutine
    !!!!!!Logistic
    subroutine DriftParameter_L(r,x,y)
        implicit none 
        real*8 r,x,y
        y=r*x*(1.0-x)
        return
      end subroutine
    !!!!!!!!Von Ber
    subroutine DriftParameter_V(kappa,linf,x,y)
        implicit none 
        real*8 kappa,linf,x,y
        y=kappa*(linf-x)
        return
    end subroutine
    !!!!!!!!DriftParam para los 3 modelos
    !!!!!!!
    subroutine DriftParam(type,beta,x,y,linf)
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
        !!!!!!!!recordemos que logaritmo es una funcion de los reales positivos
        !!!!!!!!Al hacer la simulacion podemos tener valores negativos
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
    !!!!!!!!!!!!!!!!!!!!
    subroutine SIM(type,beta,sigma,delta,x,npoints,path,linf)
        implicit none
        integer npoints,i
        character type
        real*8 delta,x,path(npoints),y1,y2
        real*8 beta,sigma,inf
        real*8 , optional :: linf
        ! ,linf
        !
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        !
        path(1)=x
        do i=2,npoints
            call DriftParam(type,beta,path(i-1),y1,inf)
            call DiffusionParam(type,sigma,path(i-1),y2,inf)
            call MilsteinStep(delta,path(i-1),y1,y2,sigma,path(i))
        enddo
        return
    end subroutine
end module simulation
