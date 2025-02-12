!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dudas:
! - En que parte del programa de debe de agregar esto?
! - Significado inputs
! - En la subrutina *diffusionBridge* no hay posibilidad 
! de ir a 2
! - En la subrutina *diffusionBridge* el parametro *y* no hace nada 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! program name
!     implicit none
!     print*,'hola'
! end program name
module em
    use usefull_func
    use simulation
    use mle
    use qua
    implicit none
contains
    subroutine choose_data(delta,npoints,path_in,nobs,delta_bridge,path_out)
        !-------------------------------------------------------------------
        !> \brief Hace una muestra de uniforme de tamaño nobs de un path de
        !>  tamaño npoints
        ! 
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] npoints(integer) Tamaño de la simulación
        !> \param[in] path_in(real*8) Matriz la cual se va a hacer la extracción
        !> \param[in] nobs(integer) De muestra 
        !> \param[out] delta_bridge(real*8) Incremento del proceso de Wiener 
        !> para la muestra
        !> \param[out] path_in(real*8) Matriz al realizar la muestra
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: npoints,nobs
        real*8, intent(in) :: delta,path_in(npoints)
        real*8, intent(out) :: delta_bridge,path_out(nobs)
        integer :: i,size_jump,label(nobs)
        size_jump=(npoints-1)/(nobs-1)   
        delta_bridge=delta*size_jump
        do i=1,nobs
            label(i)=(size_jump)*(i-1)+1
            path_out(i)=path_in(label(i))
        enddo
    return
    end subroutine
    subroutine EM_MC(type_model,beta0,sigma0,path,nobs,delta_bridge,niter,nmc,nsteps,sigma_vec,beta_vec,linf)
        !-------------------------------------------------------------------
        !> \brief Hace una muestra de uniforme de tamaño nobs de un path de
        !>  tamaño npoints
        ! 
        !> \param[in] type_model(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[in] beta0(real*8) Beta inicial
        !> \param[in] sigma0(real*8) Sigma inicial
        !> \param[in] path(real*8) Path de simulación
        !> \param[in] nobs(integer) Tamaño de path
        !> \param[in] delta_bridge(real*8) Incremento del proceso de Wiener 
        !> para el path
        !> \param[in] niter(integer) Numero de iteraciones
        !> \param[in] nmc(integer) Tamaño de vector de betas y sigmas
        !> \param[in] nsteps(integer) Numero de pasos que realiza el
        !> algoritmo CompleteBridge
        !> \param[in] linf[optional](real*8) Limite superior
        !> \param[out] beta_vec(real*8) Vector de betas
        !> \param[out] sigma_vec(real*8) Vector de sigmas
        !-------------------------------------------------------------------
        implicit none
        character, intent(in) :: type_model
        integer, intent(in) :: nobs,niter,nmc,nsteps
        real*8, intent(in) :: beta0,sigma0,path(nobs),delta_bridge
        real*8, intent(in), optional  :: linf
        real*8, intent(out) :: sigma_vec(niter),beta_vec(niter)
        integer :: i,j,npoints,numrej
        real*8 :: delta,complete_bridge((nsteps+1)*(nobs-1)+1), path_(nobs),bhs(nmc),shs(nmc)
        real*8 :: inf
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        !
        if (type_model.eq."g") then
            path_(:)=LOG(path(:))    
        else
            path_ = path
        end if
        ! 
        beta_vec(1)=beta0
        sigma_vec(1)=sigma0
        npoints=(nsteps+1)*(nobs-1)+1
        delta=delta_bridge/(nsteps+1)
        
        ! 
        do i=2,niter 
            print*,"iter_EM",i,beta_vec(i-1),sigma_vec(i-1)
            do j=1,nmc
                ! print*,j
                call CompleteBridge(type_model,beta_vec(i-1),sigma_vec(i-1),nsteps,delta_bridge, &
                    nobs,path_,complete_bridge,numrej,inf)

                call MLE_(type_model,npoints,complete_bridge,delta,bhs(j),inf)
                call Qua_Var(type_model,npoints,complete_bridge,delta,shs(j),inf)

                ! print*,j,"mle",bhs(j)
                ! print*,j,"path",complete_bridge
                print*,"-------------------------------------------------------------------"
                print*,"-------------------------------------------------------------------"
                print*,j,"path",complete_bridge
                print*,"-------------------------------------------------------------------"
                print*,"-------------------------------------------------------------------"
                print*,j,"path",path_
                print*,"-------------------------------------------------------------------"
                print*,j,"mle",bhs
                print*,"-------------------------------------------------------------------"
                print*,"-------------------------------------------------------------------"
            enddo
            print*,"-------------------------------------------------------------------"
            print*,i
            sigma_vec(i)=SUM(shs(:))/nmc
            beta_vec(i)=SUM(bhs(:))/nmc
            ! print*,"beta",beta_vec
            ! print*,"sigma_vec",sigma_vec
            print*,"-------------------------------------------------------------------"
            print*,"-------------------------------------------------------------------"
            print*,"beta:",beta_vec(i)
            print*,"-------------------------------------------------------------------"
            print*,"-------------------------------------------------------------------"
        enddo
        return
    end subroutine
    subroutine CompleteBridge(type_model,beta,sigma,nsteps,delta_bridge,nobs,path_,complete_bridge,numrej,linf)
        !-------------------------------------------------------------------
        !> \brief **********************************
        ! 
        !> \param[in] type_model(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[in] beta(real*8) Beta
        !> \param[in] sigma(real*8) Sigma
        !> \param[in] nsteps(integer) Numero de pasos
        !> \param[in] delta_bridge(real*8) Incremento del proceso de Wiener 
        !> para el path
        !> \param[in] nobs(integer) Tamaño de path
        !> \param[in] path(real*8) Path de simulación
        !> \param[in] linf[optional](real*8) Limite superior
        !> \param[out] complete_bridge(real*8) Path con algoritmo que acompleta
        !> \param[out] numrej(integer) ?????
        !-------------------------------------------------------------------
        implicit none
        character, intent(in) :: type_model
        real*8, intent(in) :: beta,sigma,delta_bridge,path_(nobs)
        integer, intent(in) :: nsteps,nobs
        integer, intent(out) :: numrej
        real*8, intent(out) :: complete_bridge(nsteps*(nobs-1)+1)
        real*8, intent(in), optional  :: linf
        real*8 :: inf
        integer i,ini,fin
        real*8 bridge(nsteps+2)
        
        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        do i=1,(nobs-1)
            ! print*,'CompleteBridge',i
            print*,'x',path_(i)
            print*,'y',path_(i+1)
            call diffusionBridge(type_model,beta,sigma,delta_bridge,nsteps,path_(i),path_(i+1),bridge,numrej,inf)
            ini=(nsteps+1)*(i-1)+1
            fin=(nsteps+1)*i+1
            complete_bridge(ini:fin)=bridge(:)
        enddo
        return
    end subroutine
    subroutine diffusionBridge(type_model,beta,sigma,delta,nsteps,x,y,bridge,numrej,linf)
        !-------------------------------------------------------------------
        !> \brief **********************************
        ! 
        !> \param[in] type_model(character) Debe de estar en los siguientes 
        !> valores:
        !>   "v" : Para el modelo Von Bert
        !>   "g" : Para el modelo Gompertz
        !>   "l" : Para el modelo Logistic
        !> \param[in] beta(real*8) Beta
        !> \param[in] sigma(real*8) Sigma
        !> \param[in] delta(real*8) Incremento del proceso de Wiener 
        !> para el path
        !> \param[in] nsteps(integer) Numero de pasos
        !> \param[in] x(real*8) ?????
        !> \param[in] y(real*8) ?????
        !> \param[in] linf[optional](real*8) Limite superior
        !> \param[out] bridge(real*8) ????
        !> \param[out] numrej(integer) ?????
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: beta,sigma,delta,x,y
        character, intent(in) :: type_model
        integer, intent(in) :: nsteps
        integer, intent(out) :: numrej
        real*8, intent(out) :: bridge(nsteps+2)
        real*8, intent(in), optional  :: linf
        integer i,j,mp
        real*8 delta_bridge,inf, points(nsteps+2),points1(nsteps+2),points2(nsteps+2),points3(nsteps+2)
        logical :: crossing_found

        if (present(linf)) then
            inf = linf
        else
            inf = 999999999999999999999999999999.00
        end if
        ! real*8 beta,sigma,delta,x,y
        ! character type_model
        ! integer nsteps
        ! integer numrej
        ! real*8 bridge(nsteps+2)
        ! real*8  linf
        ! integer i,j,mp
        ! real*8 delta_bridge,inf, points(nsteps+2),points1(nsteps+2),points2(nsteps+2),points3(nsteps+2)
        ! ! 

        inf = 999999999999999999999999999999.00
        ! 
        delta_bridge=delta/(nsteps+1)
        numrej=0
        !Nuevo
        ! 1 call SIM(type_model,beta,sigma,delta,x,(nsteps+2),points1,inf)
        !     call SIM(type_model,-beta,sigma,delta,y,(nsteps+2),points2,inf)
        ! call SIM(type_model,beta,sigma,delta,x,(nsteps+2),points1,inf)
        ! beta =0.6
        ! sigma = 0.1
        ! delta = 0.1
        ! nsteps = 3
        call SIM(type_model,beta,sigma,delta,x,(nsteps+2),points1,inf)
        call SIM(type_model,-beta,sigma,delta,y,(nsteps+2),points2,inf)
        print*,'param-------------------------------------------------------------------'
        print*,'bridge-------------------------------------------------------------------'
        print*,'pasos positivos',minval(points1),maxval(points1),points1
        print*,'pasos negativos',minval(points2),maxval(points2),points2
    
            ! Inicializamos crossing_found en falso
        crossing_found = .false.

        ! Buscamos el índice donde ocurre el cruce
        do i = 2, nsteps
            if (maxval(points1(1:nsteps)) > points2(i)) then
                points3(1:i) = points1(1:i) 
                points3((i+1):nsteps) = points2((i+1):nsteps)
                crossing_found = .true.
                exit
            endif
            if (.not.crossing_found) then
                print*,"Error no se realizo el cruce------------------------------"
            endif
        end do

        ! !Anterior
        ! ! 1 call SIM_OU(beta,sigma,delta_bridge,x,(nsteps+2),points1)
        ! ! 2 call SIM_OU(-beta,sigma,delta_bridge,y,(nsteps+2),points2)  
        !     do i=1,nsteps+2
        !         ! print*,'diffusionBridge paso "1",do',i
        !         points3(i)=points2(nsteps+3-i)
        !     enddo
        !     if (points3(1).lt.points1(1)) then
        !         do i=1,nsteps+2
        !         ! print*,'diffusionBridge paso "2",do',i
        !         if (points3(i).gt.points1(i)) then
        !                 mp=i
        !                 do j=mp,nsteps+2
        !                     ! print*,'diffusionBridge paso "3",do',j
        !                     ! print*,j
        !                     points(j)=points3(j)
        !                 enddo
        !                 ! goto 20
        !             endif
        !         enddo
        !     endif
        !     if (points3(1).gt.points1(1)) then
        !         do i=1,nsteps+2
        !             ! print*,'diffusionBridge paso "4",do',i
        !             if (points3(i).lt.points1(i)) then
        !                 mp=i
        !                 do j=mp,nsteps+2
        !                     ! print*,'diffusionBridge paso "5",do',j
        !                     points(j)=points3(j)
        !                 enddo
        !                 ! goto 20
        !             endif
        !         enddo
        !     endif
        ! ! numrej=numrej+1
        ! ! goto 1
        ! ! 
        ! ! 20 do i=1,(mp-1)
        ! do i=1,(mp-1)
        !         ! print*,'diffusionBridge paso "6",do',i
        !         points(i)=points1(i)
        !     enddo
        !     do i=1,nsteps+2
        !         ! print*,'diffusionBridge paso "7",do',i
        !         bridge(i)=points(i)
        ! enddo
        ! ! numrej=numrej+1
            
        return
    end subroutine
end module em