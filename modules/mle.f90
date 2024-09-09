module mle
  use integral
  use simulation
    implicit none
contains
      subroutine MLE_b(npoints,path,delta,bhat)
        !-------------------------------------------------------------------
        !> \brief Calcula de MLE (Maximum Likelihood Estimator) para el
        !> modelo Gompertz
        ! 
        !> \param[in] path(real*8) Observaciones de la SDE (stochastic 
        !> differential equation)
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] npoints(integer) Numero de observaciones
        !> \param[out] bhat(real*8) Estimado b
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: npoints
        real*8, intent(in) :: delta
        real*8, intent(out) :: bhat
        real*8 :: j,dem,num,path(npoints),lpath(npoints)
        lpath(:)=LOG(path(:))
        num=((npoints-1)*SUM(lpath(1:(npoints-1))*lpath(2:npoints))-SUM(lpath(2:npoints))*SUM(lpath(1:(npoints-1))))
        dem=((npoints-1)*SUM(lpath(1:(npoints-1))**2)-SUM(lpath(1:(npoints-1)))**2)
        j=num/dem
        bhat=-LOG(j)/delta
        return
    end subroutine
    subroutine MLE_r(npoints,path,delta,rhat)
        !-------------------------------------------------------------------
        !> \brief Calcula de MLE (Maximum Likelihood Estimator) para el
        !> modelo Logistico
        ! 
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] npoints(integer) Numero de observaciones
        !> \param[in] path(real*8) Observaciones de la SDE (stochastic 
        !> differential equation)
        !> \param[out] rhat(real*8) Estimador r
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: npoints
        real*8, intent(in) :: delta
        real*8, intent(out) :: rhat
        real*8 :: path(npoints),I1,I2
        call ItoIntegrate(npoints,(1.0-path(:))/path(:),path,I1)
        call Integrate(npoints,(1.0-path(:))**2,delta,I2)
        rhat=I1/I2
        return
    end subroutine
    subroutine MLE_kappa(npoints,path,delta,linf,kappahat)
        !-------------------------------------------------------------------
        !> \brief Calcula de MLE (Maximum Likelihood Estimator) para el
        !> modelo Von Bert
        ! 
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] linf(real*8) Limite superior
        !> \param[in] npoints(integer) Numero de observaciones
        !> \param[in] path(real*8) Observaciones de la SDE (stochastic 
        !> differential equation)
        !> \param[out] kappahat(real*8) Estimador kappa
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: delta, linf
        integer, intent(in) :: npoints
        real*8, intent(out) :: kappahat
        real*8 :: path(npoints),I1,pathint(npoints),time
        time=delta*(npoints-1)
        pathint(:)=1.0/(linf-path(:))
        call ItoIntegrate(npoints,pathint,path,I1)
        kappahat=I1/time
        return
    end subroutine
    subroutine MLE_times(path,delta,linf,size,npoints,kappahat,sigmahat)
        !-------------------------------------------------------------------
        !> \brief Calcula de MLE (Maximum Likelihood Estimator) de ???
        ! 
        !> \param[in] path(real*8)
        !> \param[in] delta(real*8)
        !> \param[in] linf(real*8)
        !> \param[in] size(integer)
        !> \param[in] npoints(integer)
        !> \param[in] pathlam(real*8)
        !> \param[out] kappahat(real*8)
        !> \param[out] sigmahat(real*8)
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: delta,linf
        integer, intent(in) :: npoints,size
        real*8, intent(out) :: sigmahat(size),kappahat(size)
        real*8 path(npoints),geo_brown(npoints)
        integer i,size_jump,puntos(size),n
        size_jump=(npoints-1)/(size-1) 
        do i=1,size
          puntos(i)=(size_jump)*(i)+1
        enddo
        call trans_geo_brown(linf,npoints,path,geo_brown)
        do i=1,size
          n=puntos(i)
          call MLE_GB(geo_brown(1:n),delta,n,sigmahat(i),kappahat(i))
        enddo
        return
      end subroutine
    subroutine trans_geo_brown(linf,npoints,path,geo_brown)
        !-------------------------------------------------------------------
        !> \brief ???
        ! 
        !> \param[in] linf(real*8)
        !> \param[in] npoints(integer)
        !> \param[in] path(real*8)
        !> \param[out] geo_brown(real*8)
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: npoints
        real*8, intent(in) :: linf
        real*8, intent(out) :: geo_brown(npoints)
        real*8 :: path(npoints)
        geo_brown(:)=linf-path(:)       
    return
    end subroutine
    subroutine MLE_GB(path,delta,npoints,sigmahat,kappahat)
        !-------------------------------------------------------------------
        !> \brief Calcula de MLE (Maximum Likelihood Estimator) de ???
        ! 
        !> \param[in] path(real*8)
        !> \param[in] delta(real*8)
        !> \param[in] npoints(integer)
        !> \param[out] sigmahat(real*8)
        !> \param[out] kappahat(real*8)
        !-------------------------------------------------------------------
          implicit none
          real*8, intent(in) :: delta
          integer, intent(in) :: npoints
          real*8, intent(out) :: sigmahat,kappahat
          real*8 :: path(npoints),difpath(npoints-1),a,b,c,d,e
          difpath(:)=path(2:npoints)-path(1:(npoints-1))
          a=SUM(LOG(path(2:npoints)))
          b=SUM(LOG(path(1:(npoints-1))))
          c=SUM(LOG(path(2:npoints))**2)
          d=SUM(LOG(path(2:npoints))*LOG(path(1:(npoints-1))))
          e=SUM(LOG(path(1:(npoints-1)))**2)
          sigmahat=sqrt((-a**2+2*a*b-b**2+c*npoints-2*d*npoints+e*npoints)/((npoints**2)*delta))
          kappahat=-((a-b)/(delta*(npoints-1))+(sigmahat**2)/2)
          return
    end subroutine
    subroutine MLE_(type_model,npoints,path,delta,paramhat,linf)
      !-------------------------------------------------------------------
      !> \brief Calcula MLE (Maximum Likelihood Estimator) para
      !> los modelos
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
      !> \param[out] paramhat(real*8) Valor de ??? si type_model es:
      !>   "v" : paramhat es en realidad el valor kappahat
      !>   "g" : paramhat es bhat
      !>   "l" : paramhat es en realidad el valor rhat
      !> \param[in][optional] linf(real*8) Limite superior
      !-------------------------------------------------------------------
      implicit none
      character, intent(in) :: type_model
      integer, intent(in) :: npoints
      real*8 :: path(npoints)
      real*8, intent(in) :: delta
      real*8, intent(out) :: paramhat
      real*8, intent(in) , optional :: linf
      real*8 :: inf
      if (present(linf)) then
          inf = linf
      else
          inf = 999999999999999999999999999999.00
      end if
      if (type_model.eq."g") then
          call MLE_b(npoints,path,delta,paramhat)
      else if (type_model.eq."l") then
          call MLE_r(npoints,path,delta,paramhat)
      else if (type_model.eq."v") then
          call MLE_kappa(npoints,path,delta,inf,paramhat)
      else
          print *, "Error"
      end if
      return
  end subroutine
end module mle