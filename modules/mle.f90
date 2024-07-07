module mle
  use integral
  use simulation
    implicit none
contains
    subroutine MLE_G(path,delta,npoints,sigmahat,bhat)
        !-------------------------------------------------------------------
        !> \brief Calcula de MLE (Maximum Likelihood Estimator) para el
        !> modelo Gompertz
        ! 
        !> \param[in] path(real*8) Observaciones de la SDE (stochastic 
        !> differential equation)
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] npoints(integer) Numero de observaciones
        !> \param[out] sigmahat(real*8) Estimador sigma
        !> \param[out] bhat(real*8) Estimado b
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: npoints
        real*8, intent(in) :: delta
        real*8, intent(out) :: sigmahat,bhat
        real*8 :: j,k,dem,num,h,path(npoints)
        path(:)=LOG(path(:))  
        num=((npoints-1)*SUM(path(1:(npoints-1))*path(2:npoints))-SUM(path(2:npoints))*SUM(path(1:(npoints-1))))
        dem=((npoints-1)*SUM(path(1:(npoints-1))**2)-SUM(path(1:(npoints-1)))**2)
        j=num/dem
        k=(j*SUM(path(1:(npoints-1)))-SUM(path(2:npoints)))/(npoints-1)
        h=SUM((path(2:npoints)-j*path(1:npoints-1)+k)**2)/(npoints-1)
        bhat=-LOG(j)/delta
        sigmahat=SQRT((h*2*bhat)/(1-EXP(-2*bhat*delta)))      
        return
    end subroutine
    subroutine MLE_r(delta,npoints,path,rhat)
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
    subroutine MLE_kappa(delta,linf,npoints,pathlam,kappahat)
        !-------------------------------------------------------------------
        !> \brief Calcula de MLE (Maximum Likelihood Estimator) para el
        !> modelo Von Bert
        ! 
        !> \param[in] delta(real*8) Incremento del proceso de Wiener
        !> \param[in] linf(real*8) Limite superior
        !> \param[in] npoints(integer) Numero de observaciones
        !> \param[in] pathlam(real*8) Observaciones de la SDE (stochastic 
        !> differential equation)
        !> \param[out] kappahat(real*8) Estimador kappa
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: delta, linf
        integer, intent(in) :: npoints
        real*8, intent(out) :: kappahat
        real*8 :: pathlam(npoints),I1,pathint(npoints),time
        time=delta*(npoints-1)
        pathint(:)=1.0/(linf-pathlam(:))
        call ItoIntegrate(npoints,pathint,pathlam,I1)
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
end module mle