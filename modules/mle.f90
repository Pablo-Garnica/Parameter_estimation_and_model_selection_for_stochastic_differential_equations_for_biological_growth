module mle
  use integral
  use simulation
    implicit none
contains
    !Gompertz
    ! subroutine OU(npoints,path,pOU)
    !     implicit none
    !     integer npoints
    !     real*8 path(npoints),pOU(npoints)
    !     pOU(:)=LOG(path(:))       
    !     return
    ! end subroutine
    subroutine MLE_G(path,delta,npoints,sigmahat,bhat)
        implicit none
        integer npoints
        real*8 delta,sigmahat,bhat
        real*8 path(npoints),j,k,dem,num,h
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
    !Logistic
    subroutine MLE_r(delta,npoints,path,rhat)
        implicit none
        integer npoints
        real*8 delta,path(npoints),rhat,I1,I2
        call ItoIntegrate(npoints,(1.0-path(:))/path(:),path,I1)
        call Integrate(npoints,(1.0-path(:))**2,delta,I2)
        rhat=I1/I2
        return
    end subroutine
    !Von Beb
    subroutine MLE_kappa(delta,linf,npoints,pathlam,kappahat)
        implicit none
        integer npoints
        real*8 delta,pathlam(npoints),kappahat,I1,linf,pathint(npoints),time
        time=delta*(npoints-1)
        pathint(:)=1.0/(linf-pathlam(:))
        call ItoIntegrate(npoints,pathint,pathlam,I1)
        kappahat=I1/time
        return
    end subroutine
    subroutine MLE_times(path,delta,linf,size,npoints,kappahat,sigmahat)
        implicit none
        integer npoints,i,size,size_jump,puntos(size),n
        real*8 delta,sigmahat(size),kappahat(size),linf
        real*8 path(npoints),geo_brown(npoints)
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
        implicit none
        integer npoints
        real*8 path(npoints),linf,geo_brown(npoints)
        geo_brown(:)=linf-path(:)       
    return
    end subroutine
    !_____________________________________________________________________
    !Dependencia de MLE_times
    subroutine MLE_GB(path,delta,npoints,sigmahat,kappahat)
      implicit none
      integer npoints
      real*8 delta,sigmahat,kappahat
      real*8 path(npoints),difpath(npoints-1),a,b,c,d,e
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