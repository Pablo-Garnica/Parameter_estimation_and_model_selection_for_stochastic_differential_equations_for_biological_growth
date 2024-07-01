module aic
    use integral
    implicit none
contains
    subroutine AIC_GOM(delta,b,sigma,npoints,path,aic)
        implicit none
        integer npoints
        real*8 delta,b,sigma,path(npoints),I1,I2,aic,pathint1(npoints),pathint2(npoints)
        pathint1(:)=-(b*path(:)*LOG(path(:)))/((sigma**2)*(path(:)**2))
        pathint2(:)=((b**2)*(path(:)**2))/((sigma**2)*(path(:)**2))
        call ItoIntegrate(npoints,pathint1,path,I1)
        call Integrate(npoints,pathint2,delta,I2)
        aic=-2*(I1-I2/2)
        return
    end subroutine
  !____________________________________________________  
    subroutine AIC_VON(delta,kappa,sigma,linf,npoints,path,aic)
        implicit none
        integer npoints
        real*8 delta,kappa,linf,sigma,path(npoints),I1,I2,aic,pathint1(npoints),pathint2(npoints)
        pathint1(:)=kappa/((sigma**2)*(linf-path(:)))
        pathint2(:)=(kappa**2)/(sigma**2)
        call ItoIntegrate(npoints,pathint1,path,I1)
        call Integrate(npoints,pathint2,delta,I2)
        aic=-2*(I1-I2/2)
        return
    end subroutine
  !____________________________________________________  
    subroutine AIC_LOG(delta,r,sigma,npoints,path,aic)
        implicit none
        integer npoints
        real*8 delta,r,sigma,path(npoints),I1,I2,aic,pathint1(npoints),pathint2(npoints)
        pathint1(:)=-(r*path(:)*(1.0-path(:)))/((sigma**2)*(path(:)**2))
        pathint2(:)=((r**2)*((path(:)*(1.0-path(:)))**2))/((sigma**2)*(path(:)**2))
        call ItoIntegrate(npoints,pathint1,path,I1)
        call Integrate(npoints,pathint2,delta,I2)
        aic=-2*(I1-I2/2)
    return
    end subroutine
end module aic