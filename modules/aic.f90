module aic
    use integral
    implicit none
contains
    subroutine AIC_GOM(delta,b,sigma,npoints,path,aic)
        !-------------------------------------------------------------------
        !> \brief Calcula ???
        ! 
        !> \param[in] delta(real*8)
        !> \param[in] b(real*8)
        !> \param[in] sigma(real*8)
        !> \param[in] npoints(integer)
        !> \param[in] path(real*8)
        !> \param[out] aic(real*8)
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
    subroutine AIC_VON(delta,kappa,sigma,linf,npoints,path,aic)
        !-------------------------------------------------------------------
        !> \brief Calcula ???
        ! 
        !> \param[in] delta(real*8)
        !> \param[in] kappa(real*8)
        !> \param[in] sigma(real*8)
        !> \param[in] linf(real*8)
        !> \param[in] npoints(integer)
        !> \param[in] path(real*8)
        !> \param[out] aic(real*8)
        !-------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: delta,kappa,sigma,linf
        integer, intent(in) :: npoints
        real*8, intent(out) :: aic
        real*8 :: path(npoints),I1,I2,pathint1(npoints),pathint2(npoints)
        pathint1(:)=kappa/((sigma**2)*(linf-path(:)))
        pathint2(:)=(kappa**2)/(sigma**2)
        call ItoIntegrate(npoints,pathint1,path,I1)
        call Integrate(npoints,pathint2,delta,I2)
        aic=-2*(I1-I2/2)
        return
    end subroutine
    subroutine AIC_LOG(delta,r,sigma,npoints,path,aic)
        !-------------------------------------------------------------------
        !> \brief Calcula ???
        ! 
        !> \param[in] delta(real*8)
        !> \param[in] r(real*8)
        !> \param[in] sigma(real*8)
        !> \param[in] npoints(integer)
        !> \param[in] path(real*8)
        !> \param[out] aic(real*8)
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
end module aic