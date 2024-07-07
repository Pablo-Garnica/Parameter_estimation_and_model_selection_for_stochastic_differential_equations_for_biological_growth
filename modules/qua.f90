module qua
    use integral
    implicit none
contains
    subroutine Qua_Var_L(npoints,path,delta,sigmahat)
        !-------------------------------------------------------------------
        !> \brief Calcula ???
        ! 
        !> \param[in] npoints(integer)
        !> \param[in] path(real*8)
        !> \param[in] delta(real*8)
        !> \param[out] sigmahat(real*8)
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
        !> \brief Calcula ???
        ! 
        !> \param[in] npoints(integer)
        !> \param[in] path(real*8)
        !> \param[in] delta(real*8)
        !> \param[in] linf(real*8)
        !> \param[out] sigmahat(real*8)
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
end module qua