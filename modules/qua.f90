module qua
    use integral
    implicit none
contains
    subroutine Qua_Var_L(npoints,path,delta,sigmahat)
        implicit none
        integer :: npoints,i
        real*8 :: sigmahat,path(npoints),delta,num,dem
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
        implicit none
        integer npoints,i
        real*8 sigmahat,path(npoints),delta,num,dem,linf,pathint(npoints)
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