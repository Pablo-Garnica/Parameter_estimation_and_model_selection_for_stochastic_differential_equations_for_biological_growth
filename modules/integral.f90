module integral
    implicit none
contains
    subroutine ItoIntegrate(npoints,path,W,integral)
        implicit none
        integer npoints,i
        real*8 path(npoints),integral,W(npoints)
        integral=0.0
        do i=2,npoints
        integral=integral+(W(i)-W(i-1))*(path(i-1)+path(i))/2.0
        enddo
    return
    end subroutine
    !
    subroutine Integrate(npoints,path,delta,integral)
        implicit none
        integer npoints,i
        real*8 delta,path(npoints),integral
        integral=0.0
        do i=2,npoints
            integral=integral+delta*(path(i-1)+path(i))/2.0
        enddo
        return
    end subroutine
end module integral