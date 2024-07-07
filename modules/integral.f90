module integral
    implicit none
contains
    subroutine ItoIntegrate(npoints,path,W,integral)
        !-------------------------------------------------------------------
        !> \brief Calcula la integral de Ito 
        !>
        !> \param[in] npoints(integer) Tamaño valores que en la cual se 
        !> calculará la integral
        !> \param[in] path(real*8) Matriz de valores en la cual se calculará
        !> la integral (dominio de la funcion)
        !> \param[in] W(real*8) Matriz de valores en la cual se calculará
        !> la integral (imagen de la funcion)
        !> \param[out] integral(real*8) Resultado de la integral de Ito
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in)  :: npoints
        real*8, intent(in) :: path(npoints),W(npoints)
        real*8, intent(out) :: integral
        integer :: i
        integral=0.0
        do i=2,npoints
            integral=integral+(W(i)-W(i-1))*(path(i-1)+path(i))/2.0
        enddo
    return
    end subroutine
    subroutine Integrate(npoints,path,delta,integral)
        !-------------------------------------------------------------------
        !> \brief Calcula la integral de Riemman 
        !>
        !> \param[in] npoints(integer) Tamaño valores que en la cual se 
        !> calculará la integral
        !> \param[in] path(real*8) Matriz de valores en la cual se calculará
        !> la integral (imagen de la funcion)
        !> \param[in] delta(real*8) Distancia entre las particiones (se 
        !> construyó una parición uniforme)
        !> \param[out] integral(real*8) Resultado de la integral de Riemman
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in)  :: npoints
        real*8, intent(in) :: path(npoints),delta
        real*8, intent(out) :: integral
        integer :: i
        integral=0.0
        do i=2,npoints
            integral=integral+delta*(path(i-1)+path(i))/2.0
        enddo
        return
    end subroutine
end module integral