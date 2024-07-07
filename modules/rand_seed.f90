module rand_seed
    implicit none
contains
    subroutine random_seed_set_up(seed)
        !-------------------------------------------------------------------
        !> \brief Fija la semilla aleatorea
        !>
        !> \param[in] seed(integer) Valor en el que se fija la semilla 
        !> aleatoria
        !-------------------------------------------------------------------
        implicit none
        integer, intent(in) :: seed
        integer :: n
        integer, allocatable :: seed_array(:)
        ! Obtener el tamaño necesario para la semilla
        call random_seed(size=n)
        ! Asignar espacio para la semilla
        allocate(seed_array(n))
        ! Llenar el arreglo de la semilla
        seed_array = seed
        ! Inicializar la semilla
        call random_seed(put=seed_array)
        deallocate(seed_array)
        return
    end subroutine
end module rand_seed