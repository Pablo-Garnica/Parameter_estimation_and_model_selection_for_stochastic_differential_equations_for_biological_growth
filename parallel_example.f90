! sh gfortran -fopenmp -o parallel_example parallel_example.f90
! Directiva de Compilación OpenMP:

! !$omp parallel do num_threads(num_threads) private(i) shared(a, b, c, n): Esta directiva indica al 
! compilador que el bucle do debe ser ejecutado en paralelo utilizando num_threads hilos. 
! La cláusula private(i) especifica que cada hilo tiene su propia copia de la variable i, 
! mientras que shared(a, b, c, n) indica que las matrices a, b, c y la variable n son compartidas 
! entre todos los hilos.
! !$omp end parallel do: Marca el final de la región paralela.
program parallel_example
    use omp_lib
    implicit none
    integer :: i, n
    integer, parameter :: num_threads = 4
    real :: a(100), b(100), c(100)

    n = 100
    ! Inicializar los vectores
    do i = 1, n
        a(i) = i * 1.0
        b(i) = (n - i) * 1.0
    end do

    ! Paralelizar el bucle usando OpenMP
    !$omp parallel do num_threads(num_threads) private(i) shared(a, b, c, n)
    do i = 1, n
        c(i) = a(i) + b(i)
    end do
    !$omp end parallel do

    ! Imprimir algunos resultados
    print *, "Primeros 10 resultados del vector c:"
    do i = 1, 10
        print *, c(i)
    end do
end program parallel_example