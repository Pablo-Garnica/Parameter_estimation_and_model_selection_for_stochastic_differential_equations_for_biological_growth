program main
    use rand_seed
    use menu
    use simulation
    use mle
    use qua
    use aic
    use rand_seed
    implicit none
    !Input
    character :: type_model,type_execution
    real*8 :: model_param,sigma,delta,xstart,linf
    integer :: npoints,nobs,niter,nsteps,nmc,status
    !Output
    character*20 :: name_model,model_name_param
    character*100 :: param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom
    real*8 :: paramhat, sigmahat, aic_param
    real*8 , allocatable :: path(:)
    !Program
    call menu_random_seed()
    call menu_model_type(type_model)
    call menu_type_execution(type_execution)

    call dict_models(type_model,name_model,model_name_param, &
    param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom)
    ! call menu_input_param(type_model,type_execution,model_param,sigma,delta,xstart,npoints,&
    ! nobs,niter,nsteps,linf)
    call  menu_input_param(type_model,type_execution,model_param,sigma,delta,xstart,npoints,&
    nobs,niter,nsteps,nmc,linf)
    ! call menu_input_param(type_model,type_execution,model_param,sigma,delta,xstart,npoints,nobs,linf)
    ! call menu_input_param(type_model,model_param,sigma,delta,xstart,npoints,type_model,linf)
    allocate(path(npoints))

    call menu_result(type_model,type_execution,model_param,sigma,delta,xstart,npoints,nobs,&
    niter,nsteps,nmc,linf,paramhat,sigmahat,aic_param,path)
    ! call menu_result(type_model,type_execution,model_param,sigma,delta,xstart,npoints,linf,paramhat,sigmahat,aic_param,path)
    
    call menu_export(type_model,path,npoints)

    call menu_export_params(type_model,path,npoints,model_param,sigma,delta,linf)
    PRINT *, "Ejecutando Python"
    status = SYSTEM("plot_sim")
    IF (status /= 0) THEN
        PRINT *, "Error al ejecutar el script de Python"
    ELSE
        PRINT *, "Script de Python ejecutado con éxito"
    END IF
    deallocate(path)
end program main