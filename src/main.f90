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
    character :: type_model
    real*8 :: model_param,sigma,delta,xstart,linf
    integer :: npoints
    !Output
    character*20 :: name_model,model_name_param
    character*100 :: param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom
    real*8 :: paramhat, sigmahat, aic_param
    real*8 , allocatable :: path(:)
    !Program
    call menu_random_seed()
    call menu_model_type(type_model)
    call dict_models(type_model,name_model,model_name_param, &
    param_recom,sigma_recom,delta_recom,xstart_recom,npoints_recom)
    call menu_input_param(type_model,model_param,sigma,delta,xstart,npoints,linf)
    allocate(path(npoints))
    call menu_result(type_model,model_param,sigma,delta,xstart,npoints,linf,paramhat,sigmahat,aic_param)
    deallocate(path)
end program main