program main
    use rand_seed
    use simulation
    use mle
    use qua
    use aic
    implicit none
    integer :: seed

    real*8 :: unif_val, box_m, n_var
    !Browlian step
    real*8 :: delta,startx,endx
    !MilsteinStep
    real*8 :: drix,difx,sigma,beta
    !DiffusionParam
    real*8 :: x,y
    character :: type_m
    integer :: npoints, size
    !MLE
    real*8 :: sigmahat,bhat,rhat,kappahat,linf
    !AIC
    real*8 :: b,kappa,r
    real*8 :: aic_
    real*8,allocatable :: path(:),pathkappahat(:),pathsigmahat(:)
    !Semilla aleatorea
    seed = 1234
    call random_seed_set_up(seed)
    !
    size = 5
    npoints = 10
    linf = 999999999.0
    allocate(path(npoints),pathkappahat(size),pathsigmahat(size))

    !Validación unif
    unif_val = unif()
    print*,"Val funcion unif",unif_val
    !Valiacion boxmuller
    box_m = boxmuller()
    print*,"Val funcion boxmuller",box_m
    !Validacion normalvar
    call normalvar(n_var)
    print*,"Val subrutina normalvar",n_var
    !Validacion BrownianStep
    startx = 0.0
    delta = 0.5
    call BrownianStep(delta,startx,endx)
    print*,"Val subrutina BrownianStep",endx

    !Validacion MilsteinStep
    delta = 0.5
    startx = 0.0
    drix = 0.2
    difx= 0.6
    sigma= 1.5
    call MilsteinStep(delta,startx,drix,difx,sigma,endx)
    print*,"Val subrutina MilsteinStep",endx
    !Validacion DiffusionParam(type,sigma,x,y,linf)
    type_m = "l"
    x = 1.5
    call DiffusionParam(type_m,sigma,x,y)
    print*,"Val subrutina DiffusionParam",y
    !Validacion DriftParam(type,sigma,x,y,linf)
    type_m = "v"
    call DriftParam(type_m,sigma,x,y)
    print*,"Val subrutina DriftParam",y
    !Validacion DriftParam(type,sigma,x,y,linf)
    type_m = "v"
    !Validacion SIM(type,beta,sigma,delta,x,npoints,path,linf)
    beta = 0.5
    type_m = "g"
    call SIM(type_m,beta,sigma,delta,x,npoints,path)
    print*,"Val subrutina SIM",path
    !Validacion MLE_G(path,delta,npoints,sigmahat,bhat)
    call MLE_G(path,delta,npoints,sigmahat,bhat)
    print*,"Val subrutina MLE_G",bhat
    !Validacion MLE_r(delta,npoints,path,rhat)
    call MLE_r(delta,npoints,path,rhat)
    print*,"Val subrutina MLE_r",rhat
    !Validacion MLE_kappa(delta,linf,npoints,pathlam,kappahat)
    call MLE_kappa(delta,linf,npoints,path,kappahat)
    print*,"Val subrutina MLE_kappa",kappahat
    !Validacion MLE_GB(path,delta,npoints,sigmahat,kappahat)
    call MLE_GB(path,delta,npoints,sigmahat,kappahat)
    print*,"Val subrutina MLE_GB",sigmahat,kappahat   
    !Validacion MLE_times(path,delta,linf,size,npoints,kappahat,sigmahat)
    call MLE_times(path,delta,linf,size,npoints,pathkappahat,pathsigmahat)
    print*,"Val subrutina MLE_times",pathkappahat,pathsigmahat 
    !Validacion Qua_Var_L(npoints,path,delta,sigmahat)
    call Qua_Var_L(npoints,path,delta,sigmahat)
    print*,"Val subrutina Qua_Var_L",sigmahat
    !Validacion Qua_Var_VB(npoints,path,delta,linf,sigmahat)
    call Qua_Var_VB(npoints,path,delta,linf,sigmahat)
    print*,"Val subrutina Qua_Var_VB",sigmahat
    !Validacion AIC_GOM(delta,b,sigma,npoints,path,aic)
    b = 1.5
    call AIC_GOM(delta,b,sigma,npoints,path,aic_)
    print*,"Val subrutina AIC_GOM",aic_
    !Validacion AIC_VON(delta,kappa,sigma,linf,npoints,path,aic)
    kappa = 2.3
    call AIC_VON(delta,kappa,sigma,linf,npoints,path,aic_)
    print*,"Val subrutina AIC_VON",aic_
    !Validacion AIC_LOG(delta,r,sigma,npoints,path,aic)
    r =0.4
    call AIC_LOG(delta,r,sigma,npoints,path,aic_)
    print*,"Val subrutina AIC_LOG",aic_
    deallocate(path,pathkappahat,pathsigmahat)
end program main