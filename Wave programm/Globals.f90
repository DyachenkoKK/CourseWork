MODULE Globals

    !material
    integer k_size
    real*8 alph, currentZ, kappa1, kappa2, lambd, myu, ro, omega
    complex*16 sgm1, sgm2, delt
    
    
    !study
    integer pointsNumber
    real*8 xFirst, xLast, xStep, zFirst, zLast, zStep
    complex*16, dimension(3,3) :: k_arr
    real*8, allocatable :: x(:), z(:), q(:)
    complex*16, allocatable :: u(:), tmp1(:), tmp2(:)
     
    !dinn5Settings
    real*8 t1, t2, t3, t4, tm, tp, eps, step, IntLimit
    
    namelist /material/ k_size, lambd, myu, ro, omega
    
    namelist /study/ xFirst, xLast, xStep, zFirst, zLast, zStep
    
    namelist /dinn5Settings/ t1, t2, t3, t4, tm, tp, eps, step, IntLimit
   
    CONTAINS
    
    SUBROUTINE ReadInput
        open(unit=1,file='input.txt',status='old')
        read(1, material) 
        read(1, study)
        read(1, dinn5Settings)
        close(1) 
    END SUBROUTINE ReadInput
    
    
END MODULE GLOBALS    