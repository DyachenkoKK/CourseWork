module globals
implicit none

integer i, pointsNumber

real*8 lamda, mu, rho
real*8 z, xmin, xstep, xmax, w, pi, f
real*8 kappa(2)
complex*16 sigma(2), ci
parameter (pi=3.141592653589793d0)
parameter (ci=(0d0,1d0))

real*8, allocatable:: x(:)
complex*16, allocatable:: u(:)
real*8 t1, t2, t3, t4, tm, tp, eps, step, IntLimit

    contains
    
    subroutine InitGlobals
    implicit none
        lamda = 55.5d6
        mu = 26.1d6
        rho = 2.698d0
        z = 0d0
        xmin = 0d0
        xstep = 1d0 
        xmax = 140d0
        f = 200d0
        
        
        w = 2d0*pi*f
        kappa(1) = sqrt(rho*w**2/(lamda+2d0*mu))
        kappa(2) = sqrt(rho*w**2/(mu))
        
        tm = 2d-2
        tp = 0d0 
        eps = 1d-6 
        step = 1d-2 
        IntLimit = 1d10
        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0)
        
        pointsNumber = ceiling((xMax-xMin)/xStep)
        allocate(x(pointsNumber), u(pointsNumber))
        do i = 1, pointsNumber
            x(i) = xMin + xStep*(i-1)  
        enddo
        
        
    end subroutine InitGlobals
    
 end module globals