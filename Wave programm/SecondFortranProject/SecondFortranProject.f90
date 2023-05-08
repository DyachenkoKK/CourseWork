
program main
use globals
implicit none

    



    call InitGlobals
    call dinn5(u_integrand,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,u)
    call PlotResults
    
    
    
CONTAINS
       
    FUNCTION makeSigma(kappa, alfa)
    real*8 kappa(2)
    complex*16 ci, alfa, makeSigma(2)
    parameter (ci = (0d0,1d0))
        if (Imag(alfa) == 0d0)  then 
                if (abs(alfa) < Kappa(1)) then
                    MakeSigma(1) = -ci*sqrt(Kappa(1)**2 - real(alfa)**2)
                else 
                    MakeSigma(1) = sqrt(real(alfa)**2 - Kappa(1)**2)
                endif    
                    
                if (abs(alfa) < Kappa(2)) then
                    MakeSigma(2) = -ci*sqrt(Kappa(2)**2 - real(alfa)**2)
                else 
                    MakeSigma(2) = sqrt(real(alfa)**2 - Kappa(2)**2)
                endif                
        else 
                MakeSigma(1) = sqrt(alfa**2 - cmplx(Kappa(1)**2))
                MakeSigma(2) = sqrt(alfa**2 - cmplx(Kappa(2)**2))
        endif  
    END   
    
    
    FUNCTION makeDelta(alfa, mu, kappa)
    ! single mu!!!
    real*8 mu, kappa(2)
    complex*16 alfa, makeDelta, sigma(2)
        sigma = makeSigma(kappa, alfa)
        makeDelta = 4d0*ci*mu*(-(alfa**2 - 0.5d0*kappa(2)**2)**2 + alfa**2*sigma(1)*sigma(2)) 
    END
    
    
    FUNCTION MakeP(alfa, mu, kappa)
    ! single mu!!!
    real*8 mu, kappa(2)
    complex*16 alfa, MakeP, sigma(2), delta
        sigma = makeSigma(kappa, alfa)
        delta = makeDelta(alfa, mu, kappa)
        MakeP =2d0*ci*(-(alfa**2-0.5d0*kappa(2)**2)*exp(sigma(1)*z) + sigma(1)*sigma(2)*exp(sigma(2)*z))/delta
    END
    
    
    
    SUBROUTINE u_integrand(alfa, s, n)
    implicit none;
    integer n, i
    complex*16 alfa, s(n), sigma(2), P, Pminus 
        sigma = makeSigma(kappa, alfa)
        P = MakeP(alfa, mu, kappa)
        Pminus = MakeP(-alfa, mu, kappa)
        do i = 1, n
            s(i) = -ci*alfa*P*exp(-ci*alfa*x(i)) + alfa*ci*Pminus*exp(ci*alfa*x(i))   
            s(i) = s(i)/(2d0*pi)
        enddo     
    END SUBROUTINE u_integrand
    
    
    subroutine PlotResults 
    implicit none
        open(1, file='u.txt', FORM='FORMATTED')
        do i = 1, pointsNumber
            write(1, '(2E15.6E3)'), x(i), abs(u(i))
        enddo
    end subroutine PlotResults
        
    
    

    end program main
