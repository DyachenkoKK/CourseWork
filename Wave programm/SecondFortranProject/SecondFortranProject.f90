    program SecondFortranProject
    
    use Globals
   
    
    implicit none
     real*8 pi
    complex*16 ci
    parameter (pi=3.141592653589793d0)
    parameter (ci = (0d0,1d0))
    ! Variables
    
    call ReadInput
    
    call openFirstResultRow
    
    call CalculateK
    
    call saveResults
    
    CONTAINS 
    
        function MakeSigma(kappa, alfa)
        real*8 kappa, alfa
        complex*16 MakeSigma
            if (alfa < kappa) then
                 MakeSigma = -ci*sqrt(kappa - (alfa**2))
            else 
                 MakeSigma = sqrt((alfa**2) - kappa)
            endif   
        end function MakeSigma     
        
        function MakeP(alpha, sgm1, sgm2, kappa2, delt)
        real*8 alpha, kappa2
        complex*16 MakeP, sgm1, sgm2, delt
            MakeP = 2*ci*myu*(-((alph**2) - 0.5d0 * kappa2) * exp(sgm1*currentZ) + sgm1 * sgm2 * exp(sgm2*currentZ))
        end function MakeP
        
        function MakeR(alpha, sgm1, kappa2, delt)
        real*8 alpha, kappa2
        complex*16 MakeR, sgm1, delt
            MakeR = 2*ci*myu*sgm1 * (-((alph**2) - 0.5d0 * kappa2) * exp(sgm1*currentZ) + (alph**2) * exp(sgm2*currentZ) ) / delt
        end function MakeR
        
        function MakeM(alpha, sgm2, kappa2, delt)
        real*8 alpha, kappa2
        complex*16 MakeM, sgm2, delt
            if(alph == 0 ) then 
                MakeM = 0
            else
                MakeM = 2*myu *sgm2* (-exp(sgm1*currentZ) + (1 - 0.5d0 * kappa2 / (alph**2)) * exp(sgm2*currentZ)) / delt 
            endif 
        end function MakeM
        
        function MakeS(alpha, sgm1, sgm2, kappa2, delt)
        real*8 alpha, kappa2
        complex*16 MakeS, sgm1, sgm2, delt
            MakeS = 2*myu* (-sgm1 * sgm2 * exp(sgm1*currentZ) + ((alph**2) - 0.5d0 * kappa2) * exp(sgm2*currentZ) ) / delt
        end function MakeS
        
        SUBROUTINE CalculateK
        integer i, j
            pointsNumber = ceiling((xLast-xFirst)/xStep)
            allocate(x(pointsNumber), z(pointsNumber), u(pointsNumber), q(k_size), tmp1(pointsNumber), tmp2(pointsNumber))
            call InitConst
            
            do i = 1, pointsNumber
                x(i) = xFirst+(i-1)*xStep
                z(i) = zFirst+(i-1)*zStep
            enddo
            
            do i = 1, pointsNumber
                alph=x(i)
                currentZ=z(i)
                call K(alph)
                tmp1(i)=k_arr(1,2)
                call K(-alph)
                tmp2(i)=k_arr(1,2)
            enddo
            
            t1 = kappa1/1d1; t2 = t1; t3 = t1; t4 = 1.4d0*kappa2+1d0;
            call dinn5(KIntegrand,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,u)
        END SUBROUTINE CalculateK
        
        SUBROUTINE InitConst
            q(1) = 0
            q(2) = 1
                        
            kappa1 = ro * (omega**2) / (lambd + 2 * myu)
            kappa2 = ro * (omega**2) / myu
        END SUBROUTINE InitConst
        
        SUBROUTINE KIntegrand(alfa, s, n)
        implicit none;
        integer n, i, j
        complex*16 alfa, s(n)
            do i = 1, n
                s(i) = tmp1(i) * exp(-ci*alfa*x(i)) 
                s(i) = s(i) + tmp2(i) * exp(ci*alfa*x(i))
                s(i) = s(i) / 2d0 * pi
            enddo
        END SUBROUTINE KIntegrand
        
        SUBROUTINE K(alpha)
        integer i, j
        real*8 alpha
        complex*16 P, R, M, S, N
            sgm1 = makeSigma(kappa1, alpha)
            sgm2 = makeSigma(kappa2, alpha)
            delt = 4*ci* (myu**2) * (-(((alpha**2) - 0.5d0 * kappa2)**2) + (alpha**2) * sgm1 * sgm2)
            P = makeP(alpha, sgm1, sgm2, kappa2, delt)
            R = makeR(alpha, sgm1, kappa2, delt)
            M = makeM(alpha, sgm2, kappa2, delt)
            S = makeS(alpha, sgm1, sgm2, kappa2, delt)
            k_arr(1,1) = -ci*(alpha**2)*M
            k_arr(1,2) = -ci*alpha*P
            k_arr(2,1) = alpha*S
            k_arr(2,2) = R
        END SUBROUTINE K
                
        SUBROUTINE openFirstResultRow
        integer i
            open(1, file='u.txt', FORM='FORMATTED')
            
            write(1,'(A)') "% Dinn5 settings , t1, t2, t3, t4, tm, tp, eps, step, IntLimit"
            write(1,'((A),9E15.6E3)') "% ", t1, t2, t3, t4, tm, tp, eps, step, IntLimit
            write(1,'(A)') "% 1) x, 2) z, 3) re(u), 4) im(u), 5) abs(u)"
            
            close(1); 
        END SUBROUTINE openFirstResultRow
        
        SUBROUTINE saveResults
        integer i
            open(1, file='u.txt', FORM='FORMATTED', position = "append")
            do i = 1, pointsNumber
                write(1, '(8E15.6E3)') x(i), z(i), real(u(i)), imag(u(i)), abs(u(i))
            enddo  
                
            close(1); 
        END SUBROUTINE saveResults
      
        
    end program SecondFortranProject

