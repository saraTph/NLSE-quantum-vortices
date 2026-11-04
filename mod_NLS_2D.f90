MODULE mod_NLS_2D
    ! C binding
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
    SAVE

    ! Il modulo contiene
    ! 1) FFT
    ! 2) NLS functions for RK4
    !!!!!!!!!!!!!!!!!!remember o generate another module for initialConditions and conservedQuantities

    !dichiaro variabili
    integer, parameter       :: dp=selected_real_kind(15,300)
    real(kind=dp), parameter :: pi=3.141592653589793238460_dp
    real(kind=dp), parameter :: hbar = 6.6260715d-34/(2*pi)

    CONTAINS

    function VP(Nx,Ny,x,y,x0,y0,d,n)
        !geneate a vortex pair = (Nx x Ny): a straight vortex centered in (x0,y0) and with winding number n
        integer                            :: i,j,n,Nx,Ny
        real(C_DOUBLE)                     :: x0,y0,d
        real(C_DOUBLE), dimension(Nx,Ny)   :: x,y,r1,r2
        complex(C_DOUBLE), dimension(Nx,Ny):: VP !output
        real(C_DOUBLE)                     :: a1,a2,b1,b2
        real(C_DOUBLE), dimension(Nx,Ny)   :: theta1,theta2
        complex, parameter                 :: ii =(0.0,1.0)

        !constants from Padè approximation
        a1 = real(11)/32
        b1 = (5-32*a1)/(48-192*a1)
        a2 = a1*(b1-real(1)/4)
        b2 = a2

        r1 = sqrt((x-(x0))**2+(y-(y0-d/2))**2)  !centro del primo vortice (x0+d/2,y0+d/2)
        !r1 = sqrt((x-(x0))**2+(y-(y0))**2)  !centro del primo vortice (x0+d/2,y0+d/2)
        theta1 = atan2((x-(x0)),(y-(y0))) !questo è strano pero è cosi che li vuole Proment

        r2 = sqrt((x-(x0))**2+(y-(y0+d/2))**2)  !centro del secondo vortice (x0-d/2,y0-d/2)
        !theta1 = atan2((y-(y0-d/2)),(x-(x0)))  
        !theta2 = atan2((y-(y0+d/2)),(x-(x0)))
        theta1 = atan2((x-(x0)),(y-(y0-d/2))) !questo è strano pero è cosi che li vuole Proment
        theta2 = atan2((x-(x0)),(y-(y0+d/2)))
        
        VP = sqrt(r1**2*(a1 + a2*r1**2)/(1 + b1*r1**2 + b2*r1**4))*exp(ii*n*theta1)
        VP = VP * sqrt(r2**2*(a1 + a2*r2**2)/(1 + b1*r2**2 + b2*r2**4))*exp(-ii*n*theta2)
        return
    end function VP

    function NLSE(dat,Nx,Ny,kx,ky,beta,g) 
        integer                                     :: Nx,Ny,i,j
        complex(C_DOUBLE), dimension(Nx,Ny)         :: dat,in,der2,nlse
        real(C_DOUBLE), dimension(Nx,Ny)            :: kx, ky,mod2
        ! Fourier space variables
        complex(C_DOUBLE_COMPLEX), dimension(Nx,Ny) :: out 

        real(C_DOUBLE)     :: beta,g
        complex, parameter :: ii =(0.0,1.0)
        !---FFTW plans
        type(C_PTR) :: pf, pb

        ! FFTW include
        include '/usr/local/include/fftw3.f03'

        !fftw 
        ! Go to Fourier Space
        in = dat
        pf = fftw_plan_dft_2d(Ny,Nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in,out)
        ! Derivative
        out = -(kx**2+ky**2)*out !derivata seconda rispetto a x
        ! Back to physical space
        pb = fftw_plan_dft_2d(Ny,Nx,out,in,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out,in)
        ! rescale
        der2 = in/Nx/Ny

        !nlse
        mod2= abs(dat)**2
        nlse = -ii*(-beta*der2 + g*mod2*dat)
        
        call fftw_destroy_plan(pf)
        call fftw_destroy_plan(pb)
    end function nlse

    function NLSE_imTime(dat,Nx,Ny,kx,ky,beta,g,mu) 
        integer                                     :: Nx,Ny,i,j
        complex(C_DOUBLE), dimension(Nx,Ny)         :: dat,in,der2,nlse_imTime
        real(C_DOUBLE), dimension(Nx,Ny)            :: kx, ky,mod2
        ! Fourier space variables
        complex(C_DOUBLE_COMPLEX), dimension(Nx,Ny) :: out 

        real(C_DOUBLE)     :: beta,g,mu
        complex, parameter :: ii =(0.0,1.0)
        !---FFTW plans
        type(C_PTR) :: pf, pb

        ! FFTW include
        include '/usr/local/include/fftw3.f03'

        !fftw 
        ! Go to Fourier Space
        in = dat
        pf = fftw_plan_dft_2d(Ny,Nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in,out)
        ! Derivative
        out = -(kx**2+ky**2)*out !derivata seconda rispetto a x
        ! Back to physical space
        pb = fftw_plan_dft_2d(Ny,Nx,out,in,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out,in)
        ! rescale
        der2 = in/Nx/Ny

        !nlse
        mod2= abs(dat)**2
        nlse_imTime = (1-0.01*ii)*(beta*der2-g*mod2*dat+mu*dat) !imaginary time evolution  t -> -it
        
        call fftw_destroy_plan(pf)
        call fftw_destroy_plan(pb)
    end function nlse_imTime

    subroutine velocity(dat,Nx,Ny,dx,dy,kx,ky,vx,vy,S)
        !devo calcolare le due componenti di J e poi dividerle per rho
        ! J = Im[Psi* ATF[ik * TF[Psi]]]
        ! v = J/rho = J/mod2
        integer                                       :: Nx,Ny,i,j 
        real(C_DOUBLE)                                :: dx,dy
        real(C_DOUBLE), dimension(Nx,Ny), intent(OUT) :: vx,vy,S
        real(C_DOUBLE), dimension(Nx,Ny)              :: mod2,Jx,Jy
        complex(C_DOUBLE), dimension(Nx,Ny)           :: dat
        complex(C_DOUBLE), dimension(Nx,Ny)           :: in,in_x,in_y,out,out_x,out_y
        real(C_DOUBLE), dimension(Nx,Ny)              :: kx, ky
        real, parameter                               :: beta=0.5
        complex, parameter                            :: ii =(0.0,1.0)

        type(C_PTR) :: pf, pb
        ! FFTW include
        include '/usr/local/include/fftw3.f03'

        forall(i=1:Nx,j=1:Ny) S(i,j) = atan2(imag(dat(i,j)),real(dat(i,j)))

        !fftw 
        ! Go to Fourier Space
        in_x = dat
        in_y = dat
        pf = fftw_plan_dft_2d(Ny,Nx,in_x,out,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in_x,out)
        ! derive along x and y
        out_x = ii*kx*out
        out_y = ii*ky*out
        ! Back to physical space
        pb = fftw_plan_dft_2d(Ny,Nx,out_x,in_x,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out_x,in_x)
        pb = fftw_plan_dft_2d(Ny,Nx,out_y,in_y,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out_y,in_y)
        ! rescale
        in_x = in_x/Nx/Ny !ATF[i*kx * TF[Psi]]
        in_y = in_y/Nx/Ny !ATF[i*ky * TF[Psi]]

        forall(i=1:Nx,j=1:Ny) Jx(i,j) = (2*beta)*imag( conjg(dat(i,j))*in_x(i,j) )
        forall(i=1:Nx,j=1:Ny) Jy(i,j) = (2*beta)*imag( conjg(dat(i,j))*in_y(i,j) )

        mod2 = abs(dat)**2
        vx = Jx/(mod2+10d-16)
        vy = Jy/(mod2+10d-16)

        call fftw_destroy_plan(pf)
        call fftw_destroy_plan(pb)
    end subroutine velocity

    !I can evaluate the compressible and incompressible velocity by the Helmholz degomposition.
    subroutine velocity_c(dat,Nx,Ny,dx,dy,kx,ky,vx,vy,vx_c,vy_c)
        integer                                       :: Nx,Ny,i,j 
        real(C_DOUBLE)                                :: dx,dy
        real(C_DOUBLE), dimension(Nx,Ny)              :: vx,vy
        real(C_DOUBLE), dimension(Nx,Ny), intent(OUT) :: vx_c,vy_c
        real(C_DOUBLE), dimension(Nx,Ny)              :: mod2
        complex(C_DOUBLE), dimension(Nx,Ny)           :: dat
        complex(C_DOUBLE), dimension(Nx,Ny)           :: in_x,in_y,out_x,out_y
        real(C_DOUBLE), dimension(Nx,Ny)              :: kx, ky
        real, parameter                               :: beta=0.5
        complex, parameter                            :: ii =(0.0,1.0)

        type(C_PTR) :: pf, pb, pt,ps
        ! FFTW include
        include '/usr/local/include/fftw3.f03'

        !fftw 
        ! Go to Fourier Space
        in_x = vx
        in_y = vy
        pf = fftw_plan_dft_2d(Ny,Nx,in_x,out_x,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in_x,out_x)
        pb = fftw_plan_dft_2d(Ny,Nx,in_y,out_y,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,in_y,out_y)
        ! evaluate the scalar potential Gphi
        !Gphi = ii*(kx*out_x+ky*out_y)/(kx**2+ky**2)
        out_x = kx*(kx*out_x+ky*out_y)/(kx**2+ky**2+1.0d-10)  !vx_c(k) = ii*kx*Gphi
        out_y = ky*(kx*out_x+ky*out_y)/(kx**2+ky**2+1.0d-10)
        ! Back to physical space
        pt = fftw_plan_dft_2d(Ny,Nx,out_x,in_x,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pt,out_x,in_x)
        ps = fftw_plan_dft_2d(Ny,Nx,out_y,in_y,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(ps,out_y,in_y)
        ! rescale
        vx_c = sqrt(real(in_x)**2+imag(in_x)**2)/Nx/Ny !x component of the compressible velocity
        vy_c = sqrt(real(in_x)**2+imag(in_x)**2)/Nx/Ny !y component

        call fftw_destroy_plan(pf)
        call fftw_destroy_plan(pb)
        call fftw_destroy_plan(ps)
        call fftw_destroy_plan(pt)
    end subroutine velocity_c


    subroutine energy(dat,Nx,Ny,dx,dy,kx,ky,vx,vy,Ekin,Eint,Eq,Etot,Ekin_c,Etot_class)
        !devo calcolare le due componenti di J e poi dividerle per rho
        ! J = Im[Psi* ATF[ik * TF[Psi]]]
        ! v = J/rho = J/mod2
        integer                                       :: Nx,Ny,i,j
        real(C_DOUBLE)                                :: dx,dy,N
        real(C_DOUBLE), intent(OUT)                   :: Ekin,Eint,Eq,Etot,Ekin_c,Etot_class
        real(C_DOUBLE), dimension(Nx,Ny)              :: mod2,vx,vy
        complex(C_DOUBLE), dimension(Nx,Ny)           :: dat,datStar,Jx,Jy,Jx_c,Jy_c
        complex(C_DOUBLE), dimension(Nx,Ny)           :: in_x,in_y,out,out_x,out_y
        real(C_DOUBLE), dimension(Nx,Ny)              :: kx, ky
        real, parameter                               :: beta=0.5, g=0.5
        complex, parameter                            :: ii =(0.0,1.0)

        type(C_PTR) :: pf,pb,ps,pt
        ! FFTW include
        include '/usr/local/include/fftw3.f03'

        !fftw 
        ! Go to Fourier Space
        in_x = dat
        in_y = dat
        pf = fftw_plan_dft_2d(Ny,Nx,in_x,out,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in_x,out)
        ! derive along x and y
        out_x = ii*kx*out
        out_y = ii*ky*out

        Etot_class = 0.5*sum((kx**2+ky**2)*abs(out)**2)*dx*dy
        
        ! Back to physical space
        pb = fftw_plan_dft_2d(Ny,Nx,out_x,in_x,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out_x,in_x)
        pb = fftw_plan_dft_2d(Ny,Nx,out_y,in_y,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out_y,in_y)

        ! rescale
        in_x = in_x/Nx/Ny !ATF[i*kx * TF[Psi]]
        in_y = in_y/Nx/Ny !ATF[i*ky * TF[Psi]]

        mod2 = abs(dat)**2
        forall(i=1:Nx,j=1:Ny) datStar(i,j) = conjg(dat(i,j))

        Jx = (datStar*in_x)/sqrt(datStar*dat+10d-16)                !J complex
        Jy = (datStar*in_y)/sqrt(datStar*dat+10d-16)

        Ekin = 2*beta**2 * sum((imag(Jx))**2+(imag(Jy))**2)*dx*dy !kinetic energy
        Eq = 2*beta**2 * sum((real(Jx))**2+(real(Jy))**2)*dx*dy   !quantistic energy
        Eint = g/2*sum(mod2**2)*dx*dy                             !interaction energy
        Etot = Ekin + Eq + Eint

        ! --------------------------------------------------------------------------
        ! --------------------------------------------------------------------------
        ! decompongo Jx e Jy in una parte comprimibile e non: Helmholtz dcomposition
        !fftw 
        ! Go to Fourier Space
        in_x = Jx
        in_y = Jy
        pf = fftw_plan_dft_2d(Ny,Nx,in_x,out_x,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in_x,out_x)
        pb = fftw_plan_dft_2d(Ny,Nx,in_y,out_y,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,in_y,out_y)
        ! evaluate the scalar potential Gphi
        !Gphi = ii*(kx*out_x+ky*out_y)/(kx**2+ky**2)
        out_x = kx*(kx*out_x+ky*out_y)/(kx**2+ky**2+1.0d-10)  !vx_c(k) = ii*kx*Gphi
        out_y = ky*(kx*out_x+ky*out_y)/(kx**2+ky**2+1.0d-10)
        ! Back to physical space
        pt = fftw_plan_dft_2d(Ny,Nx,out_x,in_x,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pt,out_x,in_x)
        ps = fftw_plan_dft_2d(Ny,Nx,out_y,in_y,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(ps,out_y,in_y)
        ! rescale
        Jx_c = (in_x)/Nx/Ny !x component of the compressible velocity
        Jy_c = (in_y)/Nx/Ny !y component

        Ekin_c = 2*beta**2 * sum((imag(Jx_c))**2+(imag(Jy_c))**2)*dx*dy 

        call fftw_destroy_plan(pf)
        call fftw_destroy_plan(pb)
        call fftw_destroy_plan(ps)
        call fftw_destroy_plan(pt)
    end subroutine energy

    subroutine momentum(dat,Nx,Ny,dx,dy,beta,kx,ky,px,py)
        !devo calcolare le due componenti di J e poi dividerle per rho
        ! J = Im[Psi* ATF[ik * TF[Psi]]]
        ! v = J/rho = J/mod2
        integer                                       :: Nx,Ny,i,j
        real(C_DOUBLE), intent(OUT)                   :: px,py
        complex(C_DOUBLE), dimension(Nx,Ny)           :: dat
        complex(C_DOUBLE), dimension(Nx,Ny)           :: in,in_x,in_y,out,out_x,out_y
        real(C_DOUBLE), dimension(Nx,Ny)              :: kx, ky
        real(C_DOUBLE)                                :: beta,dx,dy
        complex, parameter                            :: ii =(0.0,1.0)

        type(C_PTR) :: pf, pb
        ! FFTW include
        include '/usr/local/include/fftw3.f03'

        !fftw 
        ! Go to Fourier Space
        in_x = dat
        in_y = dat
        pf = fftw_plan_dft_2d(Ny,Nx,in_x,out,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in_x,out)
        ! px=beta*sum(kx*abs(psi)**2)*dV
        px=beta*sum(kx*abs(out)**2)*(dx*dy)
        py=beta*sum(ky*abs(out)**2)*(dx*dy)
        ! Back to physical space
        pb = fftw_plan_dft_2d(Ny,Nx,out,in_x,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out,in_x)
        ! rescale
        in_x = in_x/Nx/Ny !ATF[i*kx * TF[Psi]]

        call fftw_destroy_plan(pf)
        call fftw_destroy_plan(pb)
    end subroutine momentum


    function shift_NR(dat,Nx,Ny,dx,dy,kx,ky,deltax,deltay)
        integer                                       :: Nx,Ny
        complex(C_DOUBLE), dimension(Nx,Ny)           :: dat,shift_NR
        complex(C_DOUBLE), dimension(Nx,Ny)           :: in,out
        real(C_DOUBLE), dimension(Nx,Ny)              :: kx, ky
        real(C_DOUBLE)                                :: dx,dy,deltax,deltay
        complex, parameter                            :: ii =(0.0,1.0)

        type(C_PTR) :: pf, pb
        ! FFTW include
        include '/usr/local/include/fftw3.f03'
        !fftw 
        ! Go to Fourier Space
        in = dat
        pf = fftw_plan_dft_2d(Ny,Nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in,out)
        out=out*exp(ii*(deltax*kx+deltay*ky))
        ! Back to physical space
        pb = fftw_plan_dft_2d(Ny,Nx,out,in,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out,in)
        ! rescale
        shift_NR = in/Nx/Ny !ATF[i*kx * TF[Psi]]

        call fftw_destroy_plan(pf)
        call fftw_destroy_plan(pb)
    end function

    function vorticity_ps(dat,Nx,Ny,dx,dy,kx,ky)
        integer                                       :: Nx,Ny,i,j
        complex(C_DOUBLE), dimension(Nx,Ny)           :: dat
        real(C_DOUBLE), dimension(Nx,Ny)              :: vorticity_ps
        complex(C_DOUBLE), dimension(Nx,Ny)           :: in_x,in_y,out_x,out_y
        real(C_DOUBLE), dimension(Nx,Ny)              :: dx_re,dx_im,dy_re,dy_im
        real(C_DOUBLE), dimension(Nx,Ny)              :: vorticity_ps_x,vorticity_ps_y
        real(C_DOUBLE), dimension(Nx,Ny)              :: kx, ky
        real(C_DOUBLE)                                :: dx,dy
        complex, parameter                            :: ii =(0.0,1.0)

        type(C_PTR) :: pf, pb
        ! FFTW include
        include '/usr/local/include/fftw3.f03'
        !fftw 
        ! Go to Fourier Space
        in_x = dat
        in_y = dat
        pf = fftw_plan_dft_2d(Ny,Nx,in_x,out_x,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in_x,out_x)
        out_x=ii*kx*out_x
        pf = fftw_plan_dft_2d(Ny,Nx,in_y,out_y,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in_y,out_y)
        out_y=ii*ky*out_y
        ! Back to physical space
        pb = fftw_plan_dft_2d(Ny,Nx,out_x,in_x,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out_x,in_x)
        pb = fftw_plan_dft_2d(Ny,Nx,out_y,in_y,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out_y,in_y)
        ! rescale
        out_x=out_x/Nx/Ny
        out_y=out_y/Nx/Ny

        do i=1,Nx
            do j=1,Ny
                dx_re(i,j)=real(out_x(i,j))
                dx_im(i,j)=imag(out_x(i,j))
                dy_re(i,j)=real(out_y(i,j))
                dy_im(i,j)=imag(out_y(i,j))
            end do
        end do

        vorticity_ps_x = dx_re*dy_im
        vorticity_ps_y = -dx_im*dy_re
        vorticity_ps =sqrt(vorticity_ps_x**2+vorticity_ps_y**2)

        call fftw_destroy_plan(pf)
        call fftw_destroy_plan(pb)
    end function

    subroutine psi_part_sub(psi_part,x_cord, y_cord, complex_matrix, real_or_img)
        !integer, parameter:: GP = KIND(0.0D0)
        COMPLEX(C_DOUBLE), DIMENSION(256,256):: complex_matrix
        real(C_DOUBLE) :: psi_part
        integer :: x_cord, y_cord, final_x, final_y, real_or_img
        
        final_x = x_cord
        final_y = y_cord
        
        if(x_cord > 256) then
            final_x = 256
        else if (x_cord < 1) then
            final_x = 1
        endif
    
        if(y_cord > 256) then
            final_y = 256
        else if (y_cord < 1) then
            final_y = 1
        endif
    
        ! write(*,*) "final_x, final_y",final_x, final_y
        ! write(*,*) real_or_img
        if(real_or_img .eq. 1) then
            psi_part = real(complex_matrix(final_x, final_y))
        else
            psi_part = aimag(complex_matrix(final_x, final_y))
        endif
        
    endsubroutine psi_part_sub

    subroutine deriv(Nx,Ny,kx,ky,dat,deriv_xPsi,deriv_yPsi)

        integer, intent(IN)                            :: Nx,Ny
        complex(c_DOUBLE), DIMENSION(Nx,Ny), intent(IN):: dat
        real(C_DOUBLE), dimension(Nx,Ny),intent(IN)    :: kx,ky
        complex(C_DOUBLE), DIMENSION(Nx,Ny),intent(OUT):: deriv_xPsi,deriv_yPsi

        !fftw variables
        complex(C_DOUBLE), DIMENSION(Nx,Ny) :: in,in_x,in_y,out,out_x,out_y
        complex, parameter                  :: ij =(0.0,1.0)
        !---FFTW plans
        type(C_PTR)                         :: pf,pb,pt
        ! FFTW include
        include '/usr/local/include/fftw3.f03'

        ! -------------------------------------------------------------------
        !                 evaluate derivatives
    
        !fftw 
        ! Go to Fourier Space
        in = dat
        in_x = dat 
        in_y = dat
        pf = fftw_plan_dft_2d(Ny,Nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in,out)
        ! Derivative
        out_x = ij*kx*out 
        out_y = ij*ky*out
        ! Back to physical space
        pb = fftw_plan_dft_2d(Ny,Nx,out_x,in_x,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out_x,in_x)
        pt = fftw_plan_dft_2d(Ny,Nx,out_y,in_y,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pt,out_y,in_y)
        ! rescale
        deriv_xPsi = in_x/(Nx*Ny)
        deriv_yPsi = in_y/(Nx*Ny)
        call fftw_destroy_plan(pf)
        call fftw_destroy_plan(pb)
        call fftw_destroy_plan(pt)
        
    end subroutine

    subroutine tracking_2d(Nx,Ny,dx,dy,kx,ky,dat,x_guess,y_guess,x_vortex,y_vortex,x_v,y_v)

        !integer, parameter:: GP = KIND(0.0D0)
        integer, intent(OUT)                           :: x_vortex,y_vortex !grid point vortex coordinates
        real(C_DOUBLE), intent(OUT)                    :: x_v, y_v !real vortex coordinates
        real(C_DOUBLE)                                 :: delta_x,delta_y
        real(C_DOUBLE)                                 :: dx,dy
        integer, intent(IN)                            :: x_guess, y_guess, Nx,Ny
        complex(c_DOUBLE), DIMENSION(Nx,Ny), intent(IN):: dat
        complex(c_DOUBLE), DIMENSION(Nx,Ny)            :: dat_shifted
        real(C_DOUBLE), dimension(Nx,Ny),intent(IN)    :: kx,ky
        complex(C_DOUBLE), DIMENSION(Nx,Ny)            :: deriv_xPsi,deriv_yPsi

        !iternal variables
        integer :: h
        logical :: NR_stop, should_stop
        real(C_DOUBLE) :: dau_x_psi_r, dau_x_psi_i, dau_y_psi_r, dau_y_psi_i, psi_part,det, psi_x_guess_real, psi_x_guess_img
        real :: condensate_density, mass_of_boson, mod_psi_inf, delta
        real(C_DOUBLE) :: w_ps
        real(C_DOUBLE), dimension(2,2) :: j_inverse

        !h = 1
        NR_stop = .false.
        ! for stopping condition, TODO: Need to know proper values
        condensate_density = 1
        mass_of_boson = 1
        mod_psi_inf = (condensate_density/mass_of_boson)**0.5
        delta = 0.1

        x_vortex = x_guess
        y_vortex = y_guess
        do while(NR_stop .eqv. .FALSE.)
            ! 1 for real, 2 for imaginary
            call deriv(Nx,Ny,kx,ky,dat,deriv_xPsi,deriv_yPsi)
            dau_x_psi_r = real(deriv_xPsi(x_vortex,y_vortex))
            dau_x_psi_i = imag(deriv_xPsi(x_vortex,y_vortex))
            dau_y_psi_r = real(deriv_yPsi(x_vortex,y_vortex))  
            dau_y_psi_i = imag(deriv_yPsi(x_vortex,y_vortex))

            det = dau_x_psi_r*dau_y_psi_i - dau_x_psi_i*dau_y_psi_r

            j_inverse(1,1) = dau_y_psi_i/det
            j_inverse(1,2) = -dau_y_psi_r/det
            j_inverse(2,1) = -dau_x_psi_i/det
            j_inverse(2,2) = dau_x_psi_r/det
            !write(*,*) "Jacobian Inverse in the loop is this ", j_inverse

            call psi_part_sub(psi_x_guess_real,x_vortex,y_vortex,dat,1)
            call psi_part_sub(psi_x_guess_img,x_vortex,y_vortex,dat,2)
            !write(*,*) x_vortex
            !write(*,*) x_vortex - (j_inverse(1,1)*psi_x_guess_real + j_inverse(1,2)*psi_x_guess_img)
            !write(*,*) int(x_vortex - (j_inverse(1,1)*psi_x_guess_real + j_inverse(1,2)*psi_x_guess_img))

            x_vortex = nint(x_vortex - (j_inverse(1,1)*psi_x_guess_real + j_inverse(1,2)*psi_x_guess_img))
            y_vortex = nint(y_vortex - (j_inverse(2,1)*psi_x_guess_real + j_inverse(2,2)*psi_x_guess_img))

            write(*,*) "***end of loop GRID-x,y", x_vortex,y_vortex

            ! suppose at end of loop it jumps away
            if(x_vortex > 256 .or. x_vortex < 1) then
                x_vortex = 999
                goto 99 
            endif

            if(y_vortex > 256 .or. y_vortex < 1) then
                y_vortex = 999
                goto 99 
            endif

            !write(*,*) ">>>>mod_psi",ABS(psi(x,y)**2)
            call sleep(1)
            ! w_ps = h_bar * grad_psi_real * grad_psi_img (take h_bar = 1)
            w_ps = dau_x_psi_r*dau_y_psi_i - dau_y_psi_r*dau_x_psi_i 
            ! condition to stop: mod_psi < delta x mod_psi_inf
            if((ABS(dat(x_vortex,y_vortex)) - delta*mod_psi_inf <= 0.23) .and. w_ps /= 0) then
                write(*,*) 'Vortex grid point found'
                NR_stop = .true.
            endif
        end do

        ! ****************************************************************
        ! Fin qui ho trovato il punto griglia che contine il vortex point
        !           Now search the more precise voretex point
        NR_stop = .false.
        x_v = x_vortex !now I modify x_v and y_v
        y_v = y_vortex
        dat_shifted = dat
        do while(NR_stop .eqv. .FALSE.)

            call deriv(Nx,Ny,kx,ky,dat_shifted,deriv_xPsi,deriv_yPsi)
            dau_x_psi_r = real(deriv_xPsi(x_vortex,y_vortex))
            dau_x_psi_i = imag(deriv_xPsi(x_vortex,y_vortex))
            dau_y_psi_r = real(deriv_yPsi(x_vortex,y_vortex))  
            dau_y_psi_i = imag(deriv_yPsi(x_vortex,y_vortex))

            det = dau_x_psi_r*dau_y_psi_i - dau_x_psi_i*dau_y_psi_r

            j_inverse(1,1) = dau_y_psi_i/det
            j_inverse(1,2) = -dau_y_psi_r/det
            j_inverse(2,1) = -dau_x_psi_i/det
            j_inverse(2,2) = dau_x_psi_r/det
            !write(*,*) "Jacobian Inverse in the loop is this ", j_inverse

            call psi_part_sub(psi_x_guess_real,x_vortex,y_vortex,dat,1)
            call psi_part_sub(psi_x_guess_img,x_vortex,y_vortex,dat,2)

            delta_x = (j_inverse(1,1)*psi_x_guess_real + j_inverse(1,2)*psi_x_guess_img)
            delta_y = (j_inverse(2,1)*psi_x_guess_real + j_inverse(2,2)*psi_x_guess_img)
            x_v = x_v - delta_x
            y_v = y_v - delta_y

            ! suppose at end of loop it jumps away
            if(x_v > 256 .or. x_v < 1) then
                x_v = 999
                goto 99 
            endif

            if(y_v > 256 .or. y_v < 1) then
                y_v = 999
                goto 99 
            endif

            write(*,*) "***end of loop-x_v,y_v", x_v,y_v
            call sleep(1)

            dat_shifted = shift_NR(dat_shifted,Nx,Ny,dx,dy,kx,ky,-delta_x,-delta_y)
            !now if I evalue dat_shifted in x_vortex, y_vortex, I am really valuating dat in x_v,y_v
            w_ps = dau_x_psi_r*dau_y_psi_i - dau_y_psi_r*dau_x_psi_i 
            ! condition to stop: mod_psi < delta x mod_psi_inf

            !if((ABS(dat_shifted(x_vortex,y_vortex))- delta*mod_psi_inf <= 0.01) .and. w_ps /= 0) then
            if((ABS(dat_shifted(x_vortex,y_vortex))- delta*mod_psi_inf  <= 0.01) .and. w_ps /= 0) then
                !write(*,*) ABS(dat_shifted(x_vortex,y_vortex))**2
                NR_stop = .true.
            endif
        end do
        99 NR_stop = .true.
    end subroutine

end module mod_NLS_2D