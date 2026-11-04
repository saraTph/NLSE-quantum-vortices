MODULE BiCGStab_mod_Scatter
    ! C binding
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
    SAVE
    !          F for Newton-Raphson
    !          A for A*x = b
    !          BiCGStab subroutine

    CONTAINS
    !# Define function A
    subroutine matvec(vv,yy,dat,Nx,Ny,kx,ky,beta,g,mu,v)

        !CALL matvec(n,x,y)

        ! Define linear operator for bicgstab A
        !
        !                    dx RE[Psi]    dy Re[Psi]
        !                J = 
        !                    dx Im[Psi]    dy Im[psi]
        !
        ! A takes in input a vector vv = (nx*ny)
        ! returns in output a vector yy = (nx*ny)
        integer :: Nx,Ny,i,j,n,ios
        complex(C_DOUBLE), dimension(Nx*Ny), intent(out):: yy 
        complex(C_DOUBLE), dimension(Nx*Ny)             :: vv
        complex(C_DOUBLE), dimension(Nx,Ny)             :: dat,in,out,vv_mesh,vv_meshConj,outConj
        !real(C_DOUBLE), dimension(Nx,Ny)                :: u_adv,v_adv
        real(C_DOUBLE)                                  :: beta,g,v,mu
        real(C_DOUBLE), dimension(Nx,Ny)                :: kx,ky,mod2
        complex, parameter                              :: ii =(0.0,1.0)
        !---FFTW plans
        type(C_PTR)                                     :: pf,pb

        ! FFTW include
        include '/usr/local/include/fftw3.f03'

        !RESHAPE the real vector vv = 2*ntot = 2*(Nx*Ny) in a complex vv_mesh(Nx,Ny) 
        vv_mesh = transpose(reshape(vv,(/Nx,Ny/)))

        !evaluate derivatives: temp=-beta*(kx**2+ky**2)*temp + v*kx*temp 
        !fftw 
        ! Go to Fourier Space
        in = vv_mesh   !vv_mesh(Nx,Ny)
        pf = fftw_plan_dft_2d(Ny,Nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pf,in,out)
        ! Derivative
        forall(i=1:Nx,j=1:Ny) outConj(i,j) = conjg(out(i,j))
        out = (-beta)*(kx**2+ky**2)*out            !for single dipole: v=scalar 
        !out = (-beta)*(kx**2+ky**2)*out + (u_adv*(kx*out)) + (v_adv*(ky*out)) 
        ! Back to physical space
        pb = fftw_plan_dft_2d(Ny,Nx,out,in,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(pb,out,in)
        ! rescale
        in = in/Nx/Ny

        forall(i=1:Nx,j=1:Ny) vv_meshConj(i,j) = conjg(vv_mesh(i,j))
        in = in + (-2*g*mod2 + mu - (v*v)/(4*beta))*vv_mesh - g*(dat*dat)*vv_meshConj !sigle dipole
        !in = in + (-2*g*mod2 + mu - (sqrt(u_adv**2+v_adv**2))/(4*beta))*vv_mesh - g*(dat*dat)*vv_meshConj 

        !NOW reshape in(nx,ny) in a vector (nx*ny)
        yy = reshape(transpose(in),(/Nx*Ny/))

        call fftw_destroy_plan(pf)
        call fftw_destroy_plan(pb)
        return
    end subroutine

    function zdot_bcg(n,zx,zy)
        ! complex inner product function
        implicit none
        integer,intent(in)          :: n
        complex(C_DOUBLE)           :: zdot_bcg
        complex(C_DOUBLE),intent(in):: zx(n),zy(n)

        zdot_bcg = sum( conjg(zx) * zy )
    end function zdot_bcg

    function dnorm2_bcg(n,zx)
        ! l2 norm function
        implicit none
        integer,       intent(in)   :: n
        real(C_DOUBLE)              :: dnorm2_bcg
        complex(C_DOUBLE),intent(in):: zx(n)

        dnorm2_bcg = sqrt( abs( zdot_bcg(n, zx, zx) ))
    end function dnorm2_bcg

    subroutine zbcg2 (print_resid,l,n,x,nonzero_x,rhs,matvec,toler,mxmatvec,work,info,dat,Nx,Ny,kx,ky,bet,gnl,muu,vel)

! ------------------------------
! Description of the parameters:
! ------------------------------
!
! print_resid (input) LOGICAL. If print_resid=.true. the number of 
!            matrix-vector multiplications done so far and residual norm will
!            be printed to the standard output each iteration
!
! l          (input) INTEGER the dimension \ell of BiCGstab(\ell)
!            in this simple version it is required that l <= 2
!            l=2 is often useful for systems with nonsymmetric matrices
!
! n          (input) INTEGER size of the linear system to solve 
!
! x          (input/output) COMPLEX*16 array dimension n
!            initial guess on input, solution on output
!
! rhs        (input) COMPLEX*16 array dimension n
!            the right-hand side (r.h.s.) vector
!
! matvec     (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)
!
! nonzero_x  (input) LOGICAL tells
!            BiCGstab(\ell) if the initial guess x is zero or not. 
!            If nonzero_x is .FALSE., initial residual equals r.h.s. vector
!            and one MATVEC call is avoided
!
! toler      (input/output) DOUBLE PRECISION tolerance: the iterations are 
!            stopped as soon as || residual ||/|| initial residual|| <= toler,
!            the norm is Euclidean.  On output, if info>=0, the value of 
!            toler is set to the actually achieved residual reduction
!
! mxmatvec   (input/output) INTEGER.  On input: maximum number of matrix 
!            vector multiplications allowed to be done.  On output: 
!            if info>=0, mxmatvec is set to the actual number of matrix 
!            vector multiplications done
!
! work       (workspace) COMPLEX*16 array of dimension (n,2*l+5)
!
! info       (output) INTEGER.  info = 0 in case of succesful computations
!            and 
!            info = -m (<0) - means paramater number m has an illegal value
!            info = 1 - means no convergence achieved (stopping criterion
!            is not fulfilled)
!            info = 2 - means breakdown of the algorithm (taking a larger
!            value of parameter l usually helps)
!
! WARNING: If the iterations are ended normally (info=0 or info=1),
! the true residual norm is computed and returned as an output value 
! of the parameter toler.  The true residual norm can be slightly larger
! than the projected residual norm used by the algorithm to stop the
! iterations.  It may thus happen that on output info=0 but the value
! of toler is (slightly) larger than tolerance prescribed on input.
! ----------------------------------------------------------
        implicit none
        ! Parameters:
        logical,   intent(in)           :: print_resid,nonzero_x
        integer,   intent(in)           :: l, n
        integer,   intent(inout)        :: mxmatvec
        integer,   intent(out)          :: info
        complex(C_DOUBLE),intent(inout) :: x(n)
        complex(C_DOUBLE),intent(in)    :: rhs(n)
        real(C_DOUBLE),    intent(inout):: toler
        complex(C_DOUBLE),intent(out)   :: work(n,3+2*(l+1))

        ! Variables for matvec
        integer(C_INT)                      :: Nx,Ny
        real(C_DOUBLE)                      :: bet,gnl,muu,vel
        real(C_DOUBLE), dimension(Nx,Ny)    :: kx, ky
        complex(C_DOUBLE), dimension(Nx,Ny) :: dat
        !real(C_DOUBLE), dimension(Nx,Ny)    :: u_adv,v_adv

        ! Local variables:
        complex(C_DOUBLE)               :: matrix_z(l+1,l+1),y0(l+1),yl(l+1),zy0(l+1),zyl(l+1)
        logical                         :: rcmp, xpdt
        integer                         :: i, j, k, nmatvec
        complex(C_DOUBLE)               :: alpha, beta, omega, rho0, rho1, sigma
        complex(C_DOUBLE)               :: varrho, hatgamma
        real(C_DOUBLE)                  :: rnrm0, rnrm
        real(C_DOUBLE)                  :: mxnrmx, mxnrmr
        complex(C_DOUBLE)               :: kappa0, kappal

        ! Aliases for the parts of the work array:
        integer                     :: rr, r, u, xp, bp
        ! Constants:
        real(C_DOUBLE),parameter    :: delta = 1d-2
        complex(C_DOUBLE),parameter :: zzero = (0d0,0d0), zone = (1d0,0d0)

        ! Functions:
        external :: matvec

        info = 0
        if (l<1 .or. l>2) info = -2
        if (n<1) info = -3
        if (toler<=0d0) info = -9
        if (mxmatvec<0) info = -10

        rr = 1
        r = rr+1
        u = r+(l+1)
        xp = u+(l+1)
        bp = xp+1
        if (info/=0) return 

        ! Initialize first residual
        if (nonzero_x) then
            call matvec(x,work(1:n,r),dat,Nx,Ny,kx,ky,bet,gnl,muu,vel)  ! calcolo A*dx = A*x = work(1:n,r) 
            !call matvec(x,work(1:n,r),dat,Nx,Ny,kx,ky,bet,gnl,muu,u_adv,v_adv)
            work(1:n,r) = rhs - work(1:n,r)
            nmatvec = 1
        else
            work(1:n,r) = rhs
            nmatvec = 0
        end if
        !call precond (n, work(1:n,r))

        ! Initialize iteration loop
        work(1:n,rr) = work(1:n,r)
        work(1:n,bp) = work(1:n,r)
        work(1:n,xp) = x
        x = zzero

        rnrm0 = dnorm2_bcg(n, work(1:n,r))
        rnrm = rnrm0

        mxnrmx = rnrm0
        mxnrmr = rnrm0
        rcmp = .false.
        xpdt = .false.

        alpha = zzero
        omega = zone
        sigma = zone
        rho0  = zone

        ! Iterate
        do while (rnrm > toler*rnrm0 .and. nmatvec < mxmatvec)
            ! =====================
            ! The BiCG part ---
            ! =====================
            rho0 = -omega*rho0
            do k=1,l
                rho1 = zdot_bcg (n, work(1:n,rr), work(1:n,r+k-1))
                if (rho0.eq.zzero) then
                    info = 2
                    toler = rnrm/rnrm0
                    mxmatvec = nmatvec
                    return
                endif
                beta = alpha*(rho1/rho0)
                rho0 = rho1
                do j=0,k-1
                    work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
                end do 
                call matvec(work(1:n,u+k-1),work(1:n,u+k),dat,Nx,Ny,kx,ky,bet,gnl,muu,vel)
                !call matvec(work(1:n,u+k-1),work(1:n,u+k),dat,Nx,Ny,kx,ky,bet,gnl,muu,u_adv,v_adv)
                !call matvec (n, work(1:n,u+k-1), work(1:n,u+k))
                !call precond (n, work(1:n,u+k))
                nmatvec = nmatvec+1

                sigma = zdot_bcg (n, work(1:n,rr), work(1:n,u+k))
                if (sigma.eq.zzero) then
                    info = 2
                    toler = rnrm/rnrm0
                    mxmatvec = nmatvec
                    return
                endif
                alpha = rho1/sigma
                x(1:n) = alpha*work(1:n,u) + x(1:n)
                do j=0,k-1
                    work(1:n,r+j) = -alpha*work(1:n,u+j+1) + work(1:n,r+j)
                enddo
                call matvec(work(1:n,r+k-1),work(1:n,r+k),dat,Nx,Ny,kx,ky,bet,gnl,muu,vel)
                !call matvec(work(1:n,r+k-1),work(1:n,r+k),dat,Nx,Ny,kx,ky,bet,gnl,muu,u_adv,v_adv)
                !call matvec (n, work(1:n,r+k-1), work(1:n,r+k))
                !call precond (n, work(1:n,r+k))
                nmatvec = nmatvec+1
                rnrm = dnorm2_bcg (n, work(1,r))
                mxnrmx = max (mxnrmx, rnrm)
                mxnrmr = max (mxnrmr, rnrm)
            enddo

            ! ==================================
            ! The convex polynomial part ---
            ! ==================================

            !  --- Z = R'R

            do i=1,l+1
                do j=1,i
                    matrix_z(i,j) = conjg(zdot_bcg( n, work(1:n,r+j-1), work(1:n,r+i-1) ))
                end do
            end do
            
            !  lower triangular part of Z is computed; compute the rest knowing that Z^H=Z
            do j=2,l+1
                matrix_z(1:j-1,j) = conjg( matrix_z(j,1:j-1) )
            end do

            !  small vectors y0 and yl
            y0(1) = -zone
            y0(2) =      ( matrix_z(2,1) / matrix_z(2,2) )   ! works only for l=2
            y0(l+1) = zzero

            yl(1) = zzero
            yl(2) =      ( matrix_z(2,3) / matrix_z(2,2) )   ! works only for l=2
            yl(l+1) = -zone

            !  --- Convex combination

            ! compute Z*y0 and Z*yl
            zy0 = zzero
            zyl = zzero
            do j=1,l+1
                zy0 = zy0 + matrix_z(:,j)*y0(j)
                zyl = zyl + matrix_z(:,j)*yl(j)
            end do
            kappa0 = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )
            kappal = sqrt( abs(zdot_bcg(l+1,yl,zyl)) )

            varrho = zdot_bcg(l+1,yl,zy0) / (kappa0*kappal)
            hatgamma = varrho/abs(varrho) * max( abs(varrho),7d-1 ) 
            y0 = y0 - (hatgamma*kappa0/kappal)*yl
            
            !  --- Update
            omega = y0(l+1)
            do j=1,l
                work(1:n,u) = work(1:n,u) - y0(j+1)*work(1:n,u+j)
                x(1:n)      = x(1:n)      + y0(j+1)*work(1:n,r+j-1)
                work(1:n,r) = work(1:n,r) - y0(j+1)*work(1:n,r+j)
            enddo
            ! y0 has changed; compute Z*y0 once more
            zy0 = zzero
            do j=1,l+1
                zy0 = zy0 + matrix_z(:,j)*y0(j)
            end do
            rnrm = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )

            
            ! ================================
            ! The reliable update part ---
            ! ================================

            mxnrmx = max (mxnrmx, rnrm)
            mxnrmr = max (mxnrmr, rnrm)
            xpdt =  (rnrm < delta*rnrm0  .and. rnrm0 < mxnrmx)
            rcmp = ((rnrm < delta*mxnrmr .and. rnrm0 < mxnrmr) .or. xpdt)
            if (rcmp) then
                call matvec(x, work(1:n,r),dat,Nx,Ny,kx,ky,bet,gnl,muu,vel)
                !call matvec(x, work(1:n,r),dat,Nx,Ny,kx,ky,bet,gnl,muu,u_adv,v_adv)
                !call matvec (n, x, work(1:n,r))
                !call precond (n, work(1:n,r))
                nmatvec = nmatvec + 1
                work(1:n,r) =  work(1:n,bp) - work(1:n,r)
                mxnrmr = rnrm
                if (xpdt) then
                    work(1:n,xp) = x(1:n) + work(1:n,xp)
                    x = zzero
                    work(1:n,bp) = work(1:n,r)
                    mxnrmx = rnrm
                endif
            endif
            if (print_resid) print *,nmatvec,' ',log10(rnrm)
        enddo

        ! =========================
        ! End of iterations ---
        ! =========================

        x(1:n) = x(1:n) + work(1:n,xp)
        if (rnrm>toler*rnrm0) info = 1
        ! compute the true residual:
        ! --------------------- One matvec can be saved by commenting out this:
        call matvec(x, work(1:n,r),dat,Nx,Ny,kx,ky,bet,gnl,muu,vel)
        !call matvec(x, work(1:n,r),dat,Nx,Ny,kx,ky,bet,gnl,muu,u_adv,v_adv)
        !call matvec (n, x, work(1:n,r) )
        work(1:n,r) = rhs(1:n) - work(1:n,r) 
        !call precond (n, work(1:n,r))
        rnrm = dnorm2_bcg(n,work(1:n,r))
        nmatvec = nmatvec+1
        ! --------------------- One matvec can be saved by commenting out this^

        toler = rnrm/rnrm0
        mxmatvec = nmatvec
    end subroutine zbcg2
end MODULE BiCGStab_mod_Scatter