program scattering
    use mod_NLS_2D
    use BiCGStab_mod
    ! C binding
    use, intrinsic :: iso_c_binding
    implicit NONE

    integer(C_INT), parameter           :: Nx = 256 !DEVE ESSERE UNA POTENZA DI 2!!!!!!!!
    integer(C_INT), parameter           :: Ny = 256
    real(C_DOUBLE), parameter           :: Lx = 100, Ly = 100 !32*pi
    real(C_DOUBLE), parameter           :: dx = Lx/Nx, dy = Ly/Ny, dkx=2*pi/Lx, dky=2*pi/Ly
    real(C_DOUBLE), dimension(Nx,Ny)    :: x,y,vx,vy,Jx,Jy,S,vx_c,vy_c
    real(C_DOUBLE)                      :: px,py
    real(C_DOUBLE), parameter           :: beta = 0.5, g=0.5
    complex(C_DOUBLE), dimension(Nx,Ny) :: dat,k1,k2,k3,k4
    real(C_DOUBLE), dimension(Nx,Ny)    :: kx,ky
    complex, parameter                  :: ii =(0.0,1.0)

    real(C_DOUBLE), dimension(Nx*Ny)    :: re,im
    real(C_DOUBLE), dimension(Nx,Ny)    :: re_mesh,im_mesh
    real(C_DOUBLE)                      :: v

    ! indices
    integer        :: i, j
    CHARACTER (len=150) :: msg
    character (len=2) :: a
    integer             :: ios
    real(C_DOUBLE) :: dt, tout,ttot
    integer        :: it,nt,nout,count
    real(C_DOUBLE) :: c,Natoms,E0,Ekin,Ekin_c,Eint,Eq,Etot
  
    ! Grid
    forall(i=1:Nx,j=1:Ny)
        x(i,j) = (i-1)*dx
        y(i,j) = (j-1)*dy
    end forall

    ! Compute the wavenumbers
    forall(i=1:Nx/2+1,j=1:Ny) kx(i,j)=(i-1)*dkx
    forall(i=Nx/2+2:Nx,j=1:Ny) kx(i,j)=(i-Nx-1)*dkx
    forall(i=1:Nx,j=1:Ny/2+1) ky(i,j)=(j-1)*dky
    forall(i=1:Nx,j=Ny/2+2:Ny) ky(i,j)=(j-Ny-1)*dky

    open(60,file='scattering_IC.txt',STATUS='unknown',ACTION='read',IOSTAT=ios,iomsg=msg)
    if (ios /= 0) then
        write(*,*) msg
        stop
    end if
    open(61,file='scatteringConservedQuantities.txt',STATUS='unknown',ACTION='write',IOSTAT=ios,iomsg=msg)
    if (ios /= 0) then
        write(*,*) msg
        stop
    end if
    open(62,file='scattering_vel.txt',STATUS='unknown',ACTION='write',IOSTAT=ios,iomsg=msg)
    if (ios /= 0) then
        write(*,*) msg
        stop
    end if
    open(63,file='scattering_NR.txt',STATUS='unknown',ACTION='write',IOSTAT=ios,iomsg=msg)
    if (ios /= 0) then
        write(*,*) msg
        stop
    end if
    open(70,file='scattering_GP.txt',STATUS='unknown',ACTION='write',IOSTAT=ios,iomsg=msg)
    if (ios /= 0) then
        write(*,*) msg
        stop
    end if

    !----------------  read input data -----------------
    WRITE(*,'(a,a)',advance='no') 'Verrà aperto: scattering_GP.txt'
    WRITE(*,*) ' '
    OPEN(UNIT=60,FILE='scattering_IC.txt',status='OLD',form='FORMATTED',access='SEQUENTIAL',action='READ',iostat=ios,iomsg=msg)
    openif: IF(ios/=0) THEN
        WRITE(*,'(a,1x,a)') 'Errore apertura file input: ',msg
        WRITE(*,*)
        STOP
    END IF openif

    i=1
    DO
        READ(60,'(F20.15, a, F20.15)',iostat=ios,iomsg=msg) re(i),a,im(i)
        IF (ios/=0) THEN
           ERROR STOP msg
        END IF
        i=i+1
        IF (i>Nx*Ny) EXIT
    END DO
    CLOSE(60)

    re_mesh = transpose(reshape(re,(/Ny,Nx/)))
    im_mesh = transpose(reshape(im,(/Ny,Nx/)))
    dat = (re_mesh + ii*im_mesh )
    
    write(*,*) 'Initial condition uploaded'
    !call velocity(dat,Nx,Ny,dx,dy,kx,ky,vx,vy,S)
    !do i=1,Nx
    !    do j=1,Ny
    !        write(62,'(F20.15, ",", F20.15, ",", F20.15)') vx(i,j),vy(i,j),S(i,j)
    !    end do
    !end do
    
    !-------------------------------- GP begins ------------------------------------
    write(*,*) '----- GPE begins -----'
    dt = 1.5d-2
    ttot = 3.0*4d1 
    tout = 2.0d1 ! ogni quanto voglio output 
    nt=int(ttot/dt)  !numero passi
    it=int(tout/dt)  !ogni it passi stampo output
    nout=int(nt/it) !numero output
    write(*,*) 'N steps= ', nt
    write(*,*) 'N output= ',nout

    Natoms=0
    Natoms = sum(abs(dat)**2)*dx*dy

    do count=1,nt
        !RK4 
        k1 = nlse(dat,Nx,Ny,kx,ky,beta,g)
        k2 = nlse(dat+0.5*k1*dt,Nx,Ny,kx,ky,beta,g)
        k3 = nlse(dat+0.5*k2*dt,Nx,Ny,kx,ky,beta,g)
        k4 = nlse(dat+k3*dt,Nx,Ny,kx,ky,beta,g)
        dat = dat + dt*(k1 + 2*k2 + 2*k3 + k4)/6   
        
        ! ------------------------------------------------------------
        !                    Conserved Quantities
        !N atoms
        c = sum(abs(dat)**2)*dx*dy
        !Energy evaluated by means of J
        call energy(dat,Nx,Ny,dx,dy,kx,ky,vx,vy,Ekin,Eint,Eq,Etot,Ekin_c)
        !             compilo file "scatteringConservedQuantities"
        write(61,'(F20.15, ",", F20.15, ",", F20.15, ",", F20.15, ",", F20.15)') Etot,Ekin,Eint,Eq,Ekin_c
        
        if(count.eq.1) then 
            E0=Etot
        end if
        if (count.eq.nt)then  !check if quantities are conserved in time evolution
            write(*,*) '***** Conserved quantities *****'
            write(*,*) 'N atoms=', Natoms,'; N atoms f', c 
            write(*,*) 'DeltaN= ', (c-Natoms)*100/Natoms,'%'
            write(*,*) 'E0=',E0,'; Ef=',Etot
            write(*,*) 'DeltaE=', (Etot-E0)*100/E0,'%'
        end if

        ! ------------------------------------------------------------------
        ! ------------------------------------------------------------------
        !                      vortices traking

        !here, for every time step (if(mod(count,it).eq.0)) i want to track the vortices
        !-------------------------- Newton-Raphson begins -----------------------------

        write(*,*) '----- Newton-Raphson begins -----'
        toler=1d-6
        mxmatvec=60

        vv = reshape(transpose(dat),(/Nx*Ny/))
        forall(i=1:Nx*Ny) vv(i) = dat!0.0626         !initial guess 

        do step = initialStep+1,finalStep+1

            !b = F(dat,Nx,Ny,kx,ky,beta,g,mu,v)  !RHS of the functional of which I want the root
            b = dat
            ! Check norm as exit condition via tolerance --> if F(Psi) = 0 I am in a Energy extrema 
            norma = dnorm2_bcg(Nx*Ny,reshape(b,(/Nx*Ny/))) 
            write(*,*) 'step',step,':', 'norm(F(psi))=', norma
            if (norma.le.(toler)) then
                write(*,*) 'convergence reached at previous step'
                !calcolo vorticity
                
                do i=1,Nx !compilo file con dati NR
                    do j=1,Ny
                        write(45,'(F20.15, ",", F20.15)') real(dat(i,j)), imag(dat(i,j))
                    end do
                end do
                
                go to 500
            end if

            ! ------------------- BiCGStab method to solve A*x = b ----------------------
            rhs = reshape(transpose(b),(/Nx*Ny/))
            call zbcg2 (.false.,2,Nx*Ny,vv,.true.,rhs,matvec,toler,mxmatvec,work,info,dat,Nx,Ny,kx,ky,beta,g,mu,v)
            if (info/=0) write(*,*) 'info',info

            dPsi = transpose(reshape(vv,(/Nx,Ny/)))
            dat = dat - dPsi

        end do
        500 continue

        ! -------------------------------------------------------------------------------
        !                     compilo file "Scattering_GP "con dati output
        do i=1,Nx        
            do j=1,Ny
                if(mod(count,it).eq.0) then
                    write(70,'(F20.15, ",", F20.15)') real(dat(i,j)),imag(dat(i,j))
                    write(62,'(F20.15, ",", F20.15)') vx(i,j),vy(i,j)
                end if
            end do
        end do
    end do

    close(60)
    close(61)
    close(62)
    close(63)
    close(70)

end program 