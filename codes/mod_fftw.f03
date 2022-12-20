module mod_fftw
use mod_param
use mod_operator
use omp_lib
use mod_io
use, intrinsic :: iso_c_binding
include 'fftw3.f03'
contains

  subroutine fftw_forward_run(M,N,arrayin,arrayout)
    integer, intent(in)::M,N
    real(PP),intent(in):: arrayin(M,N)
    complex(PP),intent(out):: arrayout(M/2 + 1, N)
    integer :: plan(8)

    call dfftw_plan_dft_r2c_2d(plan,M,N,arrayin,arrayout,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan,arrayin,arrayout)
    call dfftw_destroy_plan(plan)
  end subroutine fftw_forward_run

  subroutine fftw_backward_run(M,N,arrayin,arrayout)
    integer, intent(in)::M,N
    complex(PP),intent(in):: arrayin(M/2+1,N)
    real(PP),intent(out):: arrayout(M, N)
    integer :: plan(8)
    call dfftw_plan_dft_c2r_2d(plan,M,N,arrayin,arrayout,FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(plan,arrayin,arrayout)
    call dfftw_destroy_plan(plan)
  end subroutine fftw_backward_run

  subroutine fftw3_laplace
    integer :: i,j,k
    real(PP) :: L,zv,zf
    real(PP) :: bz0mean
    real(PP), dimension(:,:), allocatable :: kxarr,kyarr,qarr,qinvarr,farr,kappa,invkappa
    real(PP), dimension(:,:), allocatable :: farrz
    complex(PP), dimension(:,:), allocatable :: fbz0
    complex(PP), dimension(:,:,:,:), allocatable :: fbik
    complex(PP), dimension(:,:,:,:), allocatable :: fbiktest
    complex(PP), dimension(:,:), allocatable :: fbzL
    integer,dimension(8)::plan_ax,plan_ay,plan_az
    integer,dimension(8)::plan_bx,plan_by,plan_bz
    real(PP)::dkx,dky
    ! Calculate potential field for given bz0, via Fourier solution 
    ! (e.g. Gary 1989, ApJ Supp. 69, 323). 
    ! Allocate memory
    allocate(kxarr(dimx/2+1,dimy),kyarr(dimx/2+1,dimy))
    allocate(qarr(dimx/2+1,dimy),qinvarr(dimx/2+1,dimy),farr(dimx/2+1,dimy))
    allocate(farrz(dimx/2+1,dimy))
    allocate(fbz0(dimx/2+1,dimy),fbik(dimx/2+1,dimy,dimz,3),fbzL(dimx/2+1,dimy))
    allocate(fbiktest(dimx/2+1,dimy,dimz,3))
    allocate(kappa(dimx/2+1,dimy))
    allocate(invkappa(dimx/2+1,dimy))
    fbik=0.0d0
    call dfftw_plan_dft_c2r_2d(plan_bx,dimx,dimy,fbik(:,:,1,1),bxyz0(:,:,1,1),FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_2d(plan_by,dimx,dimy,fbik(:,:,1,2),bxyz0(:,:,1,2),FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_2d(plan_bz,dimx,dimy,fbik(:,:,1,3),bxyz0(:,:,1,3),FFTW_ESTIMATE)

    if(Aout) then
      fap=0.0d0
      call dfftw_plan_dft_c2r_2d(plan_ax,dimx,dimy,fap(:,:,1,1),axyz(:,:,1,1),FFTW_ESTIMATE)
      call dfftw_plan_dft_c2r_2d(plan_ay,dimx,dimy,fap(:,:,1,2),axyz(:,:,1,2),FFTW_ESTIMATE)
      call dfftw_plan_dft_c2r_2d(plan_az,dimx,dimy,fap(:,:,1,3),axyz(:,:,1,3),FFTW_ESTIMATE)
    end if

    qarr=real(0,PP)
    qinvarr=real(0,PP)
    farr=real(0,PP)
    farrz=real(0,PP)
    bzL=real(0,PP)
    bz0t=real(0,PP)
    ! Define the flux-balanced top BC
    bz0mean=sum(bz0)/size(bz0)
    bzL(:,:)=bz0mean
    bz0t=bz0-bzL
    if (top_closed) then
      bz0n=bz0t 
    else
      bz0n=bz0
    endif
    ! Fourier transform boundary data
    call fftw_forward_run(dimx,dimy,bz0,fbz0)
    ! Top of the box
    L=zarr(dimz)
    ! Define wave number arrays for Fourier components in the arrays

! ------------------------- may 3 add, start ----------------------
    dkx=real(1,PP)/(dimx*dx)
    dky=real(1,PP)/(dimy*dy)

    forall(i=1:dimx/2+1,j=1:dimy) kxarr(i,j)=real(i-1,PP)*dkx
    forall(i=1:dimx/2+1,j=1:dimy/2) kyarr(i,j)=real(j-1,PP)*dky
    forall(i=1:dimx/2+1,j=dimy/2+1:dimy) kyarr(i,j)=real(j-dimy/2-1,PP)*dky-real(0.5*dimy,PP)*dky

    kxarr=-kxarr
    kyarr=-kyarr
! ------------------------- may 3 add, end ----------------------

    zf=0.0
    qarr=dsqrt(kxarr**2+kyarr**2)
    kappa=qarr
    where (qarr.eq.0)
          invkappa=0.0d0
          farr=0.0d0
          farrz=(1.d0-zf/L)
    elsewhere
         invkappa=1.d0/kappa
         farr=1.d0/(1.d0-dexp(-2.0d0*TWOPI*kappa*L))
         farrz=0.0d0
    end where
    ! Loop over layers in z, divided up among processes.
    !$omp parallel do shared(fbik)
    do k=1,dimz
      zv=zarr(k)
      zf=zarr(k)
      ! FT of field
      if (top_closed) then
        ! Form which avoids overflow for large q values...
        fbik(:,:,k,1)=cmplx(0,1)*kxarr*fbz0*invkappa*farr&
        *( exp(TWOPI*kappa*(zv-2*L))+exp(-TWOPI*kappa*zv) )
        fbik(:,:,k,2)=cmplx(0,1)*kyarr*fbz0*invkappa*farr&
        *( exp(TWOPI*kappa*(zv-2*L))+exp(-TWOPI*kappa*zv) )
        fbik(:,:,k,3)=fbz0*farr&
        *( -exp(TWOPI*kappa*(zv-2*L))+exp(-TWOPI*kappa*zv) ) &
        +farrz*fbz0
        if(Aout) then
          fap(:,:,k,1)=(cmplx(0,1)/TWOPI)*kyarr*invkappa**2*farr*fbz0&
          *(exp(TWOPI*kappa*(zv-2*L))-exp(-TWOPI*kappa*zv))
          fap(:,:,k,2)=(cmplx(0,1)/TWOPI)*kxarr*invkappa**2*farr*fbz0&
          *(-exp(TWOPI*kappa*(zv-2*L))+exp(-TWOPI*kappa*zv))
          fap(:,:,k,3)=0.0d0
        end if
      else
        fbik(:,:,k,1)=cmplx(0,1)*kxarr*fbz0*invkappa*exp(-TWOPI*kappa*zv)
        fbik(:,:,k,2)=cmplx(0,1)*kyarr*fbz0*invkappa*exp(-TWOPI*kappa*zv)
        fbik(:,:,k,3)=fbz0*exp(-TWOPI*kappa*zv)
        if(Aout) then
          fap(:,:,k,1)=(-cmplx(0,1)/TWOPI)*kyarr*invkappa**2*fbz0*exp(-TWOPI*kappa*zv)
          fap(:,:,k,2)=(cmplx(0,1)/TWOPI)*kxarr*invkappa**2*fbz0*exp(-TWOPI*kappa*zv)
          fap(:,:,k,3)=0.0d0
        end if
      end if

      ! Invert FTs
      call dfftw_execute_dft_c2r(plan_bx,fbik(:,:,k,1),bxyz0(:,:,k,1))
      call dfftw_execute_dft_c2r(plan_by,fbik(:,:,k,2),bxyz0(:,:,k,2))
      call dfftw_execute_dft_c2r(plan_bz,fbik(:,:,k,3),bxyz0(:,:,k,3))

      bxyz(1:dimx,1:dimy,k,1)=bxyz0(:,:,k,1)  
      bxyz(1:dimx,1:dimy,k,2)=bxyz0(:,:,k,2)  
      if (top_closed) then
         bxyz(1:dimx,1:dimy,k,3)=bxyz0(:,:,k,3)+bzL
      else 
         bxyz(1:dimx,1:dimy,k,3)=bxyz0(:,:,k,3)
      end if
    enddo
    !$omp end parallel do
    bxyz0=bxyz0*scale_factor
    bxyz=bxyz*scale_factor
    ! Replace exact BCs
    call dfftw_destroy_plan(plan_bx)
    call dfftw_destroy_plan(plan_by)
    call dfftw_destroy_plan(plan_bz)

    if (top_closed) then
       bxyz(1:dimx,1:dimy,1,3)=bz0+bzL
    else 
       bxyz(1:dimx,1:dimy,1,3)=bz0
    end if
    ! Save potential field in b0
    call fill_3d_vector(bxyz)

    if(Aout) then
      !$omp parallel do 
      do k=1,dimz
        call dfftw_execute_dft_c2r(plan_ax,fap(:,:,k,1),axyz(:,:,k,1))
        call dfftw_execute_dft_c2r(plan_ay,fap(:,:,k,2),axyz(:,:,k,2))
        call dfftw_execute_dft_c2r(plan_az,fap(:,:,k,3),axyz(:,:,k,3))
      enddo
      !$omp end parallel do 
    end if 
    ! Deallocate memory
    deallocate(kxarr,kyarr,qarr,qinvarr,farrz)
    deallocate(fbz0,fbik)
    if(Aout) then
      call dfftw_destroy_plan(plan_ax)
      call dfftw_destroy_plan(plan_ay)
      call dfftw_destroy_plan(plan_az)
      axyz=axyz*scale_factor
    end if
  end subroutine fftw3_laplace

  subroutine fftw3_poisson_top_closed()
    integer :: i,j,k,ii
    real(PP) :: kap,fac,zv,emkapz,L
    complex(PP), dimension(3) :: I00,I2,I4
    complex(PP), dimension(2) :: I1,I3,I5
    complex(PP) :: Ifaz_fix
    complex(PP) :: fax,fay,faz,dfaxdz,dfaydz
    real(PP), dimension(:), allocatable :: kxarr, kyarr
    real(PP) :: dkx,dky
    ! Shared parallel objects static
    complex(PP), dimension(dimz) :: arg 
    integer,dimension(8)::plan_fjx,plan_fjy,plan_fjz
    integer,dimension(8)::plan_fbx,plan_fby,plan_fbz
    integer,dimension(8)::plan_fax,plan_fay,plan_faz

    ! Fourier transform current. Note that layers above kmax have fjc(:,:,k,:)=0
    fjc=0.0d0
    fbc=0.0d0
    ! Top of the box
    L=zarr(dimz)
    ! Define arrays u, v with frequencies corresponding to Fourier
    allocate(kxarr(dimx/2+1),kyarr(dimy)) 

! ------------------------- dec 16 add, start ----------------------
    dkx=real(1,PP)/(dimx*dx)
    dky=real(1,PP)/(dimy*dy)

    forall(i=1:dimx/2+1) kxarr(i)=real(i-1,PP)*dkx
    forall(j=1:dimy/2) kyarr(j)=real(j-1,PP)*dky
    forall(j=dimy/2+1:dimy) kyarr(j)=real(j-dimy/2-1,PP)*dky-real(0.5*dimy,PP)*dky

    kxarr=-1*kxarr
    kyarr=-1*kyarr
! ------------------------- dec 16 add, end ----------------------

    call dfftw_plan_dft_c2r_2d(plan_fbx,dimx,dimy,fbc(:,:,1,1),bxyz(:,:,1,1),FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_2d(plan_fby,dimx,dimy,fbc(:,:,1,2),bxyz(:,:,1,2),FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_2d(plan_fbz,dimx,dimy,fbc(:,:,1,3),bxyz(:,:,1,3),FFTW_ESTIMATE)

    call dfftw_plan_dft_r2c_2d(plan_fjx,dimx,dimy,jxyz(:,:,1,1),fjc(:,:,1,1),FFTW_ESTIMATE)
    call dfftw_plan_dft_r2c_2d(plan_fjy,dimx,dimy,jxyz(:,:,1,2),fjc(:,:,1,2),FFTW_ESTIMATE)
    call dfftw_plan_dft_r2c_2d(plan_fjz,dimx,dimy,jxyz(:,:,1,3),fjc(:,:,1,3),FFTW_ESTIMATE)

    !$omp parallel do
    do k=1,dimz
        call dfftw_execute_dft_r2c(plan_fjx,jxyz(:,:,k,1),fjc(:,:,k,1))
        call dfftw_execute_dft_r2c(plan_fjy,jxyz(:,:,k,2),fjc(:,:,k,2))
        call dfftw_execute_dft_r2c(plan_fjz,jxyz(:,:,k,3),fjc(:,:,k,3))
    enddo
    !$omp end parallel do
    call dfftw_destroy_plan(plan_fjx)
    call dfftw_destroy_plan(plan_fjy)
    call dfftw_destroy_plan(plan_fjz)
    ! Apply solution at each grid point
    if(Aout)then
      fap=0.0d0
      axya=0.0d0
      call dfftw_plan_dft_c2r_2d(plan_fax,dimx,dimy,fap(:,:,1,1),axyz(:,:,1,1),FFTW_ESTIMATE)
      call dfftw_plan_dft_c2r_2d(plan_fay,dimx,dimy,fap(:,:,1,2),axyz(:,:,1,2),FFTW_ESTIMATE)
      call dfftw_plan_dft_c2r_2d(plan_faz,dimx,dimy,fap(:,:,1,3),axyz(:,:,1,3),FFTW_ESTIMATE)
    endif

    do i=1,dimx/2+1
      do j=1,dimy
        kap=TWOPI*sqrt(kxarr(i)**2+kyarr(j)**2)
        if (kap.ne.0) then
           fac=1.0d0/(1.0d0-exp(-2.0d0*kap*L))
        else
           fac=0.0d0
        end if
        ! Integral I0
        !$omp parallel do &
        !$omp private(ii,arg)
        do ii=1,3
          arg(1:dimz)=exp(-kap*(zarr(1:dimz)))*fjc(i,j,1:dimz,ii)
          I00(ii)=dz*(0.5d0*(arg(1)+arg(dimz))+sum(arg(2:dimz-1)))
        enddo 
        !$omp end parallel do
        !$omp parallel do &
        !$omp private(k,ii,zv,emkapz,arg,I1,I2,I3,I4,I5) &
        !$omp private(fax,fay,faz,dfaxdz,dfaydz) 
        do k=1,dimz
          zv=zarr(k)
          emkapz=exp(-kap*zv)
          ! Integrals from 0 to z
          I1=0.0d0
          I2=0.0d0
          if (k .ne. 1) then 
            ! Integral I1 
            if (top_closed) then
              do ii=1,2
                arg(1:k)=exp(-kap*(2.0*L-zv+zarr(1:k)))*fjc(i,j,1:k,ii)
                I1(ii)=dz*(0.5d0*(arg(1)+arg(k))+sum(arg(2:k-1))) 
              enddo
            endif
            ! Integral I2
            do ii=1,3
              arg(1:k)=exp(-kap*(zv-zarr(1:k)))*fjc(i,j,1:k,ii)
              I2(ii)=dz*(0.5d0*(arg(1)+arg(k))+sum(arg(2:k-1)))
            enddo
          endif 
          ! Integrals from z to L (or infinity)
          I3=0.0d0
          I4=0.0d0
          if (k .ne. dimz) then
            ! Integral I3 
            do ii=1,2
              arg(k:dimz)=exp(-kap*(2.*L+zv-zarr(k:dimz)))*fjc(i,j,k:dimz,ii)
              I3(ii)=dz*(0.5d0*(arg(k)+arg(dimz))+sum(arg(k+1:dimz-1))) 
            enddo
            ! Integral I4
            do ii=1,3
              arg(k:dimz)=exp(-kap*(zarr(k:dimz)-zv))*fjc(i,j,k:dimz,ii)
              I4(ii)=dz*(0.5d0*(arg(k)+arg(dimz))+sum(arg(k+1:dimz-1))) 
            enddo
          endif 
          ! Integral from 0 to L
          ! Integral I5
          I5=0.0d0
          do ii=1,2
            arg(1:dimz)=exp(-kap*(2.*L-zv-zarr(1:dimz)))*fjc(i,j,1:dimz,ii)
            I5(ii)=dz*(0.5d0*(arg(1)+arg(dimz))+sum(arg(2:dimz-1))) 
          enddo
          ! FTs of x and y components of Ac
          fax=0.5d0*(I1(1)+I2(1)+I3(1)+I4(1)-I5(1)-emkapz*I00(1))
          fay=0.5d0*(I1(2)+I2(2)+I3(2)+I4(2)-I5(2)-emkapz*I00(2))
          if (kap .eq. 0) then
            fax=0.0d0
            fay=0.0d0
          else 
            ! Include Factor of 1./(1-exp(-2.*kap*L))
            fax=fac*fax/kap         
            fay=fac*fay/kap 
          endif
          ! FTs of z derivatives of x and y components
          if (kap .eq. 0) then
           ! It is necessary to treat the zero frequency component 
           ! explicitly - it is an integral from z to L, which 
           ! ensures the Poisson equation is met when kappa=0
            dfaxdz=0.0d0
            dfaydz=0.0d0 
            if (k .ne. dimz) then
               arg(k:dimz)=fjc(i,j,k:dimz,1)
               dfaxdz=dz*(0.5*(arg(k)+arg(dimz))+sum(arg(k+1:dimz-1))) 
               arg(k:dimz)=fjc(i,j,k:dimz,2)
               dfaydz=dz*(0.5*(arg(k)+arg(dimz))+sum(arg(k+1:dimz-1))) 
            endif 
          else 
            dfaxdz=0.5d0*(I1(1)-I2(1)-I3(1)+I4(1)-I5(1)+emkapz*I00(1))
            dfaydz=0.5d0*(I1(2)-I2(2)-I3(2)+I4(2)-I5(2)+emkapz*I00(2))
            dfaxdz=fac*dfaxdz
            dfaydz=fac*dfaydz
          endif       
          ! FT of z component
          if (kap .eq. 0) then
            faz=0.0d0
          else
            faz=fac*( (1.0d0+exp(-2.0*kap*L))*emkapz+& 
                2.0d0*exp(-kap*(2.0*L-zv)) )*I00(3)+I2(3)+I4(3)
! --------------------------------------------------------------------
! this additional term for faz corrects the coulomb gauge
            arg(1:dimz)=exp(kap*(zv - 2.0d0*L + zarr(1:dimz))) &
                      + exp(kap*(zarr(1:dimz) - zv - 2.0d0*L)) &
                      - exp(kap*(zv - 2.0d0*L - zarr(1:dimz))) &
                      - exp(-1.0d0*kap*(zv + 2.0d0*L + zarr(1:dimz)))
            arg(1:dimz)=arg(1:dimz)*fjc(i,j,1:dimz,3)
            Ifaz_fix=fac*dz*(0.5d0*(arg(1)+arg(dimz))+sum(arg(2:dimz-1)))
            faz = faz + Ifaz_fix
! --------------------------------------------------------------------
            faz=0.5d0*faz/kap
          endif
          ! Construct FTs of components of Bc
          fbc(i,j,k,1)=-TWOPI*cmplx(0,1)*kyarr(j)*faz-dfaydz
          fbc(i,j,k,2)=dfaxdz+TWOPI*cmplx(0,1)*kxarr(i)*faz
          fbc(i,j,k,3)=TWOPI*cmplx(0,1)*(-kxarr(i)*fay+kyarr(j)*fax)
          if(Aout) then
            fap(i,j,k,1)=fax
            fap(i,j,k,2)=fay
            fap(i,j,k,3)=faz
          endif
        enddo
        !$omp end parallel do
      enddo
    enddo
    fbc=fbc*scale_factor

    !$omp parallel do
    do k=1,dimz
        call dfftw_execute_dft_c2r(plan_fbx,fbc(:,:,k,1),bxyzj(:,:,k,1))
        call dfftw_execute_dft_c2r(plan_fby,fbc(:,:,k,2),bxyzj(:,:,k,2))
        call dfftw_execute_dft_c2r(plan_fbz,fbc(:,:,k,3),bxyzj(:,:,k,3))
    enddo
    !$omp end parallel do
    call dfftw_destroy_plan(plan_fbx)
    call dfftw_destroy_plan(plan_fby)
    call dfftw_destroy_plan(plan_fbz)

    if(Aout) then
      !$omp parallel do 
      do k=1,dimz 
        call dfftw_execute_dft_c2r(plan_fax,fap(:,:,k,1),axyz(:,:,k,1))
        call dfftw_execute_dft_c2r(plan_fay,fap(:,:,k,2),axyz(:,:,k,2))
        call dfftw_execute_dft_c2r(plan_faz,fap(:,:,k,3),axyz(:,:,k,3))
      enddo
      !$omp end parallel do  
      call dfftw_destroy_plan(plan_fax)
      call dfftw_destroy_plan(plan_fay)
      call dfftw_destroy_plan(plan_faz)
      axyz=axyz*scale_factor
    endif
  end subroutine fftw3_poisson_top_closed

  subroutine fftw3_poisson_top_open()
    integer :: i,j,k,ii
    real(PP) :: kap,fac,zv,emkapz,L
    complex(PP), dimension(3) :: I00,I2,I4
    complex(PP), dimension(2) :: I1,I3,I5
    complex(PP) :: fax,fay,faz,dfaxdz,dfaydz
    real(PP), dimension(:), allocatable :: kxarr, kyarr
    real(PP) :: dkx,dky
    ! Shared parallel objects static
    complex(PP), dimension(dimz) :: arg 
    integer,dimension(8)::plan_fjx,plan_fjy,plan_fjz
    integer,dimension(8)::plan_fbx,plan_fby,plan_fbz
    integer,dimension(8)::plan_fax,plan_fay,plan_faz
    real(PP), dimension(:), allocatable :: zarr1

    ! print *,'Calculating Fourier solution to Poisson equation...'
    kmax=dimz
    ! Fourier transform current. Note that layers above kmax have fjc(:,:,k,:)=0
    fjc=0.0d0
    fbc=0.0d0
    ! Top of the box
    L=zarr(dimz)
    ! Define arrays u, v with frequencies corresponding to Fourier
    allocate(kxarr(dimx/2+1),kyarr(dimy)) 
    allocate(zarr1(dimz))
    zarr1=zarr

! ------------------------- dec 16 add, start ----------------------
    dkx=real(1,PP)/(dimx*dx)
    dky=real(1,PP)/(dimy*dy)

    forall(i=1:dimx/2+1) kxarr(i)=real(i-1,PP)*dkx
    forall(j=1:dimy/2) kyarr(j)=real(j-1,PP)*dky
    forall(j=dimy/2+1:dimy) kyarr(j)=real(j-dimy/2-1,PP)*dky-real(0.5*dimy,PP)*dky

    kxarr=-1*kxarr
    kyarr=-1*kyarr
! ------------------------- dec 16 add, end ----------------------

    call dfftw_plan_dft_c2r_2d(plan_fbx,dimx,dimy,fbc(:,:,1,1),bxyz(:,:,1,1),FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_2d(plan_fby,dimx,dimy,fbc(:,:,1,2),bxyz(:,:,1,2),FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_2d(plan_fbz,dimx,dimy,fbc(:,:,1,3),bxyz(:,:,1,3),FFTW_ESTIMATE)

    call dfftw_plan_dft_r2c_2d(plan_fjx,dimx,dimy,jxyz(:,:,1,1),fjc(:,:,1,1),FFTW_ESTIMATE)
    call dfftw_plan_dft_r2c_2d(plan_fjy,dimx,dimy,jxyz(:,:,1,2),fjc(:,:,1,2),FFTW_ESTIMATE)
    call dfftw_plan_dft_r2c_2d(plan_fjz,dimx,dimy,jxyz(:,:,1,3),fjc(:,:,1,3),FFTW_ESTIMATE)

    !$omp parallel do
    do k=1,dimz
        call dfftw_execute_dft_r2c(plan_fjx,jxyz(:,:,k,1),fjc(:,:,k,1))
        call dfftw_execute_dft_r2c(plan_fjy,jxyz(:,:,k,2),fjc(:,:,k,2))
        call dfftw_execute_dft_r2c(plan_fjz,jxyz(:,:,k,3),fjc(:,:,k,3))
    enddo
    !$omp end parallel do
    call dfftw_destroy_plan(plan_fjx)
    call dfftw_destroy_plan(plan_fjy)
    call dfftw_destroy_plan(plan_fjz)

    ! Apply solution at each grid point
    if(Aout)then
      fap=0.0d0
      axyz=0.0d0
      call dfftw_plan_dft_c2r_2d(plan_fax,dimx,dimy,fap(:,:,1,1),axyz(:,:,1,1),FFTW_ESTIMATE)
      call dfftw_plan_dft_c2r_2d(plan_fay,dimx,dimy,fap(:,:,1,2),axyz(:,:,1,2),FFTW_ESTIMATE)
      call dfftw_plan_dft_c2r_2d(plan_faz,dimx,dimy,fap(:,:,1,3),axyz(:,:,1,3),FFTW_ESTIMATE)
    endif

    do i=1,dimx/2+1
      do j=1,dimy
        kap=TWOPI*sqrt(kxarr(i)**2+kyarr(j)**2)
        if (kap.ne.0) then
           fac=1.0d0/(1.0d0-exp(-2.0d0*kap*L))
        else
           fac=0.0d0
        end if
        ! Integral I0
        !$omp parallel do &
        !$omp private(ii,arg)
        do ii=1,3
          arg(1:dimz)=exp(-kap*(zarr1(1:dimz)))*fjc(i,j,1:dimz,ii)
          I00(ii)=dz*(0.5d0*(arg(1)+arg(dimz))+sum(arg(2:dimz-1)))
        enddo 
        !$omp end parallel do
        !$omp parallel do &
        !$omp private(k,ii,zv,emkapz,arg,I1,I2,I3,I4,I5) &
        !$omp private(fax,fay,faz,dfaxdz,dfaydz) 
        do k=1,dimz
          zv=zarr1(k)
          emkapz=exp(-kap*zv)
          ! Integrals from 0 to z
          I1=0.0d0
          I2=0.0d0
          if (k .ne. 1) then
            ! Integral I2
            do ii=1,3
              arg(1:k)=exp(-kap*(zv-zarr1(1:k)))*fjc(i,j,1:k,ii)
              I2(ii)=dz*(0.5d0*(arg(1)+arg(k))+sum(arg(2:k-1)))
            enddo
          endif 
          ! Integrals from z to L (or infinity)
          I3=0.0d0
          I4=0.0d0
          if (k .ne. dimz) then
            ! Integral I4
            do ii=1,3
              arg(k:dimz)=exp(-kap*(zarr1(k:dimz)-zv))*fjc(i,j,k:dimz,ii)
              I4(ii)=dz*(0.5d0*(arg(k)+arg(dimz))+sum(arg(k+1:dimz-1))) 
            enddo
          endif 
          ! Integral from 0 to L (or infinity)
          ! Integral I5
          I5=0.0d0
          ! FTs of x and y components of Ac
          fax=0.5d0*(I1(1)+I2(1)+I3(1)+I4(1)-I5(1)-emkapz*I00(1))
          fay=0.5d0*(I1(2)+I2(2)+I3(2)+I4(2)-I5(2)-emkapz*I00(2))
          if (kap .eq. 0) then
            fax=0.0d0
            fay=0.0d0
          else 
            fax=fax/kap
            fay=fay/kap       
          endif
          ! FTs of z derivatives of x and y components
          if (kap .eq. 0) then
           ! It is necessary to treat the zero frequency component 
           ! explicitly - it is an integral from z to L, which 
           ! ensures the Poisson equation is met when kappa=0
            dfaxdz=0.0d0
            dfaydz=0.0d0 
            if (k .ne. dimz) then
               arg(k:dimz)=fjc(i,j,k:dimz,1)
               dfaxdz=dz*(0.5*(arg(k)+arg(dimz))+sum(arg(k+1:dimz-1))) 
               arg(k:dimz)=fjc(i,j,k:dimz,2)
               dfaydz=dz*(0.5*(arg(k)+arg(dimz))+sum(arg(k+1:dimz-1))) 
            endif 
          else 
            dfaxdz=0.5d0*(I1(1)-I2(1)-I3(1)+I4(1)-I5(1)+emkapz*I00(1))
            dfaydz=0.5d0*(I1(2)-I2(2)-I3(2)+I4(2)-I5(2)+emkapz*I00(2))
          endif       
          ! FT of z component
          if (kap .eq. 0) then
            faz=0.0d0
          else
            faz=I2(3)+I4(3)+emkapz*I00(3)
            faz=0.5d0*faz/kap
          endif
          ! Construct FTs of components of Bc
          fbc(i,j,k,1)=-TWOPI*cmplx(0,1)*kyarr(j)*faz-dfaydz
          fbc(i,j,k,2)=dfaxdz+TWOPI*cmplx(0,1)*kxarr(i)*faz
          fbc(i,j,k,3)=TWOPI*cmplx(0,1)*(-kxarr(i)*fay+kyarr(j)*fax)
          if(Aout) then
            fap(i,j,k,1)=fax
            fap(i,j,k,2)=fay
            fap(i,j,k,3)=faz
          endif
        enddo
        !$omp end parallel do
      enddo
    enddo
    fbc=fbc*scale_factor

    !$omp parallel do
    do k=1,dimz
        call dfftw_execute_dft_c2r(plan_fbx,fbc(:,:,k,1),bxyzj(:,:,k,1))
        call dfftw_execute_dft_c2r(plan_fby,fbc(:,:,k,2),bxyzj(:,:,k,2))
        call dfftw_execute_dft_c2r(plan_fbz,fbc(:,:,k,3),bxyzj(:,:,k,3))
    enddo
    !$omp end parallel do
    call dfftw_destroy_plan(plan_fbx)
    call dfftw_destroy_plan(plan_fby)
    call dfftw_destroy_plan(plan_fbz)

    if(Aout) then
      !$omp parallel do 
      do k=1,dimz 
        call dfftw_execute_dft_c2r(plan_fax,fap(:,:,k,1),axyz(:,:,k,1))
        call dfftw_execute_dft_c2r(plan_fay,fap(:,:,k,2),axyz(:,:,k,2))
        call dfftw_execute_dft_c2r(plan_faz,fap(:,:,k,3),axyz(:,:,k,3))
      enddo
      !$omp end parallel do  
      call dfftw_destroy_plan(plan_fax)
      call dfftw_destroy_plan(plan_fay)
      call dfftw_destroy_plan(plan_faz)
      axyz=axyz*scale_factor
    endif
  end subroutine fftw3_poisson_top_open

end module mod_fftw