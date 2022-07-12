module mod_operator
use mod_param
implicit none
contains

  function cross(a, b)
    real(PP),dimension(3) :: cross
    real(PP),dimension(3),intent(in) :: a, b
    
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross

  subroutine periodic_posi(posiin,posiout)
    real(PP)::posiin(3),posiout(3)
    posiout=posiin
    if (posiin(1) .le. 0) posiout(1)=dimx+posiin(1)
    if (posiin(1) .ge. dimx+1) posiout(1)=posiin(1)-dimx
    if (posiin(2) .le. 0) posiout(2)=dimy+posiin(2)
    if (posiin(2) .ge. dimy+1) posiout(2)=posiin(2)-dimy
  end subroutine periodic_posi

! ------- function check whether the point outside the computational domain
  function outside_boundary(posi)
    logical:: outside_boundary
    real(PP),intent(in) :: posi(3)
    
    if(Periodic) then
      if(posi(3) .ge. 1 .and. posi(3) .le. dimz ) then
        outside_boundary = .false.
      else
        outside_boundary = .true.
      end if
    else
      if(posi(1) .ge. 1 .and. posi(1) .le. dimx &
      & .and. posi(2) .ge. 1 .and. posi(2) .le. dimy &
      & .and. posi(3) .ge. 1 .and. posi(3) .le. dimz ) then
        outside_boundary = .false.
      else
        outside_boundary = .true.
      end if
    endif

  end function outside_boundary
! ---------------------------------------------------------------
! the first dimension of xv indicates the x,y,z components of B
! the location in the second dimension of xv 
! 
!               g       g
!              /       /
!             /       /
!    g- - - -8- - - -7- - - -g
!           /|      /|   
!          / |     / |  
! g- - - -5- - - -6  |
!   g- - -|- 4 - -|- 3- - - -g
!         | /     | /
!         |/      |/     
! g- - - -1- - - -2- - - -g
!         |       |
!         |       |
!         |       |
!         g- - - -g
! 
! 
! ------ calculate the Bxyz on 8 vertex -----------
  subroutine corner (nnn,xv)
    integer,intent(in)::nnn(3)
    integer::n1,n2,n3
    integer::n1p,n2p,n3p
    real(PP),intent(out)::xv(3,8)
    n1 = nnn(1)
    n2 = nnn(2)
    n3 = nnn(3)
    
    n1p = nnn(1)+1
    n2p = nnn(2)+1
    n3p = nnn(3)+1
    
    if(Periodic) then
      if(n1 .lt. 1) n1=1
      if(n2 .lt. 1) n2=1
      if(n1p .gt. dimx) n1p=n1p-dimx
      if(n2p .gt. dimy) n2p=n2p-dimy
    else
      if(n1 .lt. 1) n1=1
      if(n2 .lt. 1) n2=1
      if(n1p .gt. dimx) n1p=dimx
      if(n2p .gt. dimy) n2p=dimy
    endif    

    if(n3 .lt. 1) n3=1
    if(n3p .gt. dimz) n3p=dimz

    xv(1,1) = bxyz(n1,n2,n3,1)
    xv(1,2) = bxyz(n1p,n2,n3,1)
    xv(1,3) = bxyz(n1p,n2p,n3,1)
    xv(1,4) = bxyz(n1,n2p,n3,1)
    xv(1,5) = bxyz(n1,n2,n3p,1)
    xv(1,6) = bxyz(n1p,n2,n3p,1)
    xv(1,7) = bxyz(n1p,n2p,n3p,1)
    xv(1,8) = bxyz(n1,n2p,n3p,1)
    
    xv(2,1) = bxyz(n1,n2,n3,2)
    xv(2,2) = bxyz(n1p,n2,n3,2)
    xv(2,3) = bxyz(n1p,n2p,n3,2)
    xv(2,4) = bxyz(n1,n2p,n3,2)
    xv(2,5) = bxyz(n1,n2,n3p,2)
    xv(2,6) = bxyz(n1p,n2,n3p,2)
    xv(2,7) = bxyz(n1p,n2p,n3p,2)
    xv(2,8) = bxyz(n1,n2p,n3p,2)
    
    xv(3,1) = bxyz(n1,n2,n3,3)
    xv(3,2) = bxyz(n1p,n2,n3,3)
    xv(3,3) = bxyz(n1p,n2p,n3,3)
    xv(3,4) = bxyz(n1,n2p,n3,3)
    xv(3,5) = bxyz(n1,n2,n3p,3)
    xv(3,6) = bxyz(n1p,n2,n3p,3)
    xv(3,7) = bxyz(n1p,n2p,n3p,3)
    xv(3,8) = bxyz(n1,n2p,n3p,3)
  end subroutine corner
! --------- get the tri-linear interpolation ----------
    ! the location of each index 
    !     8- - - -7
    !    /|      /|   
    !   / |     / |  
    !  5- - - -6  |
    !  |  4 - -|- 3   
    !  | /     | /
    !  |/      |/     
    !  1- - - -2
    ! 
    ! o indicates the location of dxyz in the unit square
    ! (0,1)---------(1,1)
    !     |   |     |
    !     |   |     |
    !     |---o-----|
    !     |   |     |
    ! (0,0)---------(1,0)
    ! the left bottom one is dxyz1
    ! the right up one is dxyz2=1-dxyz
  subroutine xitp (interp,xv,dxyz1)
    real(PP),intent(out)::interp(3)
    real(PP)::dxyz2(3),weight(8)
    real(PP)::xv(3,8)
    real(PP),intent(in)::dxyz1(3)

    dxyz2 = real(1,PP) - dxyz1
    
    weight(1) = dxyz2(1)*dxyz2(2)*dxyz2(3)
    weight(2) = dxyz1(1)*dxyz2(2)*dxyz2(3)
    weight(3) = dxyz1(1)*dxyz1(2)*dxyz2(3)
    weight(4) = dxyz2(1)*dxyz1(2)*dxyz2(3)
    
    weight(5) = dxyz2(1)*dxyz2(2)*dxyz1(3)
    weight(6) = dxyz1(1)*dxyz2(2)*dxyz1(3)
    weight(7) = dxyz1(1)*dxyz1(2)*dxyz1(3)
    weight(8) = dxyz2(1)*dxyz1(2)*dxyz1(3)
    
    interp(1) = dot_product(xv(1,:),weight)
    interp(2) = dot_product(xv(2,:),weight)
    interp(3) = dot_product(xv(3,:),weight)

  end subroutine xitp
! ------ give the equation of ODE ------
  ! subroutine diffLine(Posi,Tangent)
  !   real(PP),intent(in)::Posi(neqn)
  !   real(PP),intent(out)::Tangent(neqn+1)
  !   real(PP)::Posi_tmp(neqn)
  !   real(PP)::dxyz(3)
  !   real(PP)::xv(3,8)
  !   real(PP)::bvec(3)
  !   real(PP)::binter
  !   real(PP)::eps
  !   integer :: nnn(3)
  !   Posi_tmp=Posi
  !   eps = epsilon(real(1,PP))
  !   dxyz = Posi_tmp - floor(Posi_tmp)
  !   nnn = floor(Posi_tmp)
  !   call corner(nnn,xv)
  !   call xitp(bvec,xv,dxyz)
  !   binter = sqrt(dot_product(bvec,bvec))    
  !   if (binter .lt. eps) then
  !      Tangent = real(0,PP)*(/1,1,1,1/)
  !   else
  !      Tangent(1:3) = bvec/binter
  !      Tangent(4)   = binter
  !   end if
  ! end subroutine diffLine

! ----------- R.H.S. of the ODE ----------
  ! subroutine rhs(t,y,yp)
  !   real(PP)::t
  !   real(PP)::y(neqn)
  !   real(PP)::yp(neqn)
  !   real(PP)::TangentB(4)
    
  !   t=t*real(1,PP)
  !   call diffLine(y,TangentB)
  !   yp = TangentB(1:3)
  ! end subroutine rhs

  subroutine diffLine(t,Posi,Tangent)
    real(PP) :: t
    real(PP),intent(in)::Posi(neqn)
    real(PP),intent(out)::Tangent(neqn)
    real(PP)::Posi_tmp(neqn)
    real(PP)::dxyz(3)
    real(PP)::xv(3,8)
    real(PP)::bvec(3)
    real(PP)::binter
    real(PP)::eps
    integer :: nnn(3)
    t=t*1.0d0
    Posi_tmp=Posi
    eps = epsilon(real(1,PP))
    dxyz = Posi_tmp - floor(Posi_tmp)
    nnn = floor(Posi_tmp)
    call corner(nnn,xv)
    call xitp(bvec,xv,dxyz)
    binter = sqrt(dot_product(bvec,bvec))    
    if (binter .lt. eps) then
       Tangent = real(0,PP)*(/1,1,1/)/dsqrt(3.0d0)
    else
       Tangent = bvec/binter
    end if
  end subroutine diffLine
! ---- the ODE solver runge-kutta 2 order method -----
  subroutine rk2 (f, neqn, t, t_out, y, flag)
    integer ( kind = i4 ) neqn
    external f
    real(PP) :: t,t_out,dt
    real(PP) :: y(neqn)
    real(PP) :: y_tmp(neqn)
    real(PP) :: k1(neqn)
    real(PP) :: k2(neqn)
    real(PP) :: yp(neqn)
    real(PP) :: eps
    logical :: flag
    integer :: iterNum
    
    eps = epsilon(real(1,PP))
    
    iterNum = 0
    flag = .true.
    dt = t_out - t
    y_tmp = y
    call f(t,y,yp)
    k1 = dt*yp
    
    call f(t+0.6667*dt,y + 0.6667*k1,yp)
    k2 = dt*yp
        
    y_tmp =  y + (0.25*k1+0.75*k2)
    
    do while(outside_boundary(y_tmp) .and. iterNum .le. iterMax .and. dt .gt. min_step)
       
       iterNum = iterNum + 1
       dt = dt * real(0.5,PP)
       
       call f(t,y,yp)
       k1 = dt*yp
    
       call f(t+0.6667*dt,y + 0.6667*k1,yp)
       k2 = dt*yp

       y_tmp =  y + (0.25*k1+0.75*k2)
    
    end do
    
    y = y_tmp
        
    if(outside_boundary(y_tmp)) then
      flag=.false.
    endif

    t_out = t + dt
  end subroutine rk2

! ---- the ODE solver runge-kutta 4 order method -----
  subroutine rk4 (f, neqn, t, t_out, y, flag)
    integer ( kind = i4 ) neqn
    external f
    real(PP) :: t,t_out,dt
    real(PP) :: y(neqn)
    real(PP) :: y_tmp(neqn)
    real(PP) :: k1(neqn)
    real(PP) :: k2(neqn)
    real(PP) :: k3(neqn)
    real(PP) :: k4(neqn)
    real(PP) :: yp(neqn)
    real(PP) :: eps
    logical :: flag
    integer :: iterNum
    
    eps = epsilon(real(1,PP))
    
    iterNum = 0
    flag = .true.
    dt = t_out - t
    y_tmp = y
    
    call f(t,y,yp)
    k1 = dt*yp
    
    call f(t+0.5*dt,y + 0.5*k1,yp)
    k2 = dt*yp
    
    call f(t+0.5*dt,y+0.5*k2,yp)
    k3 = dt*yp
    
    call f(t+dt,y+k3,yp)
    k4 = dt*yp
    
    y_tmp =  y &
        & +( k1 &
        & +  k2*real(2,PP) &
        & +  k3*real(2,PP) &
        & +  k4 )/real(6,PP)
    
    do while(outside_boundary(y_tmp) .and. iterNum .le. iterMax .and. dt .gt. min_step)
       
       iterNum = iterNum + 1
       dt = dt * real(0.5,PP)
       
       call f(t,y,yp)
       k1 = dt*yp
    
       call f(t+0.5*dt,y + 0.5*k1,yp)
       k2 = dt*yp
    
       call f(t+0.5*dt,y+0.5*k2,yp)
       k3 = dt*yp
    
       call f(t+dt,y+k3,yp)
       k4 = dt*yp
    
       y_tmp =  y &
           & +( k1 &
           & +  k2*real(2,PP) &
           & +  k3*real(2,PP) &
           & +  k4 )/real(6,PP)
    end do
    
    y = y_tmp

    if(outside_boundary(y_tmp)) then
      flag=.false.
    endif
    
    t_out = t + dt
  end subroutine rk4
! ----------- R.H.S. of the ODE ----------
  ! subroutine rhs(t,y,yp)
  !   real(PP)::t
  !   real(PP)::y(neqn)
  !   real(PP)::yp(neqn)
  !   real(PP)::TangentB(4)
    
  !   t=t*real(1,PP)
  !   call diffLine(y,TangentB)
  !   yp = TangentB(1:3)
  ! end subroutine rhs

  subroutine fill_3d_scalar(var)
    real(PP), dimension(0:dimx+1,0:dimy+1,0:dimz+1)::var
    if (Periodic) then
      var(0,:,:)=var(dimx,:,:)
      var(dimx+1,:,:)=var(1,:,:)
      var(:,0,:)=var(:,dimy,:)
      var(:,dimy+1,:)=var(:,1,:)
    else
      var(0,:,:)=var(1,:,:)
      var(dimx+1,:,:)=var(dimx,:,:)
      var(:,0,:)=var(:,1,:)
      var(:,dimy+1,:)=var(:,dimy,:)
    endif      
    var(:,:,0)=var(:,:,1)
    var(:,:,dimz+1)=var(:,:,dimz)
  end subroutine fill_3d_scalar

  subroutine fill_3d_vector(var)
    real(PP),dimension(0:dimx+1,0:dimy+1,0:dimz+1,3)::var
    call fill_3d_scalar(var(:,:,:,1))
    call fill_3d_scalar(var(:,:,:,2))
    call fill_3d_scalar(var(:,:,:,3))
  end subroutine fill_3d_vector

  subroutine fill_edge_vertex_scalar(var)
    real(PP), dimension(0:dimx+1,0:dimy+1)::var
    if (Periodic) then
      var(0,:)=var(dimx,:)
      var(dimx+1,:)=var(1,:)
      var(:,0)=var(:,dimy)
      var(:,dimy+1)=var(:,1)
    else
      var(0,:)=var(1,:)
      var(dimx+1,:)=var(dimx,:)
      var(:,0)=var(:,1)
      var(:,dimy+1)=var(:,dimy)
    endif      
  end subroutine fill_edge_vertex_scalar

  subroutine fill_face_edge_vertex_scalar(var)
    real(PP), dimension(0:dimx+1,0:dimy+1,0:dimz+1)::var

    var(0,1:dimy,1:dimz) = var(1,1:dimy,1:dimz)
    var(dimx+1,1:dimy,1:dimz) = var(dimx,1:dimy,1:dimz)      
    var(1:dimx,0,1:dimz) = var(1:dimx,1,1:dimz)
    var(1:dimx,dimy+1,1:dimz) = var(1:dimx,dimy,1:dimz)
    var(1:dimx,1:dimy,0) = var(1:dimx,1:dimy,1)
    var(1:dimx,1:dimy,dimz+1) = var(1:dimx,1:dimy,dimz)
    ! lines in ghost shell
    var(0,0,1:dimz) = var(1,1,1:dimz)
    var(dimx+1,0,1:dimz) = var(dimx,1,1:dimz)
    var(0,dimy+1,1:dimz) = var(1,dimy,1:dimz)
    var(dimx+1,dimy+1,1:dimz) = var(dimx,dimy,1:dimz)
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    var(1:dimx,0,0) = var(1:dimx,1,1)
    var(1:dimx,dimy+1,0) = var(1:dimx,dimy,1)
    var(1:dimx,dimy+1,dimz+1) = var(1:dimx,dimy,dimz)
    var(1:dimx,0,dimz+1) = var(1:dimx,1,dimz)
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    var(0,1:dimy,0) = var(1,1:dimy,1)
    var(dimx+1,1:dimy,0) = var(dimx,1:dimy,1)
    var(0,1:dimy,dimz+1) = var(1,1:dimy,dimz)
    var(dimx+1,1:dimy,dimz+1) = var(dimx,1:dimy,dimz)
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    ! points in ghost shell
    var(0,0,0) = var(1,1,1)
    var(dimx+1,0,0) = var(dimx,1,1)
    var(0,dimy+1,0) = var(1,dimy,1)
    var(0,0,dimz+1) = var(1,1,dimz)
    var(dimx+1,dimy+1,0) = var(dimx,dimy,1)
    var(dimx+1,0,dimz+1) = var(dimx,1,dimz)
    var(0,dimy+1,dimz+1) = var(1,dimy,dimz)
    var(dimx+1,dimy+1,dimz+1) = var(dimx,dimy,dimz)
  end subroutine fill_face_edge_vertex_scalar

  subroutine fill_face_edge_vertex_vector(var)
    real(PP), dimension(0:dimx+1,0:dimy+1,0:dimz+1,3)::var
    integer :: ix
    do ix=1,3
      call fill_face_edge_vertex_scalar(var(:,:,:,ix))
    enddo
  end subroutine fill_face_edge_vertex_vector

  subroutine alpha_bottom(posi,alpha_tmp)
    real(PP),intent(in) ::posi(3)
    real(PP),intent(out)::alpha_tmp
    real(PP)::weight(4),MaskCom(4)
    real(PP)::xindex,yindex,zindex
    integer::i,j,k
    alpha_tmp = real(0,PP)
    xindex = posi(1)
    yindex = posi(2)
    zindex = posi(3)
    i=floor(xindex)
    j=floor(yindex)
    k=floor(zindex)
    ! calculate the alpha component
    MaskCom(1) = alpha0(i,j)
    MaskCom(2) = alpha0(i+1,j)
    MaskCom(3) = alpha0(i,j+1)
    MaskCom(4) = alpha0(i+1,j+1)
    ! calculate the alpha value
    weight(1) = (real(i - xindex+1,PP) )*(real(j-yindex+1,PP) )
    weight(2) = real(xindex - i,PP)*(real(j-yindex+1,PP) )
    weight(3) = (real(i-xindex+1,PP) )*real(yindex - j,PP)
    weight(4) = real(xindex - i,PP)*real(yindex - j,PP)
    alpha_tmp = dot_product(weight,MaskCom)
  end subroutine alpha_bottom

  subroutine sigma_bottom(posi,alpha_tmp)
    real(PP),intent(in) ::posi(3)
    real(PP),intent(out)::alpha_tmp
    real(PP)::weight(4),MaskCom(4)
    real(PP)::xindex,yindex,zindex
    integer::i,j,k
    alpha_tmp = real(0,PP)
    xindex = posi(1)
    yindex = posi(2)
    zindex = posi(3)
    i=floor(xindex)
    j=floor(yindex)
    k=floor(zindex)

    ! calculate the alpha component
    MaskCom(1) = sig0(i,j)
    MaskCom(2) = sig0(i+1,j)
    MaskCom(3) = sig0(i,j+1)
    MaskCom(4) = sig0(i+1,j+1)
    ! calculate the alpha value
    weight(1) = (real(i - xindex+1,PP) )*(real(j-yindex+1,PP) )
    weight(2) = real(xindex - i,PP)*(real(j-yindex+1,PP) )
    weight(3) = (real(i-xindex+1,PP) )*real(yindex - j,PP)
    weight(4) = real(xindex - i,PP)*real(yindex - j,PP)
    alpha_tmp = dot_product(weight,MaskCom)
  end subroutine sigma_bottom

  subroutine calc_angle(var)
    real(PP)::var
    real(PP)::sint(dimx,dimy,dimz)
    real(PP)::jabs2(dimx,dimy,dimz)
    real(PP)::btmp(dimx,dimy,dimz,3)
    btmp=bxyz(1:dimx,1:dimy,1:dimz,:)

    jabs2=sum(jxyz**2,dim=4)
    sint=sqrt(jabs2-sum(jxyz*btmp,dim=4)/sum(btmp**2,dim=4))
    var=sum(sint(1:dimx,1:dimy,1:dimz))/sum(jabs2(1:dimx,1:dimy,1:dimz))
  end subroutine calc_angle

  subroutine calc_energy(var)
    real(PP)::var
    var=sum(bxyz(1:dimx,1:dimy,1:dimz,3)**2)
    var=var
  end subroutine calc_energy

  subroutine combine_alpha
    print *,'---'
    print *,'Self-consistency cycle complete'
    print *,'Constructing new alpha boundary values'
    print *,'---'
    alpha0=0.0
    ! sig0=0.0 ! Changed MSW Aug 2017
    if (alpha_error) then
      ! Put maximum uncertainties at points with sig = 0
      where (sig0p .eq. 0)
        sig0p=max_sig0
      end where
      where (sig0n .eq. 0)
        sig0n=max_sig0
      end where
      alpha0(1:dimx,1:dimy)=(alpha_p(:,:,1)/sig0p**2+alpha_n(:,:,1)/sig0n**2)/ &
             (1./sig0p**2+1./sig0n**2) 
      ! Change MSW Aug 2017
      ! sig00=sqrt(2.)/sqrt(1./sig0p**2+1./sig0n**2)
      ! Reverted MSW Aug 2017
      sig00=1./sqrt(1./sig0p**2+1./sig0n**2)
      print *,'max(alpha0)=',maxval(alpha0)
      print *,'min(alpha0)=',minval(alpha0)
      print *,'max_sig0=',max_sig0 ! Changed MSW Aug 2017
      print *,'max(sig00)=',maxval(sig00) ! Changed MSW Aug 2017
      print *,'min(sig00)=',minval(sig00) ! Changed MSW Aug 2017
    else
       alpha0(1:dimx,1:dimy)=0.5*(alpha_p(:,:,1)+alpha_n(:,:,1))
    endif
    call fill_edge_vertex_scalar(alpha0)
  end subroutine combine_alpha

  subroutine std(array,std_n0)
    real(PP), dimension(:,:) :: array
    real(PP) :: avg_n0,std_n0
  
    avg_n0 = sum(array, MASK= array /= 0.0)/count(array /= 0.0)
    std_n0 = sqrt(sum((array-avg_n0)**2, MASK= array /= 0.0) &
             /count(array /= 0.0))
  end subroutine std

  subroutine curvature(array,carray)
    integer :: i,j,nx,ny
    real(PP), dimension(:,:) :: array
    real(PP), dimension(:,:) :: carray
  
    nx=size(array,1)
    ny=size(array,2)
  
    carray=0.0
    !$omp parallel do private (j)
    do i=2,nx-1
      do j=2,ny-1
          carray(i,j)=array(i,j)-0.25*(array(i+1,j)+array(i-1,j) &
                                      +array(i,j+1)+array(i,j-1))
      enddo
    enddo
    !$omp end parallel do
  end subroutine curvature

  real(PP) function xderiv(array,i,j,k)
    real(PP), dimension(:,:,:) :: array
    integer :: i,j,k
    real(PP) :: idx
    idx=1.0d0/dx
    if ((i /= 1) .and. (i /= dimx)) then
      xderiv=0.5d0*idx*(array(i+1,j,k)-array(i-1,j,k))
    else if (i .eq. 1) then
      xderiv=(-3.0d0*array(1,j,k)+4.0d0*array(2,j,k)-array(3,j,k))/(2.0d0*dx)
    else if (i .eq. dimx) then
      xderiv=(3.0d0*array(dimx,j,k)-4.0d0*array(dimx-1,j,k)+array(dimx-2,j,k))/(2.0d0*dx)
    end if
  end function xderiv

  !Derivative in y
  real(PP) function yderiv(array,i,j,k)
    real(PP), dimension(:,:,:) :: array
    integer :: i,j,k
    real(PP) :: idy
    idy=1.0d0/dy
    if ((j /= 1) .and. (j /= dimy)) then
      yderiv=0.5d0*idy*(array(i,j+1,k)-array(i,j-1,k))
    else if (j .eq. 1) then
      yderiv=(-3.0d0*array(i,1,k)+4.0d0*array(i,2,k)-array(i,3,k))/(2.0d0*dy)
    else if (j .eq. dimy) then
      yderiv=(3.0d0*array(i,dimy,k)-4.0d0*array(i,dimy-1,k)+array(i,dimy-2,k))/(2.0d0*dy)
    end if
  end function yderiv

  ! Derivative in z
  real(PP) function zderiv(array,i,j,k)
    real(PP), dimension(:,:,:) :: array
    integer :: i,j,k
    real(PP) :: idz
    idz=1.0d0/dz
    if ((k /= 1) .and. (k /= dimz)) then
      zderiv=0.5d0*idz*(array(i,j,k+1)-array(i,j,k-1))
    else if (k .eq. 1) then
      zderiv=(-3.0d0*array(i,j,1)+4.0d0*array(i,j,2)-array(i,j,3))/(2.0d0*dz)
    else if (k .eq. dimz) then
      zderiv=(3.0d0*array(i,j,dimz)-4.0d0*array(i,j,dimz-1)+array(i,j,dimz-2))/(2.0d0*dz)
    end if
  end function zderiv

subroutine ff_check(new_th) 
  real(PP),intent(out) :: new_th
  real(PP) :: jb,jmag,bmag,sin_theta,ratio ! Add 26 Nov 2014
  real(PP) :: sum_w,sum_jmag ! The variables used for storing the sumation are double precesion
  real(PP),dimension(3) :: Jc,Bc
  integer :: i,j,k
    
  ! Initalise sums
  sum_jmag = 0.0d0
  sum_w = 0.0d0
  
  ! Current-weighted average angle. This loop computes sum(|J|*sin(theta)),
  ! and sum(|J|) over the volume. The summation is performed in parallel 
  ! using the OpenMP reduction derective.
  !$omp parallel do private(i,j,k,jb,bmag,jmag,Jc,Bc,sin_theta,ratio) default(none) shared(dimx,dimy,dimz,bxyz) reduction(+:sum_w,sum_jmag)
  do i=1,dimx
  do j=1,dimy
  do k=1,dimz
    ! Compute Jc = curl(B) using centred differenes
    Jc(1) = yderiv(bxyz(1:dimx,1:dimy,1:dimz,3),i,j,k)-zderiv(bxyz(1:dimx,1:dimy,1:dimz,2),i,j,k)
    Jc(2) = zderiv(bxyz(1:dimx,1:dimy,1:dimz,1),i,j,k)-xderiv(bxyz(1:dimx,1:dimy,1:dimz,3),i,j,k)
    Jc(3) = xderiv(bxyz(1:dimx,1:dimy,1:dimz,2),i,j,k)-yderiv(bxyz(1:dimx,1:dimy,1:dimz,1),i,j,k)
    ! Extract B vector 
    Bc = bxyz(i,j,k,:)
    ! Compute dot product squared (J.dot.B)^2
    jb = dot_product(Bc,Jc)**2
    ! Compute magnitudes of |J| and |B|
    bmag = sqrt(dot_product(Bc,Bc))
    jmag = sqrt(dot_product(Jc,Jc))
    ! Compute angle between J and B. If (|J|*|B|)^2 is zero then
    ! sin(theta) is set to zero
    if((bmag*jmag)**2 .ne. 0) then
      ! Compute ratio (J.B/|J||B|)^2
      ratio = jb/(bmag*jmag)**2
      ! Check that ratio is not greater than 1. This can 
      ! can happend due to rounding error, and produces a NAN.
      ! Added 26 Nov 2014
      if(ratio .gt. 1) then
        sin_theta = 0.0d0
      else
        sin_theta = sqrt(1.0d0-ratio) 
      endif
    else
      sin_theta = 0.0d0
    endif
    ! Compute summation on this thread. The results for each
    ! thread are combined at the end of the loop because the 
    ! reduction command has been specified in the $omp directive.
    sum_jmag = sum_jmag + jmag
    sum_w = sum_w + sin_theta*jmag
  enddo
  enddo
  enddo  
  !$omp end parallel do

  ! Return with value. If sum(|Jc|) is zero (i.e. B is potential),
  ! then new_th=0 is returned. 
  if(sum_jmag .gt. 0) then
    new_th = (180.0d0/PI)*asin(sum_w/sum_jmag)
    return
  else
    new_th = 0.0d0
    return
  endif
end subroutine ff_check

  subroutine update_b  
    integer :: M
    real(PP) :: aab_chg
    real(PP), dimension(:,:,:,:), allocatable :: btmp
   
    allocate(btmp(dimx,dimy,dimz,3))
    btmp=bxyzj
    ! Update with under-relaxation
    bxyz(1:dimx,1:dimy,1:dimz,:)=(1.0-factor)*bxyz(1:dimx,1:dimy,1:dimz,:)+factor*(bxyz0+bxyzj)
    ! Calculate standard deviation of change
    M=dimx*dimy*dimz
    call fill_3d_vector(bxyz)
    aab_chg=sum(abs(bxyz(1:dimx,1:dimy,1:dimz,:)-btmp))/M
    write(*,'(A,g10.4)')'|                        Average absolute change in field: ',aab_chg
    deallocate(btmp)
  end subroutine update_b

  subroutine reduce_spike
    integer :: nspikes
    real(PP) :: std_cbcx0,std_cbcy0
    real(PP), dimension(:,:), allocatable :: bx0,by0
    real(PP), dimension(:,:), allocatable :: bcx0,bcy0,cbcx0,cbcy0
    ! Apply de-PIking, if needed
    ! Allocate memory
    allocate(bcx0(dimx,dimy),bcy0(dimx,dimy),cbcx0(dimx,dimy),cbcy0(dimx,dimy))
    allocate(bx0(dimx,dimy),by0(dimx,dimy))
    ! Define arrays
    bcx0=bxyzj(:,:,1,1)
    bcy0=bxyzj(:,:,1,2)
    bx0=bxyz(1:dimx,1:dimy,1,1)
    by0=bxyz(1:dimx,1:dimy,1,2)
    ! Calculate curvature
    call curvature(bcx0,cbcx0)
    call curvature(bcy0,cbcy0)
    ! Determine standard deviations of curvature over non-zero
    ! elements 
    call std(cbcx0,std_cbcx0)
    call std(cbcy0,std_cbcy0)      
    if (maxval(abs(cbcx0)) .gt. SPK*std_cbcx0 .or. &
        maxval(abs(cbcy0)) .gt. SPK*std_cbcy0) then 
      nspikes=count(abs(cbcx0) .gt. SPK*std_cbcx0)
      nspikes=nspikes+count(abs(cbcy0) .gt. SPK*std_cbcy0)
      write(*,'(A)')'De-PIking is required...' 
      write(*,'(A,g10.4)')'Number of points with PIkes = ',nspikes 
    endif
    ! Remove PIkes in bcx0, bcy0
    where (abs(cbcx0) .gt. SPK*std_cbcx0)
      ! bx0=0 ! Needed if under-relaxation is used
      bcx0=0
    end where
    where (abs(cbcy0) .gt. SPK*std_cbcy0)
      ! Needed if under-relaxation is used
      bcy0=0 
    end where
    ! Replace de-PIked values in bc
    bxyzj(:,:,1,1)=bcx0
    bxyzj(:,:,1,2)=bcy0
    bxyz(1:dimx,1:dimy,1,1)=bx0
    bxyz(1:dimx,1:dimy,1,2)=by0
    ! De-allocate memory
    deallocate(bcx0,bcy0,cbcx0,cbcy0)
  end subroutine reduce_spike

end module mod_operator
