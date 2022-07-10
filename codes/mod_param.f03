module mod_param
implicit none
  save
  integer, parameter :: i4= SELECTED_INT_KIND(R=9) 
  integer, parameter :: i8= SELECTED_INT_KIND(R=18) 
  integer, parameter :: r8= KIND(1.0)
  integer, parameter :: SP= KIND(1.0)
  integer, parameter :: DP= KIND(1.0D0)
  integer, parameter :: PP= DP
  integer, parameter :: neqn = 3


  character(len=144) :: par
  character(len=144) :: indataformat
  character(len=144) :: OutFileName, AlphaName
  character(len=144) :: AlphaErrName  
  character(len=144) :: Bz0Name,RestartName
  character(len=144) :: method,stepmethod
  character(len=144) :: precision_flag

  integer :: dimx,dimy,dimz,kmax
  integer :: NumThreads
  integer :: isn=10000
  integer :: iterMax
  integer :: ncycle,savNum
  integer :: Polarity,ERRCOUNT,nloop,startcyc,startloop

  ! here bxyz is the total magnetic filed, bxyz0 is the potential field
  ! bxyzj is the magnetic field only comes from the current
  ! axyz is the vector potential field for the potential field when calculating the potential field.
  ! axyz is the vector potential field for the current associated magnetic field when calculating this field.
  ! jxyz is the electric current induced by line tracing.
  real(PP),dimension(:,:,:,:),allocatable :: bxyz
  real(PP),dimension(:,:,:,:),allocatable :: bxyz0
  real(PP),dimension(:,:,:,:),allocatable :: bxyzj
  real(PP),dimension(:,:,:,:),allocatable :: axyz
  real(PP),dimension(:,:,:,:),allocatable :: jxyz
  real(PP),dimension(:,:,:),allocatable :: x,y,z
  real(PP),dimension(:,:,:),allocatable :: alpha,left_box
  real(PP),dimension(:,:,:),allocatable :: alpha_p,alpha_n
  real(PP),dimension(:,:,:),allocatable :: weight_p,weight_n


  real(PP),dimension(:,:),allocatable :: sig0,sig00,sig0p,sig0n
  real(PP),dimension(:,:),allocatable :: alpha0,bz0
  real(PP),dimension(:,:),allocatable :: bzL,bz0n,bz0t
  real(PP),dimension(:),allocatable :: xarr 
  real(PP),dimension(:),allocatable :: yarr
  real(PP),dimension(:),allocatable :: zarr

  complex(PP), dimension(:,:,:,:), allocatable :: fjc,fbc
  complex(PP), dimension(:,:,:,:), allocatable :: fap


  real(PP) :: PI=4*ATAN(real(1,PP))
  real(PP) :: TWOPI=8*ATAN(real(1,PP))
  real(PP) :: delta_s
  real(PP) :: E0,dx,dy,dz
  real(PP) :: BoundaryEps,factor
  real(PP) :: eps_B,min_step
  real(PP) :: SPK,max_sig0
  real(PP) :: scale_factor

  logical :: POT,LSQ,CHK,alpha_error,Aout,SPK_flag
  logical :: check,restart,top_closed,Periodic
  logical :: logopened,np_lenght_weight

  namelist /filename_par/ OutFileName, Bz0Name,&
   indataformat,AlphaName,AlphaErrName,RestartName
  namelist /cal_par/ NumThreads,dimx,dimy,dimz,delta_s,method,&
  stepmethod,ncycle,check,restart,savNum,factor,top_closed,&
  SPK_flag,SPK,Periodic,Aout,Polarity,nloop,startcyc,startloop,&
  alpha_error,min_step,np_lenght_weight,precision_flag

end module mod_param