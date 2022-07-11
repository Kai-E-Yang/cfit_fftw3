module mod_io
use mod_param
use mod_operator
implicit none
contains
  subroutine read_parameter()
    implicit none
    ! set default value    
    Bz0Name      = 'bz0.dat'
    AlphaName    = 'alpha0.dat'
    AlphaErrName = 'alpha0err.dat'
    OutFileName  = './result/'
    indataformat = 'binary'
    method       = 'fft'
    RestartName  = 'restartB.dat'


    NumThreads = 1
    ncycle     = 1
    dimx       = 128
    dimy       = 128
    dimz       = 128
    delta_s    = 0.5
    min_step   = 0.05
    check      = .false.
    factor     = 0.5
    restart    = .false.
    savNum     = 0
    top_closed = .true.
    SPK_flag   = .false.
    SPK        = 1e-5
    Periodic   = .false.
    stepmethod = 'rk2'
    Aout       = .true.
    Polarity   = 1
    nloop      = 1
    startcyc   = 1
    startloop  = 1
    alpha_error=.false.
    np_lenght_weight=.false.
    precision_flag = 'double'

    ! read parameters   
    write(*,'(A)')'| = = = = = = = = = = = = = = = = = = = = = = = = = = = = ='
    write(*,'(A)')'| Reading parameters from file: '//trim(par)
    open(Unit=1,file=par,STATUS='OLD',ACTION='READ')
    read(1,nml=filename_par,end=3)
3   rewind(1)
    read(1,nml=cal_par,end=5)
5   close(1)

    open(Unit=1,file=trim(OutFileName)//'input.par',status='replace',action='write')
    write(Unit=1, nml=filename_par)
    write(Unit=1, nml=cal_par)
    close(1)

    if(savNum.eq.0) savNum=ncycle
    logopened=.false.
  end subroutine read_parameter

  subroutine grid    
    integer :: i,j,k
    dx=real(2*PI,PP)/real(dimx,PP)
    dy = dx
    dz = dx 
  
    do i = 1,dimx
      xarr(i) = dx*(i-real(1,PP))
    enddo
    do j = 1,dimy
      yarr(j) = dy*(j-real(1,PP))
    enddo
    do k = 1,dimz
      zarr(k) = dz*(k-real(1,PP))
    enddo

  end subroutine grid

  subroutine allocate_var()

    allocate(axyz(dimx,dimy,dimz,3))
    allocate(bxyz(0:dimx+1,0:dimy+1,0:dimz+1,3))
    allocate(bxyz0(dimx,dimy,dimz,3))
    allocate(bxyzj(dimx,dimy,dimz,3))
    allocate(jxyz(dimx,dimy,dimz,3))

    allocate(alpha_p(dimx,dimy,dimz))
    allocate(alpha_n(dimx,dimy,dimz))

    allocate(alpha(dimx,dimy,dimz))    
    allocate(left_box(dimx,dimy,dimz))

    allocate(bz0(dimx,dimy))
    allocate(bzL(dimx,dimy))
    allocate(bz0n(dimx,dimy))
    allocate(bz0t(dimx,dimy))
    allocate(xarr(dimx))
    allocate(yarr(dimy))
    allocate(zarr(dimz))

    allocate(alpha0(0:dimx+1,0:dimy+1))
    allocate(sig0(dimx,dimy))
    allocate(sig00(dimx,dimy))
    allocate(sig0p(dimx,dimy))
    allocate(sig0n(dimx,dimy))

    allocate(fjc(dimx/2+1,dimy,dimz,3))
    allocate(fbc(dimx/2+1,dimy,dimz,3))

    if(np_lenght_weight)then
      allocate(weight_p(dimx,dimy,dimz))
      allocate(weight_n(dimx,dimy,dimz))
    endif

    if(Aout) then
      allocate(fap(dimx/2+1,dimy,dimz,3))
    endif

    iterMax = 5
    BoundaryEps = real(1,PP)+min_step
    call grid    
    scale_factor=real(1,PP)/real(dimx,PP)/real(dimy,PP)

  end subroutine allocate_var

  subroutine deallocate_var()
    deallocate(axyz)
    deallocate(bxyz)
    deallocate(bxyz0)
    deallocate(bxyzj)
    deallocate(jxyz)

    deallocate(alpha_p)
    deallocate(alpha_n)
    deallocate(alpha)    
    deallocate(left_box)

    deallocate(bz0)
    deallocate(bzL)
    deallocate(bz0n)
    deallocate(bz0t)
    deallocate(xarr)
    deallocate(yarr)
    deallocate(zarr)

    deallocate(alpha0)
    deallocate(sig0)
    deallocate(sig00)
    deallocate(sig0p)
    deallocate(sig0n)

    deallocate(fjc)
    deallocate(fbc)
    if(Aout) then
      deallocate(fap)
    endif
  end subroutine deallocate_var

  subroutine show_start()
    write(*,'(A)')'| = = = = = = = = = = = = = = = = = = = = = = = = = = = = ='
    write(*,'(A)')'|               _____   _____   ___   _____               |'
    write(*,'(A)')'|              |  ___| |  ___| |   | |_   _|              |'
    write(*,'(A)')'|              | |     | |__    | |    | |                |'
    write(*,'(A)')'|              | |     |  __|   | |    | |                |'
    write(*,'(A)')'|              | |___  | |      | |    | |                |'
    write(*,'(A)')'|              |_____| |_|     |___|   |_|                |'
    write(*,'(A)')'|                                                         |'

  end subroutine show_start

  subroutine show_information()
    write(*,'(A)')'|--------------------------------------------------------------------'
    write(*,'(A)')'| Mission: Calculate Non-linear force-free field by Grad-Rubin method'
    write(*,'(A)')'| The method based on Wheatland 2006, 2007'
    if (restart)then
    write(*,'(A)')'| 3D B0 field data comes from: '
    write(*,'(A)')'|                    '//trim(RestartName)
    endif
    write(*,'(A)')'| 2D Bn data comes from: '
    write(*,'(A)')'|                    '//trim(Bz0Name)
    write(*,'(A)')'| 2D Alpha data comes from: '
    write(*,'(A)')'|                    '//trim(AlphaName)
    write(*,'(A)')'| The result is saved in file: '
    write(*,'(A)')'|                    '//trim(OutFileName)
    write(*,'(A)')'| '
    write(*,'(A)')'|--------------------------------------------------------'
    write(*,'(A)')'| The NDIM of the Magnetic data is: '
    write(*,'(A,I3)')'|                  x dimension: ',dimx
    write(*,'(A,I3)')'|                  y dimension: ',dimy
    write(*,'(A,I3)')'|                  z dimension: ',dimz
    write(*,'(A)')'| '
    write(*,'(A,I4,A)')'| This mission use ',NumThreads,' threads'
    write(*,'(A)')'|--------------------------------------------------------'
  end subroutine show_information

  subroutine read_b0
    logical          :: alive
    character(len=144):: readname
    real(SP),dimension(:,:,:,:),allocatable :: bxyz0_tmp
    write(readname,'(A,i4.4,A,i4.4,''.dat'')') trim(OutFileName)//"B_",0,'_',0
    inquire(file=readname,exist=alive)
    if (alive) then
      if (trim(precision_flag).eq.'float')then
        allocate(bxyz0_tmp(1:dimx,1:dimy,1:dimz,1:3))
        write(*,'(A)')'| Loading B0 ...'
        open(4,File=readname,Access="stream",Form = "unformatted" )
        write(4) bxyz0_tmp(1:dimx,1:dimy,1:dimz,1:3)
        close(4)
        write(*,'(A)')'| Data sucessfully loaded !'
        bxyz0(1:dimx,1:dimy,1:dimz,1:3)=real(bxyz0_tmp(1:dimx,1:dimy,1:dimz,1:3),PP)
        deallocate(bxyz0_tmp)
      else   
        write(*,'(A)')'| Loading B0 ...'
        open(4,File=readname,Access="stream",Form = "unformatted" )
        write(4) bxyz0(1:dimx,1:dimy,1:dimz,1:3)
        close(4)
        write(*,'(A)')'| Data sucessfully loaded !'
      endif
    else
      write(*,'(A)')'| File '//trim(readname)//' does not exist'
      stop
    endif
  end subroutine read_b0

  subroutine read_b_restart
    logical          :: alive
    real(SP),dimension(:,:,:,:),allocatable :: bxyz0_tmp
    if(restart) then
      inquire(file=RestartName,exist=alive)
      if(alive) then 
        if (trim(precision_flag).eq.'float')then
          allocate(bxyz0_tmp(1:dimx,1:dimy,1:dimz,1:3))
          write(*,'(A)')'| Loading restart B0...'//trim(RestartName)
          open(3,File=trim(RestartName),Access="stream",Form = "unformatted" )
          read(3) bxyz0_tmp(1:dimx,1:dimy,1:dimz,1:3)
          close(3)
          write(*,'(A)')'| Data sucessfully loaded !'
          bxyz0(1:dimx,1:dimy,1:dimz,1:3)=real(bxyz0_tmp(1:dimx,1:dimy,1:dimz,1:3),PP)
          deallocate(bxyz0_tmp)
        else
          write(*,'(A)')'| Loading restart B0...'//trim(RestartName)
          open(3,File=trim(RestartName),Access="stream",Form = "unformatted" )
          read(3) bxyz0(1:dimx,1:dimy,1:dimz,1:3)
          close(3)
          write(*,'(A)')'| Data sucessfully loaded !'
        endif
      else
        write(*,'(A)')'| File '//trim(RestartName)//' does not exist'
        stop
      end if
    end if
  end subroutine read_b_restart

  subroutine read_bcs()
    logical          :: alive
    real(SP),dimension(:,:),allocatable :: alpha0_tmp,bz0_tmp,sig0_tmp

    inquire(file=AlphaName,exist=alive)
    if(alive) then
      alpha0=alpha0*real(0,PP)
      if (trim(precision_flag).eq.'float')then
        allocate(alpha0_tmp(1:dimx,1:dimy))
        write(*,'(A)')'| Loading Alpha field data...  '//trim(AlphaName)
        open(3,File=AlphaName,Access="stream",Form = "unformatted" )
        read(3) alpha0_tmp(1:dimx,1:dimy)
        close(3)
        write(*,'(A)')'| Data sucessfully loaded !'
        alpha0(1:dimx,1:dimy)=real(alpha0_tmp(1:dimx,1:dimy),PP)
        deallocate(alpha0_tmp)
      else
        write(*,'(A)')'| Loading Alpha field data...  '//trim(AlphaName)
        open(3,File=AlphaName,Access="stream",Form = "unformatted" )
        read(3) alpha0(1:dimx,1:dimy)
        close(3)
        write(*,'(A)')'| Data sucessfully loaded !'
      endif
    else
      write(*,'(A)')'| File '//AlphaName//' does not exist'
      stop
    end if

    inquire(file=Bz0Name,exist=alive)
    if(alive) then
      if (trim(precision_flag).eq.'float')then
        allocate(bz0_tmp(dimx,dimy))
        write(*,'(A)')'| Loading Bz bottom field data...  '//trim(Bz0Name)
        open(3,File=Bz0Name,Access="stream",Form = "unformatted" )
        read(3) bz0_tmp
        close(3)
        write(*,'(A)')'| Data sucessfully loaded !'
        bz0=real(bz0_tmp,PP)
        deallocate(bz0_tmp)
      else
        write(*,'(A)')'| Loading Bz bottom field data...  '//trim(Bz0Name)
        open(3,File=Bz0Name,Access="stream",Form = "unformatted" )
        read(3) bz0
        close(3)
        write(*,'(A)')'| Data sucessfully loaded !'
      endif
    else
      write(*,'(A)')'| File '//Bz0Name//' does not exist'
      stop
    end if
    call fill_edge_vertex_scalar(alpha0)

    if(alpha_error) then
      inquire(file=AlphaErrName,exist=alive)
      if(alive) then
        sig0=sig0*real(0,PP)
        if (trim(precision_flag).eq.'float')then
          allocate(sig0_tmp(1:dimx,1:dimy))
          write(*,'(A)')'| Loading Alpha field data...  '//trim(AlphaErrName)
          open(3,File=AlphaErrName,Access="stream",Form = "unformatted" )
          read(3) sig0(1:dimx,1:dimy)
          close(3)
          write(*,'(A)')'| Data sucessfully loaded !'
          sig0(1:dimx,1:dimy)=real(sig0_tmp(1:dimx,1:dimy),PP)
          deallocate(sig0_tmp)
        else
          write(*,'(A)')'| Loading Alpha field data...  '//trim(AlphaErrName)
          open(3,File=AlphaErrName,Access="stream",Form = "unformatted" )
          read(3) sig0(1:dimx,1:dimy)
          close(3)
          write(*,'(A)')'| Data sucessfully loaded !'
        endif
        max_sig0=maxval(sig00)

        if (alpha_error) then
          where(bz0 .gt. 0)
            sig0p=sig0
          elsewhere
            sig0p=max_sig0
          endwhere
          where(bz0 .lt. 0)
            sig0n=sig0
          elsewhere
            sig0n=max_sig0
          endwhere
        endif

      else
        write(*,'(A)')'| File '//AlphaErrName//' does not exist'
        stop
      end if
    endif

  end subroutine read_bcs

  ! subroutine write_data(icyc)
  !   integer :: icyc
  !   character(len=144):: savename

  !   write(savename,'(A,i4.4,''.dat'')') trim(OutFileName)//"_",icyc
  !   write(*,'(A)')'| Writing the result ...'
  !   open(4,File=savename,Access="stream",STATUS="REPLACE",&
  !   &Form = "unformatted" )
  !   write(4) bxyz(1:dimx,1:dimy,1:dimz,3)
  !   close(4)
  ! end subroutine write_data

  subroutine write_a(icyc,iloop)
    integer :: icyc,iloop
    character(len=144):: savename

    write(savename,'(A,i4.4,A,i4.4,''.dat'')') trim(OutFileName)//"A_",icyc,'_',iloop
    write(*,'(A)')'| Writing the result ...'
    open(4,File=savename,Access="stream",STATUS="REPLACE",&
    &Form = "unformatted" )
    if (trim(precision_flag).eq.'float')then
      write(4) real(axyz(1:dimx,1:dimy,1:dimz,:),SP)
    else
      write(4) axyz(1:dimx,1:dimy,1:dimz,:)
    endif
    close(4)
  end subroutine write_a

  subroutine write_b(icyc,iloop)
    integer :: icyc,iloop
    character(len=144):: savename

    write(savename,'(A,i4.4,A,i4.4,''.dat'')') trim(OutFileName)//"B_",icyc,'_',iloop
    write(*,'(A)')'| saving the B ...'
    open(4,File=savename,Access="stream",STATUS="REPLACE",&
    &Form = "unformatted" )
    if (trim(precision_flag).eq.'float')then
      write(4) real(bxyz(1:dimx,1:dimy,1:dimz,1:3),SP)
    else
      write(4) bxyz(1:dimx,1:dimy,1:dimz,1:3)
    endif    
    close(4)
  end subroutine write_b

  subroutine write_bxyzj(icyc,iloop)
    integer :: icyc,iloop
    character(len=144) :: savename

    write(savename,'(A,i4.4,A,i4.4,''.dat'')') trim(OutFileName)//"Bc_",icyc,'_',iloop
    write(*,'(A)')'| saving the Bc ...'
    open(4,File=savename,Access="stream",STATUS="REPLACE",&
    &Form = "unformatted" )
    if (trim(precision_flag).eq.'float')then
      write(4) real(bxyzj,SP)
    else
      write(4) bxyzj
    endif
    close(4)
  end subroutine write_bxyzj

  subroutine write_alpha(icyc,iloop)
    integer :: icyc,iloop
    character(len=144) :: savename

    write(savename,'(A,i4.4,A,i4.4,''.dat'')') trim(OutFileName)//"Alpha3d_",icyc,'_',iloop
    write(*,'(A)')'| saving the alpha ...'
    open(4,File=savename,Access="stream",STATUS="REPLACE",&
    &Form = "unformatted" )
    if (trim(precision_flag).eq.'float')then
      write(4) real(alpha,SP)
    else
      write(4) alpha
    endif
    close(4)
  end subroutine write_alpha

  subroutine write_alpha0(icyc,iloop)
    integer :: icyc,iloop
    character(len=144) :: savename

    write(savename,'(A,i4.4,A,i4.4,''.dat'')') trim(OutFileName)//"Alpha2d_",icyc,'_',iloop
    write(*,'(A)')'| saving the bottom alpha ...'
    open(4,File=savename,Access="stream",STATUS="REPLACE",&
    &Form = "unformatted" )
    if (trim(precision_flag).eq.'float')then
      write(4) real(alpha0,SP)
    else
      write(4) alpha0
    endif
    close(4)
  end subroutine write_alpha0

  subroutine write_leftbox(icyc,iloop)
    integer :: icyc,iloop
    character(len=144) :: savename

    write(savename,'(A,i4.4,A,i4.4,''.dat'')') trim(OutFileName)//"LeftBox_",icyc,'_',iloop
    write(*,'(A)')'| saving the left box ...'
    open(4,File=savename,Access="stream",STATUS="REPLACE",&
    &Form = "unformatted" )
    if (trim(precision_flag).eq.'float')then
      write(4) real(left_box,SP)
    else
      write(4) left_box
    endif
    close(4)
  end subroutine write_leftbox

  subroutine write_jxyz(icyc,iloop)
    integer :: icyc,iloop
    character(len=144) :: savename

    write(savename,'(A,i4.4,A,i4.4,''.dat'')') trim(OutFileName)//"Jc_",icyc,'_',iloop
    write(*,'(A)')'| saving the Jxyz ...'
    open(4,File=savename,Access="stream",STATUS="REPLACE",&
    &Form = "unformatted" )
    if (trim(precision_flag).eq.'float')then
      write(4) real(jxyz,SP)
    else
      write(4) jxyz
    endif
    close(4)
  end subroutine write_jxyz

  subroutine printlog(icyc,new_th,new_en)
    integer :: icyc
    real(PP) :: new_th,new_en
    character(len=800) :: filename,filehead
    character(len=2048) :: line,datastr

    write(filename,'(a)') trim(OutFileName)//"/cfit_log.csv"
    if(.not.logopened) then
      ! generate filename
      logopened=.true.
      filehead="iter,theta,energy"
      line=''
      write(datastr,'(i8,a)') icyc,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6,a)') new_th,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6)') new_en
      line=trim(line)//trim(datastr)

      open(7,file=trim(filename),access='sequential',position='append',status="replace")
      write(7,'(A17)') trim(filehead)
      write(7,'(A37)') line
      close(7)
    else
      line=''
      write(datastr,'(i8,a)') icyc,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6,a)') new_th,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6)') new_en
      line=trim(line)//trim(datastr)
      open(7,file=trim(filename),access='sequential',position='append',status="old")
      write(7,'(A37)') line
      close(7)
    endif
  end subroutine printlog

end module mod_io
