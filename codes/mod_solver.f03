module mod_solver
use mod_param
use mod_operator
use mod_io
use omp_lib
use mod_fftw
implicit none
contains
! ------------- master -----------
  subroutine master()
    integer :: ClockStart,ClockEnd,LoopClockStart,LoopClockEnd,ClockRate,ClockMax,Ticks
    integer :: IterCYC
    real(PP):: Seconds,NewTheta,NewEnergy
    real(PP):: NewDiv

    call system_clock(count_max=ClockMax,count_rate=ClockRate)
    call system_clock(ClockStart)
    call allocate_var()
    call read_bcs()
    call show_information()

    write(*,'(A)')'| Initialization ...'
    if(restart) then
      write(*,'(A)')'| Restart from file: '//trim(RestartName)
      write(*,'(A,i4.4,A,i4.4)')'|         restart cycle: ',startcyc,', loop',startloop
      call read_b_restart
      bxyz(1:dimx,1:dimy,1:dimz,:)=bxyz0
      call fftw3_laplace
      call fill_3d_vector(bxyz)
      call calc_energy(NewEnergy)
      call ff_check(NewTheta)
      call div_check(NewDiv)

    else
      write(*,'(A)')'| Starting a new calculation'
      call fftw3_laplace
      call write_b(0,0)
      call calc_energy(NewEnergy)
      call ff_check(NewTheta)
      call div_check(NewDiv)

      if(check) then
        call write_bxyzj(0,0)
        call write_leftbox(0,0)
        call write_alpha(0,0)
        call write_jxyz(0,0)
      endif
      if(Aout) call write_a(0,0)
    endif

    E0=NewEnergy

    if(np_lenght_weight) then
      if(restart) then
        call read_b_restart
        bxyz(1:dimx,1:dimy,1:dimz,:)=bxyz0
      else
        call read_b0
        bxyz(1:dimx,1:dimy,1:dimz,:)=bxyz0
      endif
      call do_weight_np_loop(0)
      call system_clock(LoopClockStart)
    else
      do IterCYC=startcyc,ncycle
        if(IterCYC .gt. 1)then
          if(restart) then
            call read_b_restart
            bxyz(1:dimx,1:dimy,1:dimz,:)=bxyz0
          else
            call read_b0
            bxyz(1:dimx,1:dimy,1:dimz,:)=bxyz0
          endif
        endif
        call system_clock(LoopClockStart)
        if (ncycle.eq.1) then
          write(*,'(A,i5.5,g10.4,A)')'| Cycle: ',IterCYC,100.0*IterCYC/ncycle,'%'
        endif
        ! if (ncycle.gt.1) Polarity=Polarity*(-1)**(IterCYC+1)
        if (IterCYC.gt.1) Polarity=Polarity*(-1)

        if(IterCYC.eq.1 .or. restart)then
          call read_bcs
        endif
        if (Polarity.gt.0)then
          call do_positive_loop(IterCYC)
        else
          call do_negative_loop(IterCYC)
        endif
        if(IterCYC.ne.ncycle)then
          call system_clock(LoopClockEnd)
          Ticks=LoopClockEnd-LoopClockStart
          Seconds=float(Ticks)/float(ClockRate)
          write(*,'(A,g10.4,A)')'|  This polarity cycle tooks',Seconds,' sec'
        endif
        if( mod(IterCYC,2).eq.0)then
          call combine_alpha
          call write_alpha0(IterCYC)
        endif
      enddo
    endif
    call system_clock(ClockEnd)
    Ticks=ClockEnd-ClockStart
    Seconds=float(Ticks)/float(ClockRate)
    write(*,'(A,g10.4,A)')'| Total calculation tooks',Seconds,' sec'
  end subroutine master

  subroutine do_positive_loop(IterCYC)
    integer,intent(in)::IterCYC
    integer :: IterLoop
    real(PP):: NewTheta,NewEnergy
    real(PP):: NewDiv
    integer :: IterNow
    real(PP):: OldTheta,OldEnergy,DeltaTheta,DeltaEnergy

    OldTheta=0.0d0
    OldEnergy=0.0d0
    IterNow=(IterCYC-1)*nloop
    do IterLoop=1,nloop
      write(*,'(A)')'.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .' 
      write(*,'(A)')'|  ' 
      write(*,'(A)')'| Positive solution ... cycle and loop'
      write(*,'(A,i4.4,A,i4.4)')'|                        ',IterCYC,', ',IterLoop
      call calc_alpha()
      alpha=alpha_p

      call construct_jc()
      if(top_closed)then
        call fftw3_poisson_top_closed()
      else
        call fftw3_poisson_top_open()
      endif
      if(SPK_flag) call reduce_spike
      call update_b

      if(mod(IterLoop,savNum).eq.0 .or. IterLoop .eq. nloop)then
        call write_b(IterCYC,IterLoop)
        if(check) then
          call write_bxyzj(IterCYC,IterLoop)
          call write_leftbox(IterCYC,IterLoop)
          call write_alpha(IterCYC,IterLoop)
          call write_jxyz(IterCYC,IterLoop)
        endif
        if(Aout) call write_a(IterCYC,IterLoop)
      endif
      IterNow=IterNow+1
      call ff_check(NewTheta)
      call div_check(NewDiv)
      call calc_energy(NewEnergy)
      call printlog(IterNow,NewTheta,NewEnergy/E0,NewDiv)
      if(IterLoop.eq.1)then
        DeltaTheta=0.0d0
        DeltaEnergy=0.0d0
      else
        DeltaTheta=NewTheta-OldTheta
        DeltaEnergy=NewEnergy-OldEnergy
      endif
      OldTheta=NewTheta
      OldEnergy=NewEnergy
      write(*,'(A,g10.4)')'|                        Normalized Energy = ',NewEnergy/E0
      write(*,'(A,g10.4)')'|                        Average angle = ',NewTheta
      write(*,'(A,g10.4)')'|                        Integral div = ',NewDiv
      write(*,'(A,g10.4)')'|                        delta theta = ',DeltaTheta
      write(*,'(A,g10.4)')'|                        delta E = ',DeltaEnergy
      write(*,'(A)')' '
    enddo
  end subroutine do_positive_loop

  subroutine do_negative_loop(IterCYC)
    integer,intent(in)::IterCYC
    integer :: IterLoop
    real(PP):: NewTheta,NewEnergy,NewDiv
    integer :: IterNow
    real(PP):: OldTheta,OldEnergy,DeltaTheta,DeltaEnergy

    OldTheta=0.0d0
    OldEnergy=0.0d0
    IterNow=(IterCYC-1)*nloop
    do IterLoop=1,nloop
      write(*,'(A)')'.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .'
      write(*,'(A)')'|  '  
      write(*,'(A)')'| Negative solution ... cycle and loop'
      write(*,'(A,i4.4,A,i4.4)')'|                        ',IterCYC,', ',IterLoop
      call calc_alpha()
      alpha=alpha_n

      call construct_jc()
      if(top_closed)then
        call fftw3_poisson_top_closed()
      else
        call fftw3_poisson_top_open()
      endif
      if(SPK_flag) call reduce_spike
      call update_b

      if(mod(IterLoop,savNum).eq.0 .or. IterLoop .eq. nloop)then
        call write_b(IterCYC,IterLoop)
        if(check) then
          call write_bxyzj(IterCYC,IterLoop)
          call write_leftbox(IterCYC,IterLoop)
          call write_alpha(IterCYC,IterLoop)
          call write_jxyz(IterCYC,IterLoop)
        endif
        if(Aout) call write_a(IterCYC,IterLoop)
      endif
      IterNow=IterNow+1
      call ff_check(NewTheta)
      call div_check(NewDiv)
      call calc_energy(NewEnergy)
      call printlog(IterNow,NewTheta,NewEnergy/E0,NewDiv)
      if(IterLoop.eq.1)then
        DeltaTheta=0.0d0
        DeltaEnergy=0.0d0
      else
        DeltaTheta=NewTheta-OldTheta
        DeltaEnergy=NewEnergy-OldEnergy
      endif
      OldTheta=NewTheta
      OldEnergy=NewEnergy
      write(*,'(A,g10.4)')'|                        Normalized Energy = ',NewEnergy/E0
      write(*,'(A,g10.4)')'|                        Average angle = ',NewTheta
      write(*,'(A,g10.4)')'|                        Integral div = ',NewDiv
      write(*,'(A,g10.4)')'|                        delta theta = ',DeltaTheta
      write(*,'(A,g10.4)')'|                        delta E = ',DeltaEnergy
      write(*,'(A)')' '
    enddo
  end subroutine do_negative_loop

  subroutine do_weight_np_loop(IterCYC)
    integer,intent(in)::IterCYC
    integer :: IterLoop
    real(PP):: NewTheta,NewEnergy,NewDiv
    integer :: IterNow
    real(PP):: OldTheta,OldEnergy,DeltaTheta,DeltaEnergy

    OldTheta=0.0d0
    OldEnergy=0.0d0
    IterNow=(IterCYC-1)*nloop
    do IterLoop=1,nloop
      write(*,'(A)')'.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .' 
      write(*,'(A)')'|  ' 
      write(*,'(A)')'| Positive solution ... cycle and loop'
      write(*,'(A,i4.4,A,i4.4)')'|                        ',IterCYC,', ',IterLoop
      call calc_alpha()
      alpha=alpha_p*weight_p+alpha_n*weight_n

      call construct_jc()
      if(top_closed)then
        call fftw3_poisson_top_closed()
      else
        call fftw3_poisson_top_open()
      endif
      if(SPK_flag) call reduce_spike
      call update_b

      if(mod(IterLoop,savNum).eq.0 .or. IterLoop .eq. nloop)then
        call write_b(IterCYC,IterLoop)
        if(check) then
          call write_bxyzj(IterCYC,IterLoop)
          call write_leftbox(IterCYC,IterLoop)
          call write_alpha(IterCYC,IterLoop)
          call write_jxyz(IterCYC,IterLoop)
        endif
        if(Aout) call write_a(IterCYC,IterLoop)
      endif
      IterNow=IterNow+1
      call ff_check(NewTheta)
      call div_check(NewDiv)
      call calc_energy(NewEnergy)
      call printlog(IterNow,NewTheta,NewEnergy/E0,NewDiv)
      if(IterLoop.eq.1)then
        DeltaTheta=0.0d0
        DeltaEnergy=0.0d0
      else
        DeltaTheta=NewTheta-OldTheta
        DeltaEnergy=NewEnergy-OldEnergy
      endif
      OldTheta=NewTheta
      OldEnergy=NewEnergy
      write(*,'(A,g10.4)')'|                        Normalized Energy = ',NewEnergy/E0
      write(*,'(A,g10.4)')'|                        Average angle = ',NewTheta
      write(*,'(A,g10.4)')'|                        Integral div = ',NewDiv
      write(*,'(A,g10.4)')'|                        delta theta = ',DeltaTheta
      write(*,'(A,g10.4)')'|                        delta E = ',DeltaEnergy
      write(*,'(A)')' '
    enddo
  end subroutine do_weight_np_loop

!-------- the main loop for calculation --------
  subroutine do_loop_test()
    integer :: iterN
    real(PP) :: NewTheta,NewEnergy,NewDiv

    call fftw3_laplace
    call ff_check(NewTheta)
    call div_check(NewDiv)
    ! call calc_angle(NewTheta)
    call calc_energy(E0)
    NewEnergy=E0
    print *,'               IterN',0
    print *,'               Average angle (weighted) = ',NewTheta
    print *,'               Energy/E0 = ',NewEnergy/E0

    call printlog(0,NewTheta,1.0d0,NewDiv)

    do iterN=1,ncycle
      write(*,'(A)')'| '
      print *,'               IterN',IterN
      call calc_alpha()
      alpha=alpha_n
      call construct_jc()
      if(top_closed)then
        call fftw3_poisson_top_closed()
      else
        call fftw3_poisson_top_open()
      endif
      if(SPK_flag) call reduce_spike
      call update_b

      call ff_check(NewTheta)
      call div_check(NewDiv)
      call calc_energy(NewEnergy)
      call printlog(iterN,NewTheta,NewEnergy/E0,NewDiv)
      print *,'               Average angle (weighted) = ',NewTheta
      print *,'               Energy/E0 = ',NewEnergy/E0
      call write_b(1,iterN)
      call write_bxyzj(1,iterN)
    enddo
  end subroutine do_loop_test

  subroutine calc_alpha()
    integer :: k
    integer,dimension(8) :: time_begin,time_end
    real(PP) :: time_delta(8)

    alpha=alpha*real(0,PP)
    call date_and_time(VALUES=time_begin)
    call OMP_SET_NUM_THREADS(NumThreads)
    !$omp parallel do
    do k=1,dimz
      call calc_plane(k)
    end do
    !$omp end parallel do
    call date_and_time(VALUES=time_end)
    time_delta=time_end-time_begin
    write(*,'(A,g10.4)')'| field line tracing costs (min): '&
         & ,(time_delta(7)/3600.0 + time_delta(6)/60. &
         & + time_delta(5) + time_delta(3)*24.0)*60
  end subroutine calc_alpha
! ------- do parallel cal on each layer ------
  subroutine calc_plane(k)
    integer,intent(in)::k
    integer :: i,j
    real(PP) :: ap_slice(dimx,dimy),an_slice(dimx,dimy),m_slice(dimx,dimy)
    real(PP) :: wp_slice(dimx,dimy),wn_slice(dimx,dimy)
    real(PP) :: xindex,yindex,zindex

    zindex = real(k,PP)
    do j=1,dimy
    do i=1,dimx
      xindex = real(i,PP)
      yindex = real(j,PP)
      if (np_lenght_weight) then
        call calc_point_weight((/xindex,yindex,zindex/),ap_slice(i,j),an_slice(i,j),wp_slice(i,j),wn_slice(i,j),m_slice(i,j))
      else        
        call calc_point((/xindex,yindex,zindex/),ap_slice(i,j),an_slice(i,j),m_slice(i,j))
      endif
    end do
    end do

    left_box(1:dimx,1:dimy,k)=m_slice
    alpha_p(1:dimx,1:dimy,k)=ap_slice
    alpha_n(1:dimx,1:dimy,k)=an_slice

    if(np_lenght_weight)then
      weight_p(1:dimx,1:dimy,k)=wp_slice
      weight_n(1:dimx,1:dimy,k)=wn_slice
    endif
  end subroutine calc_plane

! ------ do calculation at each point ------
  subroutine calc_point(Posi,ap,an,mark)
    real(PP),intent(in) :: Posi(3)
    real(PP),intent(out) :: ap,an
    real(PP),dimension(neqn) :: LineF,LineB
    real(PP) :: lengthF,lengthB
    real(PP) :: markF,markB,mark

    call integralLine(Posi,LineF,real(1,PP),lengthF) 
    call integralLine(Posi,LineB,real(-1,PP),lengthB)
    ! ------------ get the mark flag for the two end ---------
    ! meaning of the value of end_mark
    ! 1: close field line
    ! 2: open field line, positive end roots in bottom
    ! 3: open field line, negative end roots in bottom
    ! 4: field line without end roots in the bottom boundary
    if(abs(LineF(3)) .lt. BoundaryEps .and. &
      abs(LineB(3)) .lt. BoundaryEps) then
      mark=real(1,PP)
      markF=1
      markB=1
    end if
    if(abs(LineB(3)) .lt. BoundaryEps .and. &
      abs(LineF(3)) .gt. BoundaryEps) then
      mark=real(2,PP)
      markF=0
      markB=1
    end if
    if(abs(LineF(3)) .lt. BoundaryEps .and. &
      abs(LineB(3)) .gt. BoundaryEps) then
      mark=real(3,PP)
      markF=1
      markB=0
    end if
    if(abs(LineF(3)) .gt. BoundaryEps .and. &
      abs(LineB(3)) .gt. BoundaryEps) then
      mark=real(4,PP)
      markF=0
      markB=0      
    end if

    call alpha_bottom(LineF,an)
    call alpha_bottom(LineB,ap)
    ! if(markF.lt.0.5) an=0.0d0
    ! if(markB.lt.0.5) ap=0.0d0

    if (.not.Periodic) then
      if (mark>1.5)then
        an=0.0d0
        ap=0.0d0
      endif
    endif

  end subroutine calc_point

  subroutine calc_point_weight(Posi,ap,an,wp,wn,mark)
    real(PP),intent(in) :: Posi(3)
    real(PP),intent(out) :: ap,an
    real(PP),dimension(neqn) :: LineF,LineB
    real(PP) :: lengthF,lengthB
    real(PP) :: markF,markB,mark
    real(PP) :: wp,wn

    call integralLine(Posi,LineF,real(1,PP),lengthF) 
    call integralLine(Posi,LineB,real(-1,PP),lengthB)
    ! ------------ get the mark flag for the two end ---------
    ! meaning of the value of end_mark
    ! 1: close field line
    ! 2: open field line, positive end roots in bottom
    ! 3: open field line, negative end roots in bottom
    ! 4: field line without end roots in the bottom boundary
    if(abs(LineF(3)) .lt. BoundaryEps .and. &
      abs(LineB(3)) .lt. BoundaryEps) then
      mark=real(1,PP)
      markF=1
      markB=1
      wp=lengthF/(lengthB+lengthF)
      wn=lengthB/(lengthB+lengthF)
    end if
    if(abs(LineB(3)) .lt. BoundaryEps .and. &
      abs(LineF(3)) .gt. BoundaryEps) then
      mark=real(2,PP)
      markF=0
      markB=1
      wp=real(1,PP)
      wn=real(0,PP)
    end if
    if(abs(LineF(3)) .lt. BoundaryEps .and. &
      abs(LineB(3)) .gt. BoundaryEps) then
      mark=real(3,PP)
      markF=1
      markB=0
      wp=real(0,PP)
      wn=real(1,PP)
    end if
    if(abs(LineF(3)) .gt. BoundaryEps .and. &
      abs(LineB(3)) .gt. BoundaryEps) then
      mark=real(4,PP)
      markF=0
      markB=0  
      wp=real(0,PP)
      wn=real(0,PP)    
    end if

    call alpha_bottom(LineF,an)
    call alpha_bottom(LineB,ap)
    if(markF.lt.0.5) an=0.0d0
    if(markB.lt.0.5) ap=0.0d0

  end subroutine calc_point_weight

  subroutine fieldline(LineP,LineF,LineB,lengthF,lengthB)
    real(PP)::LineP(neqn)
    real(PP)::LineF(neqn)
    real(PP)::LineB(neqn)
    real(PP)::sig,lengthF,lengthB
    sig = real(1,PP)
    call integralLine(LineP,LineF,sig,lengthF)    
    sig = real(-1,PP)
    call integralLine(LineP,LineB,sig,lengthB)
    lengthB=abs(lengthB)
    lengthF=abs(lengthF)
  end subroutine fieldline

  subroutine integralLine(LineP,LineI,sig,s_end)
    integer::n_step
    logical :: flag
    real(PP),intent(in)::sig
    real(PP)::LineP(neqn)
    real(PP)::LineI(neqn)
    real(PP)::tmp(neqn)
    real(PP)::s_start,s_end,ds
    
    s_start = real(0,PP)
    s_end = real(0,PP)
    
    tmp = LineP
    ds = real(delta_s,PP)*sig
    n_step = 0
    flag = .true.
    do while (n_step .le. isn .and. flag)
       s_end = s_start + ds
       select case(trim(stepmethod))
         case('rk2')
           call rk2(diffLine,neqn,s_start,s_end,tmp,flag)
         case('rk4')
           call rk4(diffLine,neqn,s_start,s_end,tmp,flag)
       end select
       if(Periodic) call periodic_posi(tmp,tmp)
       s_start = s_end
       n_step  = n_step+1
    end do
    LineI(1:3) = tmp
  end subroutine integralLine

  subroutine construct_jc()
    integer :: k
    !$omp parallel do
    do k=1,dimz
      jxyz(:,:,k,1)=alpha(:,:,k)*bxyz(1:dimx,1:dimy,k,1) 
      jxyz(:,:,k,2)=alpha(:,:,k)*bxyz(1:dimx,1:dimy,k,2) 
      jxyz(:,:,k,3)=alpha(:,:,k)*bxyz(1:dimx,1:dimy,k,3) 
    enddo
    !$omp end parallel do
  end subroutine construct_jc

end module mod_solver
