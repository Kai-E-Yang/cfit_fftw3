!=============================================================
! This program is written in Fortran
! @Date:   2019-01-11 13:06:36
! @Last Modified time: 2021-04-26 13:52:07
! @Version: V26
! @Licensing: Copy right from Kai E. Yang
! @Purpose: Modify the original CFIT code to the new fortran style
!=============================================================
program main
  use mod_param
  use mod_io
  use mod_solver
  use mod_operator
  implicit none

  call get_command_argument(1,par)
  call read_parameter()
  call show_start()
  call OMP_SET_NUM_THREADS(NumThreads)
  call master()
  write(*,'(A)')'| The calculation finished!'
  call deallocate_var()

end program main