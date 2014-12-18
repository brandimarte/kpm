!  *******************************************************************  !
!  KPM Fortran Code 2014                                                !
!                                                                       !
!  Written by Eric de Castro e Andrade (eandrade@ift.unesp.br) and      !
!             Pedro Brandimarte (brandimarte@gmail.com).                !
!                                                                       !
!  Copyright (c), All Rights Reserved                                   !
!                                                                       !
!  This program is free software. You can redistribute it and/or        !
!  modify it under the terms of the GNU General Public License          !
!  (version 3 or later) as published by the Free Software Foundation    !
!  <http://fsf.org/>.                                                   !
!                                                                       !
!  This program is distributed in the hope that it will be useful, but  !
!  WITHOUT ANY WARRANTY, without even the implied warranty of           !
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     !
!  General Public License for more details (file 'LICENSE_GPL'          !
!  distributed along with this program or at                            !
!  <http://www.gnu.org/licenses/gpl.html>).                             !
!  *******************************************************************  !
!                              MODULE init                              !
!  *******************************************************************  !
!  Description: intialize the program properly and read input options.  !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !

MODULE init

!
! Modules
!
  use precision,       only: dp
  use io,              only: 
  use options,         only: 
  use fdf

  implicit none

  PUBLIC  :: initialize, time_begin
  PRIVATE :: header, initread

! Initial processor time in seconds (processor-dependent approximation).
  real(dp) :: time_begin


CONTAINS


!  *******************************************************************  !
!                              initialize                               !
!  *******************************************************************  !
!  Description: initializes timer and the MPI execution environment.    !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ***************************** OUTPUT ******************************  !
!  real*8 time_begin           : Initial processor time in seconds      !
!  *******************************************************************  !
  subroutine initialize

!
! Modules
!
    use io,              only: IOinit
    use options,         only: OPTread, NumThreads

!   Print version information.
    call header

!   Initial time.
    call cpu_time (time_begin)

!   Init logical units control.
    call IOinit

!   Initialise read.
    call initread

!   Read simulation data.
    call OPTread

!   Set number of threads for parallel (multi-threaded) MKL subroutines.
    call mkl_set_num_threads (NumThreads)


  end subroutine initialize


!  *******************************************************************  !
!                               initread                                !
!  *******************************************************************  !
!  Description: initialize the reading of the data for KPM.             !
!                                                                       !
!  Uses FDF (Flexible Data Fromat) package of J.M.Soler and A.Garcia.   !
!                                                                       !
!  Taken from 'reinit' subroutine from SIESTA package.                  !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !
  subroutine initread

!
! Modules
!
    use io,              only: IOassign, IOclose
    use fdf

!   Local variables.
    character string*20
    character filein*20, fileout*20, line*150 
    integer :: count, length, lun, lun_tmp
    logical :: debug_input, file_exists

    write (6,'(/,27(1h*),a,27(1h*))') ' Dump of input data file '

!   Dump data file to output file and generate scratch file
!   for FDF to read from (except if INPUT_DEBUG exists).
    inquire (file='INPUT_DEBUG', exist=debug_input)
    if (debug_input) then
       write (6,'(a)') 'WARNING: ' //                                   &
            'KPM is reading its input from file INPUT_DEBUG'
           
       call IOassign(lun)
       filein = 'INPUT_DEBUG'
       open (lun, file='INPUT_DEBUG', form='formatted', status='old')
       rewind (lun)
    else
       write (6,'(a,/)') 'initread: Reading from standard input'
       lun = 5
       call IOassign (lun_tmp)
       do
          call system_clock (count)
          write (string,*) count
          filein = 'INPUT_TMP.' // adjustl(string)
          inquire (file=filein, exist=file_exists)
          if (.not. file_exists) exit
       end do
       open (lun_tmp, file=filein, form='formatted', status='replace')
       rewind (lun_tmp)
    endif

10  continue
    read (lun, err=20, end=20, fmt='(a)') line
    call chrlen (line, 0, length)
    if (length /= 0) then
       write(6,'(a,a)') 'initread: ', line(1:length)
       if (.not. debug_input) write (lun_tmp,'(a)') line(1:length)
    endif
    goto 10
20  continue

    write (6,'(2a)') 'initread: ', repeat('*', 69)

!   Choose proper file for fdf processing.
    if (debug_input) then
       call IOclose (lun)
    else
       call IOclose (lun_tmp)
    endif

!   Set up fdf.
    fileout = 'fdf.log'
    call fdf_init (filein, fileout)


  end subroutine initread


!  *******************************************************************  !
!                                header                                 !
!  *******************************************************************  !
!  Description: prints welcome message, date and copyright.             !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !
  subroutine header

!   Local variables.
    integer, dimension(8) :: values

    call date_and_time (VALUES=values)

    write (6,'(/,a,73(1h*),/)') '   '
    write (6,'(a,/)')                                                   &
         '                   *  WELCOME TO KPM CODE v2014.01  *'
    write (6,'(a,i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')               &
         '                            ',                                &
         values(1), '.', values(2), '.', values(3), ' , ',              &
         values(5), ':', values(6), ':', values(7)
    write (6,'(/,a)') '      Written by Eric de Castro e Andrade '      // &
         '(eandrade@ift.unesp.br) and'
    write (6,'(a)') '                 Pedro Brandimarte '            // &
         '(brandimarte@gmail.com).'
    write (6,'(/,a)') '      Copyright (c), All Rights Reserved'
    write (6,'(/,a)') '      This program is free software. '        // &
         'You can redistribute it and/or'
    write (6,'(a)') '      modify it under the terms of the '        // &
         'GNU General Public License'
    write (6,'(a)') '      (version 3 or later) as published '       // &
         'by the Free Software Foundation'
    write (6,'(a)') '      <http://fsf.org/>. See the GNU '          // &
         'General Public License for details.'
    write (6,'(/,a,73(1h*))') '   '


  end subroutine header


!  *******************************************************************  !


END MODULE init
