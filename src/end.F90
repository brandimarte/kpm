!  *******************************************************************  !
!  KPM Fortran Code 2014                                                !
!                                                                       !
!  By Pedro Brandimarte (brandimarte@gmail.com) and                     !
!     Eric de Castro e Andrade (eandrade@ift.unesp.br)                  !
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
!                              MODULE end                               !
!  *******************************************************************  !
!  Description: closes the program properly.                            !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !

MODULE end

!
! Modules
!
  use precision,       only: 
  use parallel,        only: 
  use init,            only: 
  use io,              only: 
  use hsparse,         only: 
  use moment,          only: 
  use kernel,          only: 
  use dct,             only: 

  implicit none

  PUBLIC  :: finalize
  PRIVATE ! default is private


CONTAINS


!  *******************************************************************  !
!                               finalize                                !
!  *******************************************************************  !
!  Description: closes the program properly.                            !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  real*8 time_begin           : Initial processor time in seconds      !
!  *******************************************************************  !
  subroutine finalize

!
! Modules
!
    use precision,       only: dp
    use parallel,        only: IOnode
    use init,            only: time_begin
    use io,              only: IOfree
    use hsparse,         only: Hfree
    use moment,          only: Mfree
    use kernel,          only: KERfree
    use dct,             only: DCTfree

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    real(dp) :: time_end
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines

    call MPI_Barrier (MPI_Comm_World, MPIerror)
#endif

    if (IOnode) write (6,'(/,34("*"),a,33("*"))') ' Ending KPM '

!   Free memory.
    if (IOnode) write (6,'(/,a)', ADVANCE='no')                         &
         'finalize: Freeing memory...'
    call IOfree
    call Hfree
    call Mfree
    call KERfree
    call DCTfree

    if (IOnode) write (6,'(a,/)') ' done!'
 
!   Finalizes MPI.
#ifdef MPI
    call MPI_Finalize (MPIerror)
#endif

    if (IOnode) then

!      Final time.
       call cpu_time (time_end)

       write (6,'(a,f12.4,a)') "Time of calculation was ",              &
            time_end - time_begin, " seconds"
       write (6,'(/,a,/)') "End of program KPM"

    endif


  end subroutine finalize


!  *******************************************************************  !


END MODULE end

