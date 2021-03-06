!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Green's functions.        !
!! Copyright (C) 2012                                                       !
!!                                                                          ! 
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Function calculation                             ! 
!!                                                                          !
!! Developers: Alessandro Pecchia, Gabriele Penazzi                         !
!! Former Conctributors: Luca Latessa, Aldo Di Carlo                        !
!!                                                                          !
!! libNEGF is free software: you can redistribute it and/or modify          !
!! it under the terms of the GNU Lesser General Public License as published !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !  
!!--------------------------------------------------------------------------!


                !$#define MEMLOG

module ln_allocation

  use ln_precision
  implicit none
  private

  integer, save :: iolog
  INTEGER(8), SAVE :: alloc_mem, peak_mem
  
  public :: log_allocate, log_deallocate
  public :: log_allocatep, log_deallocatep
  public :: writeMemInfo, writePeakInfo, getMem
  public :: resetMemLog, openMemLog, writeMemLog, closeMemLog
  public :: memstr

  interface log_allocatep
     module procedure allocate_pl, allocate_pd, allocate_pi, allocate_pz
     module procedure allocate_pd2, allocate_pi2
  end interface

  interface log_allocate
     module procedure allocate_l
     module procedure allocate_d, allocate_i, allocate_z
     module procedure allocate_d2, allocate_i2, allocate_z2
     module procedure allocate_d3, allocate_i3, allocate_z3
     module procedure allocate_d4
     module procedure allocate_c2
     module procedure allocate_c3
  end interface

  interface log_deallocatep
     module procedure deallocate_pl, deallocate_pd, deallocate_pi, deallocate_pz
     module procedure deallocate_pd2, deallocate_pi2
  end interface

  interface log_deallocate
     module procedure deallocate_l
     module procedure deallocate_d, deallocate_i, deallocate_z
     module procedure deallocate_d2, deallocate_i2, deallocate_z2
     module procedure deallocate_d3, deallocate_i3, deallocate_z3
     module procedure deallocate_d4
     module procedure deallocate_c2
     module procedure deallocate_c3
  end interface


  !---------------------------------------------------------------
  !---------------------------------------------------------------
contains


   subroutine allocate_pl(array,length)
    logical, DIMENSION(:), POINTER :: array
    integer(4) :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
!#         ifdef MEMLOG
!          call writeMemLog
!#         endif
       endif
    endif
  end subroutine allocate_pl



  !---------------------------------------------------------------
  subroutine allocate_pi(array,length)
    integer(4), DIMENSION(:), POINTER :: array
    integer(4) :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_pi

  !---------------------------------------------------------------
  subroutine allocate_pi2(array,row,col)
    integer(4), DIMENSION(:,:), POINTER :: array
    integer(4) :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#	          ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif

  end subroutine allocate_pi2

!---------------------------------------------------------------
  subroutine allocate_pd(array,length)
    real(kind=dp), DIMENSION(:), POINTER :: array
    integer(4) :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef  MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_pd
  !---------------------------------------------------------------
  subroutine allocate_pd2(array,row,col)
    real(dp), DIMENSION(:,:), POINTER :: array
    integer(4) :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#	          ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif

  end subroutine allocate_pd2


  subroutine allocate_pz(array,length)
    complex(kind=dp), DIMENSION(:), POINTER :: array
    integer(4) :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*2*dp    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_pz
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  subroutine allocate_i(array,length)
    integer(4), DIMENSION(:), ALLOCATABLE :: array
    integer(4) :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_i

  subroutine allocate_d(array,length)
    real(kind=dp), DIMENSION(:), ALLOCATABLE :: array
    integer(4) :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (ALLOCATED(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not.ALLOCATED(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef  MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_d

  subroutine allocate_z(array,length)
    complex(kind=dp), DIMENSION(:), ALLOCATABLE :: array
    integer(4) :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (ALLOCATED(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not.ALLOCATED(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*2*dp    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_z
  !---------------------------------------------------------------   
 
  subroutine allocate_i2(array,row,col)
    integer(4), DIMENSION(:,:), ALLOCATABLE :: array
    integer(4) :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#	          ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_i2

  subroutine allocate_d2(array,row,col)
    real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array
    integer(4) :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG	
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_d2

  subroutine allocate_z2(array,row,col)
    complex(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array
    integer(4) :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*2*dp    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_z2
  
  subroutine allocate_c2(array,row,col)
    complex(kind=sp), DIMENSION(:,:), ALLOCATABLE :: array
    integer(4) :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*2*sp    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_c2
  
  subroutine allocate_c3(array,row,col,np)
    complex(kind=sp), DIMENSION(:,:,:), ALLOCATABLE :: array
    integer(4) :: row,col,np,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col,np),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*2*sp    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_c3
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine allocate_i3(array,row,col,dep)
    integer(4), DIMENSION(:,:,:), ALLOCATABLE :: array
    integer(4) :: row,col,dep,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col,dep),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#	          ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_i3

  subroutine allocate_d3(array,row,col,dep)
    real(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array
    integer(4) :: row,col,dep,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col,dep),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp   
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG	
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_d3

  subroutine allocate_z3(array,row,col,dep)
    complex(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array
    integer(4) :: row,col,dep,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col,dep),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*2*dp   
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_z3

  subroutine allocate_d4(array,row,col,dep,qep)
    real(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: array
    integer(4) :: row,col,dep,qep,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col,dep,qep),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp   
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG	
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_d4

  !---------------------------------------------------------------
  !---------------------------------------------------------------
  subroutine allocate_l(array,length)
    logical, DIMENSION(:), ALLOCATABLE :: array
    integer(4) :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4    
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem 
          endif
!#		  ifdef MEMLOG
!          call writeMemLog  
!#		  endif
       endif
    endif
  end subroutine allocate_l
  
  !---------------------------------------------------------------
  !--------------------------------------------------------------- 
  
  subroutine deallocate_pl(array)
    logical, DIMENSION(:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
!#       ifdef MEMLOG
!       call writeMemLog
!#       endif
    else
       write(*,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_pl

  !---------------------------------------------------------------

  subroutine deallocate_pi(array)
    integer(4), DIMENSION(:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pi
  !---------------------------------------------------------------
  subroutine deallocate_pi2(array)
    integer(4), DIMENSION(:,:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pi2
  !---------------------------------------------------------------
  subroutine deallocate_pd(array)
    real(kind=dp), DIMENSION(:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*dp    
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif

    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pd
  !---------------------------------------------------------------
  subroutine deallocate_pd2(array)
    real(dp), DIMENSION(:,:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*dp
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pd2
  !---------------------------------------------------------------
  subroutine deallocate_pz(array)
    complex(kind=dp), DIMENSION(:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*2*dp
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pz
  !---------------------------------------------------------------
  subroutine deallocate_l(array)
    logical, DIMENSION(:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_l
  
  subroutine deallocate_i(array)
    integer(4), DIMENSION(:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_i

  subroutine deallocate_d(array)
    real(kind=dp), DIMENSION(:), allocatable :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*dp   
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif

    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_d

  subroutine deallocate_z(array)
    complex(kind=dp), DIMENSION(:), allocatable :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*2*dp
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_z

  ! ------------------------------------------------------------
  subroutine deallocate_i2(array)
    integer(4), DIMENSION(:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_i2

  subroutine deallocate_d2(array)
    real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*dp
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_d2

  subroutine deallocate_z2(array)
    complex(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*2*dp
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_z2
  
  subroutine deallocate_c2(array)
    complex(kind=sp), DIMENSION(:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*2*sp
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_c2
  
  subroutine deallocate_c3(array)
    complex(kind=sp), DIMENSION(:,:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*2*sp
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_c3
  ! ------------------------------------------------------------

  subroutine deallocate_i3(array)
    integer(4), DIMENSION(:,:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_i3

  subroutine deallocate_d3(array)
    real(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*dp
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_d3

  subroutine deallocate_z3(array)
    complex(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*2*dp
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_z3

  subroutine deallocate_d4(array)
    real(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*dp
       deallocate(array)
!#		ifdef MEMLOG
!       call writeMemLog  
!#		endif
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_d4

  ! ------------------------------------------------------------
  subroutine resetMemLog
     alloc_mem=0
     peak_mem=0
  end subroutine resetMemLog
  ! ------------------------------------------------------------
  subroutine getMem(mem)
    integer(8) :: mem

    mem = alloc_mem
  end subroutine getMem
  ! ------------------------------------------------------------
  subroutine writeMemInfo(iofile)

    integer iofile
    character(3) :: str
    integer :: dec

    call memstr(alloc_mem,dec,str)
    write(iofile,'(A26,F8.2,A3)') 'current memory allocated: ',alloc_mem*1.0/dec,str

  end subroutine writeMemInfo
  ! ------------------------------------------------------------
  subroutine writePeakInfo(iofile)

    integer iofile
    character(3) :: str
    integer :: dec

    call memstr(peak_mem,dec,str)
    write(iofile,'(A26,F8.2,A3)') 'peak memory allocated: ',peak_mem*1.0/dec,str

  end subroutine writePeakInfo

  ! ------------------------------------------------------------
  subroutine openMemLog(iofile)
    integer iofile,err

    if(iofile.ne.6) then
       open(iofile,file='memory.log',iostat=err)
       if (err.ne.0) then
          write(*,*) 'Cannot open memory log-file'
          stop
       endif
    endif
    iolog=iofile

  end subroutine openMemLog

  ! ------------------------------------------------------------
  subroutine writeMemLog
    call writeMemInfo(iolog)
  end subroutine writeMemLog

  ! ------------------------------------------------------------
  subroutine closeMemLog
    call writePeakInfo(iolog)
    close(iolog)
  end subroutine closeMemLog

  ! ------------------------------------------------------------
  subroutine memstr(mem,dec,str)
 
    character(3) :: str
    integer(8) :: mem
    integer :: dec

    if(mem.lt.1000) then
       str=' bt'; dec=1
       return      
    endif

    if(mem.lt.1e7) then
       str=' kb'; dec=1000
       return      
    endif

    if(mem.ge.1e7) then
       str=' Mb'; dec=1000000
       return      
    endif

  end subroutine memstr

end module ln_allocation
