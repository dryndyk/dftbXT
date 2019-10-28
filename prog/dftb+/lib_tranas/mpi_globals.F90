!--------------------------------------------------------------------------------------------------!
!  DFTB+XT open software package for quantum nanoscale modeling                                    !
!  Copyright (C) 2018 Dmitry A. Ryndyk                                                             !
!--------------------------------------------------------------------------------------------------!
!  GNU Lesser General Public License version 3 or (at your option) any later version.              !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
!  This file is part of the TraNaS library for quantum transport at nanoscale.                     !
!  Developer: Dmitry A. Ryndyk.                                                                    !
!  Based on the LibNEGF library developed by                                                       !
!  Alessandro Pecchia, Gabriele Penazzi, Luca Latessa, Aldo Di Carlo.                              !
!--------------------------------------------------------------------------------------------------!

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

#:include 'common.fypp'

module mpi_globals

#:if WITH_MPI == 1
  
  use libmpifx_module
  
  INTEGER, SAVE ::  mpi_comm
  INTEGER, SAVE ::  numprocs
  INTEGER, SAVE ::  id
  LOGICAL, SAVE ::  id0

  contains
    
  subroutine negf_mpi_init(comm)
    
    type(mpifx_comm) :: comm
        
    id = comm%rank
    numprocs = comm%size
      
    id0 = .false.
    if (id.eq.0) id0 = .true.  

    !print*, 'INIT MPI-NEGF ON',numprocs,'NODES' 
    !print*, 'CPU',id,'READY'
    !print*, 'PRINTING CPU:',id0
      
  end subroutine negf_mpi_init

#:else
  
  logical, parameter ::  id0 = .true.
  integer, parameter :: id = 0, numprocs = 1
         
#:endif

end module mpi_globals
