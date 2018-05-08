!--------------------------------------------------------------------------------------------------!                              
! DFTB+XT: DFTB+ eXTended version for model and atomistic quantum transport at nanoscale.          !
!                                                                                                  !
! Copyright (C) 2017-2018 DFTB+ developers group.                                                  !
! Copyright (C) 2018 Dmitry A. Ryndyk.                                                             !
!                                                                                                  !
! GNU Lesser General Public License version 3 or (at your option) any later version.               !
! See the LICENSE file for terms of usage and distribution.                                        !
!--------------------------------------------------------------------------------------------------!
! This file is part of the TraNaS library for quantum transport at nanoscale.                      !
!                                                                                                  !
! Developer: Dmitry A. Ryndyk.                                                                     !
!                                                                                                  !
! Based on the LibNEGF library developed by                                                        !
! Alessandro Pecchia, Gabriele Penazzi, Luca Latessa, Aldo Di Carlo.                               !
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

!> The module implements the classes for many-body interactions.

!> The module implements an abstract class to interface different
!! many-body interactions. 

module interactions

  use globals, only : LST
  use ln_precision, only : dp
  use mat_def, only : z_dns,z_csr
  use ln_structure, only : TStruct_info

  implicit none
  private

  public :: interaction, Tmbngf

  !-----------------------------------------------------------------------------

  type, abstract :: Interaction

     character(len=LST) :: descriptor
     !> Maximum number of SCBA iterations. 
     !! corresponds to no iterations (self energy is not calculated)
     integer :: scba_niter = 0
     real(kind=dp) :: Mixing = 0.5
     !> Keep track of SCBA iteration 
     integer :: scba_iter = 0
     !> SCBA Tolerance
     real(dp) :: scba_tol = 1.0d-7
     !> Number of energy points from integration grid
     !integer :: en_npoints = 0
     !> Buffer for Gr

     !> Local diagonal representation of retarded self energy
!     complex(dp), allocatable, dimension(:) :: sigma_r

     !> System partitioning (as in TNEGF)
     type(TStruct_info) :: struct

   contains

     procedure(abst_add_sigma_r), deferred :: add_sigma_r
     procedure(abst_get_sigma_n), deferred :: get_sigma_n
     procedure(abst_set_Gr), deferred :: set_Gr
     procedure(abst_set_Gn), deferred :: set_Gn

  end type Interaction

  abstract interface

     !> This interface should append
     !! the retarded self energy to ESH
     subroutine abst_add_sigma_r(this, esh)
       import :: interaction
       import :: z_dns
       class(interaction) :: this
       type(z_dns), dimension(:,:), allocatable, intent(inout) :: esh
     end subroutine abst_add_sigma_r

     !> Returns the lesser (n) Self Energy in block format
     !! @param [in] this: calling instance
     !! @param [in] struct: system structure
     !! @param [inout] blk_sigma_n: block dense sigma_n
     !! @param [in] ie: index of energy point
     subroutine abst_get_sigma_n(this, blk_sigma_n, en_index)
       import :: interaction
       import :: z_dns
       class(interaction) :: this
       type(z_dns), dimension(:,:), allocatable, intent(inout) :: blk_sigma_n
       integer, intent(in) :: en_index
     end subroutine abst_get_sigma_n

     !> Give the Gr at given energy point to the interaction
     subroutine abst_set_Gr(this, Gr, en_index, mix)
       import :: interaction
       import :: z_dns
       import :: dp
       class(interaction) :: this
       type(z_dns), dimension(:,:), allocatable, intent(in) :: Gr
       integer :: en_index
       real(dp) :: mix
     end subroutine abst_set_Gr

     !> Give the Gn at given energy point to the interaction
     subroutine abst_set_Gn(this, Gn, en_index, mix)
       import :: interaction
       import :: z_dns
       import :: dp
       class(interaction) :: this
       type(z_dns), dimension(:,:), allocatable, intent(in) :: Gn
       integer :: en_index
       real(dp) :: mix
     end subroutine abst_set_Gn

  end interface

  !contains

  !> Initialize information needed for buffering G on memory or disk
  !  Now it only pass the number of energy grid points but it 
  !  could turn in something more complicated (e.g. an energy path object)
!!$    subroutine init_Gbuffer(this, en_npoints)
!!$      class(interaction) :: this
!!$      integer, intent(in) :: en_npoints
!!$
!!$      this%en_npoints = en_npoints
!!$      
!!$      
!!$    end subroutine init_Gbuffer

  !-----------------------------------------------------------------------------
  !DAR begin - class Tmbngf
  !-----------------------------------------------------------------------------
  type :: Tmbngf

    !> System partitioning (as in Tnegf%str)
    type(TStruct_info) :: str

    !> Block retarded self-energy for HF approximation
    type(z_DNS), dimension(:,:), allocatable :: SelfEnergyR_HF

    !> Retarded Green functions and self-energies
    type(z_CSR), dimension(:), allocatable :: GreenFunctionR
    type(z_CSR), dimension(:), allocatable :: PolarizationOperatorR

    !> Lesser Green functions and self-energies

    !> SCC Tolerance
    real(dp) :: scc_tol = 1.0d-7

    logical :: tHartreeFock = .false.
    logical :: tRPA = .false.

  contains

    procedure :: add_SelfEnergyR_HF
    procedure :: get_SelfEnergyR_HF
    procedure :: get_SelfEnergyR_RHF

  end type Tmbngf

contains

  !-----------------------------------------------------------------------------  

  !> Append the retarded self-energy to ESH
  !! Hartree-Fock approximation
  subroutine add_SelfEnergyR_HF(mbngf,ESH)

    class(Tmbngf) :: mbngf
    type(z_DNS), dimension(:,:), allocatable, intent(inout) :: ESH
    integer :: n, nbl, ii

    nbl = mbngf%str%num_PLs

    do n=1,nbl
       ESH(n,n)%val = ESH(n,n)%val - mbngf%SelfEnergyR_HF(n,n)%val
    end do

    do n=2,nbl
       ESH(n-1,n)%val = ESH(n-1,n)%val - mbngf%SelfEnergyR_HF(n-1,n)%val
       ESH(n,n-1)%val = ESH(n,n-1)%val - mbngf%SelfEnergyR_HF(n,n-1)%val
    end do

  end subroutine add_SelfEnergyR_HF

  !-----------------------------------------------------------------------------  

  !> Calculates the retarded self-energy
  !! Hartree-Fock approximation
  subroutine get_SelfEnergyR_HF(mbngf,rho_dense,U)

    class(Tmbngf) :: mbngf
    type(z_DNS), intent(in) :: rho_dense
    real(kind=dp), dimension(:,:), intent(in) :: U
    integer :: n,m, nbl, ii,jj,kk
    
!DAR! Add later allocation for non-diagonal elements.
!DAR! At the moment only interactions inside every PL are taken into account. 

    !debug begin
    !do ii=1,mbngf%str%central_dim
    !   write(*,"(10000F9.3)")U(ii,1:mbngf%str%central_dim)
    !end do
    !debug end
    
    nbl = mbngf%str%num_PLs   
    do n=1,nbl
    !do m=1,nbl   
       associate(pl_start1=>mbngf%str%mat_PL_start(n),pl_end1=>mbngf%str%mat_PL_end(n), &
                 pl_start2=>mbngf%str%mat_PL_start(n),pl_end2=>mbngf%str%mat_PL_end(n))
         do ii = 1, pl_end1 - pl_start1 + 1
            forall(jj = 1:pl_end2 - pl_start2 + 1)
               mbngf%SelfEnergyR_HF(n,n)%val(ii,jj)=-U(ii,jj)*rho_dense%val(pl_start1+ii-1,pl_start1+jj-1)       
            end forall
            do kk = 1, pl_end1 - pl_start1 + 1
               mbngf%SelfEnergyR_HF(n,n)%val(ii,ii)=mbngf%SelfEnergyR_HF(n,n)%val(ii,ii)+ &
                    U(ii,kk)*rho_dense%val(pl_start1+kk-1,pl_start1+kk-1)          
            end do
         end do
       end associate
    !end do
    end do

  end subroutine get_SelfEnergyR_HF

  !-----------------------------------------------------------------------------  

  !> Calculates the retarded self-energy
  !! RESTRICTED Hartree-Fock approximation for SPINLESS Hamiltonian
  subroutine get_SelfEnergyR_RHF(mbngf,rho_dense,U)

    class(Tmbngf) :: mbngf
    type(z_DNS), intent(in) :: rho_dense
    real(kind=dp), dimension(:,:), intent(in) :: U
    integer :: n,m, nbl, ii,jj,kk
    
!DAR! Add later allocation for non-diagonal elements.
!DAR! At the moment only interactions inside every PL are taken into account. 

    !debug begin
    !do ii=1,mbngf%str%central_dim
    !   write(*,"(10000F9.3)")U(ii,1:mbngf%str%central_dim)
    !end do
    !debug end
    
    nbl = mbngf%str%num_PLs   
    do n=1,nbl
    !do m=1,nbl   
       associate(pl_start1=>mbngf%str%mat_PL_start(n),pl_end1=>mbngf%str%mat_PL_end(n), &
                 pl_start2=>mbngf%str%mat_PL_start(n),pl_end2=>mbngf%str%mat_PL_end(n))
         do ii = 1, pl_end1 - pl_start1 + 1
            forall(jj = 1:pl_end2 - pl_start2 + 1)
               mbngf%SelfEnergyR_HF(n,n)%val(ii,jj)=-U(ii,jj)*rho_dense%val(pl_start1+ii-1,pl_start1+jj-1)
            end forall
            do kk = 1, pl_end1 - pl_start1 + 1
               mbngf%SelfEnergyR_HF(n,n)%val(ii,ii)=mbngf%SelfEnergyR_HF(n,n)%val(ii,ii)+ &
                    2.*U(ii,kk)*rho_dense%val(pl_start1+kk-1,pl_start1+kk-1)
            end do
         end do
       end associate
    !end do
    end do
    
  end subroutine get_SelfEnergyR_RHF

  !-----------------------------------------------------------------------------
  !DAR end
  !-----------------------------------------------------------------------------

end module interactions
