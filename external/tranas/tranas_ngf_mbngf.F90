!--------------------------------------------------------------------------------------------------!
! DFTB+XT open software package for quantum nanoscale modeling                                     !
! Copyright (C) 2017-2018 DFTB+ developers group                                                   !
! Copyright (C) 2018 Dmitry A. Ryndyk                                                              !
!--------------------------------------------------------------------------------------------------!
! GNU Lesser General Public License version 3 or (at your option) any later version.               !
! See the LICENSE file for terms of usage and distribution.                                        !
!--------------------------------------------------------------------------------------------------!
! This file is a part of the TraNaS library for quantum transport at nanoscale.                    !
! Developer: Dmitry A. Ryndyk.                                                                     !
! Partially based on the LibNEGF library developed by                                              !
! Alessandro Pecchia, Gabriele Penazzi, Luca Latessa, Aldo Di Carlo.                               !
!--------------------------------------------------------------------------------------------------!
  
module tranas_ngf_mbngf

  use mpi_globals, only : id, numprocs, id0  
  use libmpifx_module
  use ln_precision
  use sparsekit_drv  
  use mat_def

  use tranas_types_main
  use tranas_types_mbngf, only : TMBNGF
  use tranas_ngf_integrations

  implicit none
  private

  public :: mbngfInit
  public :: mbngfCompute
  public :: mbngfDestroy
  
!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

  !------------------------------------------------------------------------------------------------!
  !> Allocation of the many-body self-energies.
  !------------------------------------------------------------------------------------------------!  
  subroutine mbngfInit(tranas)

    type(TTraNaS), target :: tranas
    type(Tnegf), pointer :: negf

    integer :: n,nbl,ierr

    negf => tranas%negf
    tranas%ngf%mbngf%str=negf%str
  
    ! Allocation of self-energies
    
    if(negf%mbngf%tHartreeFock) then
      nbl = negf%str%num_PLs
      allocate(negf%mbngf%SelfEnergyR_HF(nbl,nbl),stat=ierr)
      if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate negf%mbngf%SelfEnergyR_HF(nbl,nbl)'
      do n=1,nbl
        associate(pl_number => negf%str%mat_PL_end(n)-negf%str%mat_PL_start(n)+1)
          allocate(negf%mbngf%SelfEnergyR_HF(n,n)%val(pl_number,pl_number))
        end associate
        negf%mbngf%SelfEnergyR_HF(n,n)%val = 0.0_dp
      end do
      !DAR! Add later the allocation for non-diagonal block elements.
      !DAR! At the moment only interactions inside every PL are taken into account.     
    
      ! Computation of the density matrix and the HF self-energy
      !negf%readOldSGF=2
      negf%iteration=1
      negf%DorE = 'D'
      negf%outer=0
    end if

    if(negf%mbngf%tRPA) then
    
    end if
    
  end subroutine mbngfInit
    
  !------------------------------------------------------------------------------------------------!
  !> Calculation of the many-body self-energies.
  !------------------------------------------------------------------------------------------------!  
  subroutine mbngfCompute(tranas)

    type(TTraNaS), target :: tranas
    type(Tnegf), pointer :: negf
    type(TMBNGF), pointer :: mbngf
    type(z_DNS) :: rho_dense,rho_dense_previous,sigma_dense

    integer :: n,m,ii,jj,nbl,ierr,scba_iter
    real(dp) :: scba_error

    negf => tranas%negf
    mbngf => tranas%ngf%mbngf

    scba_error = 0.0_dp
  
    if (id0.and.negf%verbose.gt.50) then
      write(*,*)   
      write(*,"('>>> The MBNGF calculation is started.')")
      if (negf%MaxIter==0) write(*,"('>> Calculation is skipped because MaxNumIter = 0.')")
      if (negf%MaxIter==1) write(*,"('>> First-order calculation is started.')")
      if (negf%MaxIter>1) then
        write(*,"('>> Self-Consistent cycle is started.')")
        write(*,"('Maximum number of iterations  =   ',I0)")negf%MaxIter
      end if   
    end if

    call allocation
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Here the SC cycle is started. !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    do scba_iter = 1,negf%MaxIter
       
      if (id0.and.negf%verbose.gt.50) write(*,"('> SC iteration',I6)")scba_iter

      !------------------------------------------------------------------------
      ! Density calculation
      !------------------------------------------------------------------------

      if (negf%mbngf%tHartreeFock) call density

      !------------------------------------------------------------------------
      ! Self-energy calculation
      !------------------------------------------------------------------------

      call sigma

      !------------------------------------------------------------------------
      ! Convergency check
      !------------------------------------------------------------------------

      call error
       
      if (id0.and.negf%verbose.gt.50) write(*,"('Error at SC iteration ',I6,' : ',E12.6)") scba_iter, scba_error
      if (scba_error .lt. negf%Tolerance) then !negf%inter%scba_tol) then 
        if (id0.and.negf%verbose.gt.50) write(*,"('SC exit succesfully after ',I6,' iterations')")scba_iter
        exit
      end if

      call new_cycle

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Here the SC cycle is finished. !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (scba_error.gt.negf%Tolerance) then 
      if (id0.and.negf%verbose.gt.0) write(*,"('WARNING : SC exit with error ',E12.6, &
           & ' above reference tolerance ',E12.6)") scba_error, negf%Tolerance
    end if

    call deallocation

    if (id0.and.negf%verbose.gt.50.and.(negf%MaxIter.gt.0)) then
      !write(*,*)   
      write(*,"('>> Self-Consistent cycle is finished.')")
    end if
    
    if (id0.and.negf%verbose.gt.50) then
      !write(*,*)
      write(*,"('>>> The MBNGF calculation is ended.')")
    end if
  
  CONTAINS  

  !------------------------------------------------------------------------------------------------!

  subroutine allocation

    nbl = negf%str%num_PLs

    if (negf%mbngf%tHartreeFock) then
      allocate(rho_dense%val(negf%str%total_dim,negf%str%total_dim),stat=ierr)
      if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate rho_dense%val(negf%str%total_dim,negf%str%total_dim)'
      allocate(rho_dense_previous%val(negf%str%total_dim,negf%str%total_dim))
      rho_dense%val = 0.0_dp
      rho_dense%nrow = negf%str%total_dim
      rho_dense%ncol = negf%str%total_dim
      rho_dense_previous%nrow = negf%str%total_dim
      rho_dense_previous%ncol = negf%str%total_dim
      rho_dense_previous = rho_dense
    end if

  end subroutine allocation

  !------------------------------------------------------------------------------------------------!

  subroutine density
  
    if (id0.and.negf%verbose.gt.50) write(*,"('Density is started')")
       
    rho_dense_previous=rho_dense  

    !! Did anyone passed externally allocated DM? If not, create it
    call create_DM(negf)
    
    ! Reference contact for contour/real axis separation
    call set_ref_cont(negf)
   
    !Decide what to do with surface GFs.
    !sets readOldSGF: if it is 0 or 1 it is left so 
    if (negf%readOldSGF.eq.2) then
      if(negf%iteration.eq.1) then        
        negf%readOldSGF=2  ! compute and save SGF on files
      else
        negf%readOldSGF=0  ! read from files
      endif
    endif

    if (negf%Np_n(1)+negf%Np_n(2)+negf%n_poles.gt.0) then
      call contour_int_def(negf)
      call contour_int(negf)
    endif

    if (negf%Np_real(1).gt.0) then
      call real_axis_int_def(negf)
      call real_axis_int(tranas)
    endif
       
    call mpifx_allreduceip(negf%mpicomm, negf%rho%nzval, MPI_SUM)

    call csr2dns(negf%rho,rho_dense)

    rho_dense%val=negf%Mixing*rho_dense%val+(1-negf%Mixing)*rho_dense_previous%val    

    !debug begin
    if (id0.and.negf%verbose.gt.70) then
      print *, 'real(rho_dense)'
      do ii=1,negf%str%central_dim
        print *, real(rho_dense%val(ii,1:negf%str%central_dim))
      end do
      !do ii=1,negf%str%central_dim 
      !  write(*,"(10000F12.8)")negf%U(ii,1:negf%str%central_dim)
      !end do
    end if  
    !debug end

    if (id0.and.negf%verbose.gt.50) write(*,"('Density is calculated')")

  end subroutine density

  !------------------------------------------------------------------------------------------------!

  subroutine sigma

    if (negf%mbngf%tHartreeFock) then
 
      if (id0.and.negf%verbose.gt.50) write(*,"('Calculation of HF self-energy is started')")

      call negf%mbngf%get_SelfEnergyR_HF(rho_dense,negf%U)
      !call negf%mbngf%get_SelfEnergyR_RHF(rho_dense,negf%U)   
      !#!DAR Change to negf%mbngf%get_SelfEnergyR_HF(negf%rho,negf%U)
      deallocate(negf%rho) ! It is important to deallocate negf%rho.

      !debug begin
      if (id0.and.negf%verbose.gt.70) then
        allocate(sigma_dense%val(negf%str%total_dim,negf%str%total_dim))
        sigma_dense%nrow=negf%str%total_dim
        sigma_dense%ncol=negf%str%total_dim
        print *, 'real(SelfEnergyR_HF)'
        sigma_dense%val=0.0_dp                
        do n=1,nbl
        !do m=1,nbl   
          associate(pl_start1=>negf%str%mat_PL_start(n),pl_end1=>negf%str%mat_PL_end(n), &
               pl_start2=>negf%str%mat_PL_start(n),pl_end2=>negf%str%mat_PL_end(n))
            forall(ii = 1:pl_end1 - pl_start1 + 1)
            forall(jj = 1:pl_end2 - pl_start2 + 1)
              sigma_dense%val(pl_start1+ii-1,pl_start2+jj-1) = negf%mbngf%SelfEnergyR_HF(n,n)%val(ii,jj)
            end forall
            end forall
          end associate
          !end do
        end do
        do ii=1,negf%str%central_dim
          print *, real(sigma_dense%val(ii,1:negf%str%central_dim))
        end do
        deallocate(sigma_dense%val)
      end if
      !debug end

      if (id0.and.negf%verbose.gt.50) write(*,"('HF self-energy is calculated')")

    end if

    if (tranas%input%tPhotons) then
 
      if (id0.and.negf%verbose.gt.50) write(*,"('Calculation of PHOTON Born self-energy is started.')")

      call transmission_int_def(negf)
      call integrationsSelfEnergies(tranas)



      

      if (id0.and.negf%verbose.gt.50) write(*,"('Photon Born self-energy is calculated.')")

    end if
      
  end subroutine sigma

  !------------------------------------------------------------------------------------------------!

  subroutine error
    
    if (negf%mbngf%tHartreeFock) scba_error = maxval(abs(rho_dense%val - rho_dense_previous%val))/maxval(abs(rho_dense%val))

  end subroutine error

  !------------------------------------------------------------------------------------------------!

  subroutine new_cycle

    if (negf%mbngf%tHartreeFock) call destroy(rho_dense_previous)

  end subroutine new_cycle

  !------------------------------------------------------------------------------------------------!

  subroutine deallocation
  
    if (negf%mbngf%tHartreeFock) deallocate(rho_dense%val)

  end subroutine deallocation

  !------------------------------------------------------------------------------------------------!
  
  end subroutine mbngfCompute

  !------------------------------------------------------------------------------------------------!
  !> Deallocation of the many-body self-energies.
  !------------------------------------------------------------------------------------------------!
  subroutine mbngfDestroy(tranas)

    use tranas_ngf_iterative, only : destroy_blk

    type(TTraNaS) :: tranas
    type(Tnegf) :: negf

    negf = tranas%negf
    
    if (negf%mbngf%tHartreeFock) then
      call destroy_blk(negf%mbngf%SelfEnergyR_HF)
      deallocate(negf%mbngf%SelfEnergyR_HF)
    end if

    tranas%negf = negf
    
  end subroutine mbngfDestroy

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!  



!--------------------------------------------------------------------------------------------------!
!> Create density matrix (copy from tranas.F90).
!--------------------------------------------------------------------------------------------------!

subroutine create_DM(negf)

  type(Tnegf) :: negf   
  
  if (negf%intDM) then
    if (.not.associated(negf%rho)) allocate(negf%rho)
    if (.not.associated(negf%rho_eps)) allocate(negf%rho_eps)
  end if  

end subroutine create_DM

!--------------------------------------------------------------------------------------------------!
!> Destroy density matrix (copy from tranas.F90).
!--------------------------------------------------------------------------------------------------!
  
subroutine destroy_DM(negf)

  type(Tnegf) :: negf   

  if (negf%intDM) then
      
    if (associated(negf%rho)) then
      if (allocated(negf%rho%nzval)) then
        !print*,'(destroy) deallocate negf%rho',%LOC(negf%rho%nzval)
        call destroy(negf%rho) 
      end if
      deallocate(negf%rho)
      nullify(negf%rho)
    end if
      
    if (associated(negf%rho_eps)) then
      if (allocated(negf%rho_eps%nzval)) then
        !print*,'(destroy) deallocate negf%rho_eps',%LOC(negf%rho_eps%nzval)
        call destroy(negf%rho_eps) 
      end if
      deallocate(negf%rho_eps)
      nullify(negf%rho_eps) 
    end if

  end if  

end subroutine destroy_DM

!--------------------------------------------------------------------------------------------------!
!> Destroy matrices created runtime (copy from tranas.F90).
!--------------------------------------------------------------------------------------------------!  

subroutine destroy_matrices(negf)
    
  type(Tnegf) :: negf   
  integer :: i

  do i=1,negf%str%num_conts
    if (allocated(negf%HC(i)%val)) call destroy(negf%HC(i))
    if (allocated(negf%SC(i)%val)) call destroy(negf%SC(i))
    if (allocated(negf%HMC(i)%val)) call destroy(negf%HMC(i))
    if (allocated(negf%SMC(i)%val)) call destroy(negf%SMC(i))
  end do
       
end subroutine destroy_matrices

!--------------------------------------------------------------------------------------------------!
!> Sets the Reference contact for non-eq calculations (copy from tranas.F90).
!! 
!! The behaviour depends on how negf%min_or_max has been set.
!!
!! min_or_max = 0 : refcont is chosen at the minimum   mu
!! min_or_max = 1 : refcont is chosen at the maximum   mu 
!--------------------------------------------------------------------------------------------------!

subroutine set_ref_cont(negf)

  type(TNegf) :: negf

  integer :: nc_vec(1), ncont

  ncont = negf%str%num_conts

  if (ncont > 0) then
    if (negf%min_or_max .eq. 0) then
      negf%muref = minval(negf%mu(1:ncont))
      nc_vec = minloc(negf%mu(1:ncont))  
    else
      negf%muref = maxval(negf%mu(1:ncont))
      nc_vec = maxloc(negf%mu(1:ncont))
    end if
    negf%refcont = nc_vec(1)
  else
    negf%muref = negf%mu(1)
    negf%refcont = 1  
  end if
     
end subroutine set_ref_cont

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
  
end module tranas_ngf_mbngf
  
