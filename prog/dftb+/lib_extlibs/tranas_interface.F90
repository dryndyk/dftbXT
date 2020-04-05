!--------------------------------------------------------------------------------------------------!
! DFTB+XT open software package for quantum nanoscale modeling (TraNaS OpenSuite)                  !
! Copyright (C) 2018-2019 Dmitry A. Ryndyk                                                         !
! DFTB+: general package for performing fast atomistic simulations                                 !
! Copyright (C) 2017-2019 DFTB+ developers group                                                   !
!--------------------------------------------------------------------------------------------------!
! GNU Lesser General Public License version 3 or (at your option) any later version.               !
! See the LICENSE file for terms of usage and distribution.                                        !
!--------------------------------------------------------------------------------------------------!
! This file is a part of the DFTB+XT/TraNaS interface.                                             !
! Developer: Dmitry A. Ryndyk.                                                                     !
! Based on the DFTB+/LibNEGF interface.                                                            !
!--------------------------------------------------------------------------------------------------!
! TraNaS is the library for quantum transport at nanoscale.                                        !
! Developer: Dmitry A. Ryndyk.                                                                     !
! Based on the LibNEGF library developed by                                                        !
! Alessandro Pecchia, Gabriele Penazzi, Luca Latessa, Aldo Di Carlo.                               !
!--------------------------------------------------------------------------------------------------!

module tranas_interface

  use dftbp_accuracy, only : mc, lc
  !use dftbp_constants
  use dftbp_matconv
  use dftbp_sparse2dense
  use dftbp_densedescr
  use dftbp_commontypes, only : TOrbitals
  use dftbp_mpifx
  use dftbp_formatout
  use dftbp_globalenv
  use dftbp_message
  use dftbp_elecsolvertypes, only : electronicSolverTypes
  
  use ln_precision
  use ln_constants
  use ln_allocation
  use tranas_vars   !, only : TGDFTBstructure, TGDFTBTunDos, &
                    !       & TGDFTBGreenDensInfo, TTransPar, Telph
  use ln_structure
  use mpi_globals
  use mat_def 
  use lib_param    !, only : Tnegf, set_defaults,  & 
                   !     set_elph_dephasing, set_elph_block_dephasing, &
                   !     set_elph_s_dephasing, destroy_elph_model
  use tranas_types_main                
  use tranas
  use ln_extract
  use libmpifx_module

  use dftbp_taggedoutput                                                         
  
  implicit none
  private

  Type(Tnegf), target, save, public :: negf

  public :: negf_init      ! general library initializations
  public :: negf_init_nogeom                                               
  public :: negf_init_str  ! passing structure parameters
  public :: negf_init_elph                                                 
  public :: negf_destroy
  
  ! wrapped functions passing dftb matrices. Needed for parallel
  public :: calcdensity_green
  public :: calcEdensity_green
  public :: calcPDOS_green
  public :: calc_current
  public :: local_currents
  public :: tranasNoGeom                                                    !DAR
 
  !!$ public :: negf_init_elph, negf_destroy_elph
  
  ! interface csr matrices. The pattering must be predefined 
  ! using negf_init_csr 
  public :: negf_init_csr
  type(z_CSR), SAVE :: csrHam, csrOver

  ! non wrapped direct calls  
  private :: negf_density, negf_current, negf_ldos

  contains

  !------------------------------------------------------------------------------------------------!
  !> Initialize QT environment and variables (public)
  !------------------------------------------------------------------------------------------------!
  subroutine negf_init(structure, transpar, greendens, tundos, mpicomm, initinfo)
          
    Type(TTranspar), intent(IN) :: transpar
    Type(TGDFTBstructure), intent(IN) :: structure
    Type(TGDFTBGreenDensInfo), intent(IN) :: greendens
    Type(TGDFTBTunDos), intent(IN) :: tundos
    Type(mpifx_comm), intent(in) :: mpicomm
    logical, intent(OUT) :: initinfo
        
    ! local variables
    real(dp), dimension(:), allocatable :: pot, eFermi
    integer :: i,error, l, ncont, nc_vec(1), j, nldos
    integer, dimension(:), allocatable :: sizes  
    ! string needed to hold processor name
    character(:), allocatable :: hostname
    type(lnParams) :: params
    
    error = 0 
    initinfo = .true.       

    call negf_mpi_init(mpicomm)
     
    if (transpar%defined) then
      ncont = transpar%ncont 
    else
      ncont = 0
    endif
    
    ! ------------------------------------------------------------------------------
    ! check that GS and contacts are created  (one day this will be automatic)
    ! ------------------------------------------------------------------------------
    if(transpar%defined .and. ncont.gt.0) then
                                
      open(199,file='GS/test',IOSTAT=error)
      if (error.ne.0) then 
        call mpifx_get_processor_name(hostname)
        write(*,*) 'ERROR: please create a directory called "GS"'
        initinfo = .false.; stop  
      end if   
      close(199)
      open(199,file='contacts/test',IOSTAT=error)
      if (error.ne.0) then 
        call mpifx_get_processor_name(hostname)
        write(*,*) 'ERROR: please create a directory called "contacts"'
        initinfo = .false.; stop  
      end if         
      close(199)
    
    end if
    ! ------------------------------------------------------------------------------
    !! Set defaults and fill up the parameter structure with them
    call init_negf(negf)
    call get_params(negf, params)
    
    ! ------------------------------------------------------------------------------
    !This must be different for different initialisations, to be separated
    !Higher between transport and greendens is taken, temporary
    !DAR begin - params%verbose
    if (tundos%defined .and. greendens%defined) then                                          
       if (tundos%verbose.gt.greendens%verbose) then 
          params%verbose = tundos%verbose
       else
          params%verbose = greendens%verbose
       endif
    else
       if (tundos%defined) params%verbose = tundos%verbose
       if (greendens%defined) params%verbose = greendens%verbose
    end if
    !DAR end 
    ! ------------------------------------------------------------------------------
    ! This parameter is used to set the averall drop threshold in libnegf
    ! It affects especially transmission that is not accurate more than 
    ! this value.
    call set_drop(1.d-20)


    ! ------------------------------------------------------------------------------
    ! Assign spin degenracy and check consistency between different input blocks
    if (tundos%defined .and. greendens%defined) then
      if (tundos%gSpin .ne. greendens%gSpin) then
        write(*,*) "ERROR: spin degeneracy is not consistent between different input blocks"
        initinfo = .false.; return 
      else
        params%g_spin = real(tundos%gSpin) ! Spin degeneracy
      end if
    else if (tundos%defined) then
        params%g_spin = real(tundos%gSpin) ! Spin degeneracy
    else if (greendens%defined) then
        params%g_spin = real(greendens%gSpin) ! Spin degeneracy
    end if



    ! ------------------------------------------------------------------------------
    !            SETTING ELECTOCHEMICAL POTENTIALS INCLUDING BUILT-IN 
    ! ------------------------------------------------------------------------------
    ! Fermi level is given by the contacts. If no contacts => no transport,
    ! Then Fermi is defined by the Green solver
    if (transpar%defined) then
      pot = transpar%contacts(1:ncont)%potential
      eFermi = transpar%contacts(1:ncont)%eFermi(1) 
      do i = 1,ncont
        ! Built-in potential to equilibrate Fermi levels
        pot(i) = pot(i) + eFermi(i) - minval(eFermi(1:ncont))

        ! Setting contact temperatures
        if (transpar%contacts(i)%kbT .ge. 0.0_dp) then
          params%kbT(i) = transpar%contacts(i)%kbT
        else
          params%kbT(i) = structure%tempElec
          ! Make sure low temperatures (< 10K) converted to 0.
          ! This avoid numerical problems with contour integration
          if (params%kbT(i) < 3.0e-5_dp) then
            params%kbT(i) = 0.0_dp
          end if
        end if

        ! set parameters for wide band approximations
        params%FictCont(i) = transpar%contacts(i)%wideBand
        params%contact_DOS(i) = transpar%contacts(i)%wideBandDOS
      
        if (id0.and.transpar%verbose.gt.50) then
          write(*,*) '(negf_init) CONTACT INFO #',i
            if (params%FictCont(i)) then
              write(*,*) 'FICTICIOUS CONTACT '
              write(*,*) 'DOS: ', params%contact_DOS(i)
            end if
          write(*,*) 'Temperature: ', params%kbT(i)
          write(*,*) 'Potential (with built-in): ', pot(i)
          write(*,*) 'eFermi: ', eFermi(i) 
          write(*,*) 
        endif 

      enddo

      ! Define electrochemical potentials
      params%mu(1:ncont) = eFermi(1:ncont) - pot(1:ncont)
      if (id0.and.transpar%verbose.gt.50) write(*,*) 'Electro-chemical potentials: ', params%mu(1:ncont)
      deallocate(pot)

    else !transpar not defined 
      params%mu(1) = greendens%oneFermi(1) 
    end if


    ! ------------------------------------------------------------------------------
    !                  SETTING COUNTOUR INTEGRATION PARAMETERS
    ! ------------------------------------------------------------------------------
    if (greendens%defined) then
      params%Ec = greendens%enLow           ! lowest energy
      params%Np_n(1:2) = greendens%nP(1:2)  ! contour npoints
      params%Np_real = greendens%nP(3)      ! real axis points
      params%n_kt = greendens%nkt           ! n*kT for Fermi

      !Read G.F. from very first iter
      if (greendens%readSGF .and. .not.greendens%saveSGF) params%readOldSGF=0  
      !compute G.F. at every iteration 
      if (.not.greendens%readSGF .and. .not.greendens%saveSGF) params%readOldSGF=1  
      !Default Write on first iter
      if (.not.greendens%readSGF .and. greendens%saveSGF) params%readOldSGF=2  

      if(any(params%kbT.gt.0) .and. greendens%nPoles.eq.0) then
         STOP 'ERROR: Number of Poles = 0 but T > 0' 
      else
         params%n_poles = greendens%nPoles
      end if
      if(all(params%kbT.eq.0)) then 
        params%n_poles = 0
      end if

    end if

    ! ------------------------------------------------------------------------------
    !! Setting the delta: priority on Green Solver, if present
    !! dos_delta is used by libnegf to smoothen T(E) and DOS(E)
    !! and is currently set in transmission
    params%dos_delta = tundos%broadeningDelta
    if (tundos%defined) then
      params%delta = tundos%delta      ! delta for G.F.
    end if
    if (greendens%defined) then
      params%delta = greendens%delta   ! delta for G.F.
    end if
    
    ! ------------------------------------------------------------------------------
    !                    SETTING TRANSMISSION PARAMETERS
    ! ------------------------------------------------------------------------------

    if (tundos%defined) then

      l = size(tundos%ni)    
      params%ni(1:l) = tundos%ni(1:l)
      params%nf(1:l) = tundos%nf(1:l)

      ! setting of intervals and indeces for projected DOS
      nldos = size(tundos%dosOrbitals)
      call init_ldos(negf, nldos)
      do i = 1, nldos
         call set_ldos_indexes(negf, i, tundos%dosOrbitals(i)%data)   
      end do 
      
      nldos = size(tundos%dosOrbitals)      
      call init_ldos(negf, nldos)
      do j = 1, nldos
        call set_ldos_indexes(negf, j, tundos%dosOrbitals(j)%data)
      end do

      params%Emin =  tundos%Emin
      params%Emax =  tundos%Emax
      params%Estep = tundos%Estep

    endif
    
    ! Energy conversion only affects output units. 
    ! The library writes energies as (E * negf%eneconv) 
    params%eneconv = HAR 

    if (allocated(sizes)) call log_deallocate(sizes)

    call set_params(negf,params)

    !--------------------------------------------------------------------------
    !DAR begin - negf_init - TransPar to negf
    !--------------------------------------------------------------------------
    if (transpar%defined) then
      call tp2negf(transpar)
    else
      negf%tDephasingVE=transpar%tDephasingVE
      negf%tDephasingBP=transpar%tDephasingBP  
    end if
  
    !if(negf%tSpinDegeneracy) then !! If decomment - breaks autotest.
    !   negf%g_spin=2.
    !else
    !   negf%g_spin=1.
    !end if
    !--------------------------------------------------------------------------
    !DAR end
    !--------------------------------------------------------------------------

  end subroutine negf_init

  !------------------------------------------------------------------------------------------------!
  !> Initialize QT environment and variables without geometry (public)
  !------------------------------------------------------------------------------------------------!
  subroutine negf_init_nogeom(transpar, greendens, tundos, mpicomm, initinfo)
         
    Type(TTranspar), intent(IN) :: transpar
    Type(TGDFTBstructure) :: structure
    Type(TGDFTBGreenDensInfo) :: greendens
    Type(TGDFTBTunDos) :: tundos
    Type(mpifx_comm), intent(in) :: mpicomm
    logical, intent(OUT) :: initinfo
       
    ! local variables
    real(dp), dimension(:), allocatable :: pot, eFermi
    integer :: i,error, l, ncont, nc_vec(1), j, nldos
    integer, dimension(:), allocatable :: sizes  
    ! string needed to hold processor name
    character(:), allocatable :: hostname
    type(lnParams) :: params
    
    error = 0 
    initinfo = .true.

    call negf_mpi_init(mpicomm)
         
    if (transpar%defined) then
      ncont = transpar%ncont 
    else
      ncont = 0
    endif
      
    ! ------------------------------------------------------------------------------
    ! check that GS and contacts are created  (one day this will be automatic)
    ! ------------------------------------------------------------------------------
    if(transpar%defined .and. ncont.gt.0) then
                                
      open(199,file='GS/test',IOSTAT=error)
      if (error.ne.0) then 
        call mpifx_get_processor_name(hostname)
        write(*,*) "ERROR: please create a directory called GS on "//trim(hostname)
        initinfo = .false.; return  
      end if   
      close(199)
      open(199,file='contacts/test',IOSTAT=error)
      if (error.ne.0) then 
        call mpifx_get_processor_name(hostname)
        write(*,*) "ERROR: please create a directory called contacts "//trim(hostname)
        initinfo = .false.; return  
      end if         
      close(199)
    
    end if
    ! ------------------------------------------------------------------------------
    !! Set defaults and fill up the parameter structure with them
    
    call init_negf(negf)
    
    call get_params(negf, params)

    ! ------------------------------------------------------------------------------
    !This must be different for different initialisations, to be separated
    !Higher between transport and greendens is taken, temporary
    if (tundos%defined .and. greendens%defined) then                                          
       if (tundos%verbose.gt.greendens%verbose) then 
          params%verbose = tundos%verbose
       else
          params%verbose = greendens%verbose
       endif
    else
       if (tundos%defined) params%verbose = tundos%verbose
       if (greendens%defined) params%verbose = greendens%verbose
    end if
    ! ------------------------------------------------------------------------------
    ! This parameter is used to set the averall drop threshold in libnegf
    ! It affects especially transmission that is not accurate more than 
    ! this value.
    call set_drop(1.d-20)

 
    ! ------------------------------------------------------------------------------
    ! Assign spin degenracy and check consistency between different input blocks
    if (tundos%defined .and. greendens%defined) then
      if (tundos%gSpin .ne. greendens%gSpin) then
        write(*,*) "ERROR: spin degeneracy is not consistent between different input blocks"
        initinfo = .false.; return 
      else
        params%g_spin = real(tundos%gSpin) ! Spin degeneracy
      end if
    else if (tundos%defined) then
        params%g_spin = real(tundos%gSpin) ! Spin degeneracy
    else if (greendens%defined) then
        params%g_spin = real(greendens%gSpin) ! Spin degeneracy
    end if
 

    ! ------------------------------------------------------------------------------
    !            SETTING ELECTOCHEMICAL POTENTIALS INCLUDING BUILT-IN 
    ! ------------------------------------------------------------------------------
    ! Fermi level is given by the contacts. If no contacts => no transport,
    ! Then Fermi is defined by the Green solver
    if (transpar%defined) then
      pot = transpar%contacts(1:ncont)%potential
      eFermi = transpar%contacts(1:ncont)%eFermi(1) 
      do i = 1,ncont
        ! Built-in potential to equilibrate Fermi levels
        pot(i) = pot(i) + eFermi(i) - minval(eFermi(1:ncont))

        ! Setting contact temperatures
       
        if (transpar%contacts(i)%kbT .ge. 0.0_dp) then
          params%kbT(i) = transpar%contacts(i)%kbT
        else
          params%kbT(i) = 0.0_dp !structure%tempElec
          ! Make sure low temperatures (< 10K) converted to 0.
          ! This avoid numerical problems with contour integration
          if (params%kbT(i) < 3.0e-5_dp) then
            params%kbT(i) = 0.0_dp
          end if
        end if

        ! set parameters for wide band approximations
        params%FictCont(i) = transpar%contacts(i)%wideBand
        params%contact_DOS(i) = transpar%contacts(i)%wideBandDOS
      
        if (id0.and.transpar%verbose.gt.50) then
          write(*,"(/,'(negf_init) CONTACT INFO #',I0)")i
            if (params%FictCont(i)) then
              write(*,*) 'FICTICIOUS CONTACT '
              write(*,*) 'DOS: ', params%contact_DOS(i)
            end if
          write(*,*) '  Temperature: ', params%kbT(i)
          write(*,*) '  Potential (with built-in): ', pot(i)
          write(*,*) '  eFermi: ', eFermi(i) 
        endif 

      enddo

      ! Define electrochemical potentials
      params%mu(1:ncont) = eFermi(1:ncont) - pot(1:ncont)
      if (id0.and.transpar%verbose.gt.50) write(*,*)
      if (id0.and.transpar%verbose.gt.50) write(*,*) '  Electro-chemical potentials: ', params%mu(1:ncont)
      deallocate(pot)

    else !transpar not defined 
      params%mu(1) = greendens%oneFermi(1) 
    end if


    ! ------------------------------------------------------------------------------
    !                  SETTING COUNTOUR INTEGRATION PARAMETERS
    ! ------------------------------------------------------------------------------
    if (greendens%defined) then
      params%Ec = greendens%enLow           ! lowest energy
      params%Np_n(1:2) = greendens%nP(1:2)  ! contour npoints
      params%Np_real = greendens%nP(3)      ! real axis points
      params%n_kt = greendens%nkt           ! n*kT for Fermi
 
      !Read G.F. from very first iter
      if (greendens%readSGF .and. .not.greendens%saveSGF) params%readOldSGF=0  
      !compute G.F. at every iteration 
      if (.not.greendens%readSGF .and. .not.greendens%saveSGF) params%readOldSGF=1  
      !Default Write on first iter
      if (.not.greendens%readSGF .and. greendens%saveSGF) params%readOldSGF=2  

      if(any(params%kbT.gt.0) .and. greendens%nPoles.eq.0) then
         STOP 'ERROR: Number of Poles = 0 but T > 0' 
      else
         params%n_poles = greendens%nPoles
      end if
      if(all(params%kbT.eq.0)) then 
        params%n_poles = 0
      end if

    end if

    ! ------------------------------------------------------------------------------
    !! Setting the delta: priority on Green Solver, if present
    !! dos_delta is used by libnegf to smoothen T(E) and DOS(E)
    !! and is currently set in transmission
    params%dos_delta = tundos%broadeningDelta
    if (tundos%defined) then
      params%delta = tundos%delta      ! delta for G.F.
    end if
    if (greendens%defined) then
      params%delta = greendens%delta   ! delta for G.F.
    end if

    ! ------------------------------------------------------------------------------
    !                    SETTING TRANSMISSION PARAMETERS
    ! ------------------------------------------------------------------------------

    if (tundos%defined) then
 
      l = size(tundos%ni)    
      params%ni(1:l) = tundos%ni(1:l)
      params%nf(1:l) = tundos%nf(1:l)
 
      ! setting of intervals and indeces for projected DOS
      nldos = size(tundos%dosOrbitals)

      call init_ldos(negf, nldos)
 
      do i = 1, nldos
         call set_ldos_indexes(negf, i, tundos%dosOrbitals(i)%data)   
      end do 
       
      nldos = size(tundos%dosOrbitals)      
      call init_ldos(negf, nldos)
      do j = 1, nldos
  
        call set_ldos_indexes(negf, j, tundos%dosOrbitals(j)%data)
 
     end do

      params%Emin =  tundos%Emin
      params%Emax =  tundos%Emax
      params%Estep = tundos%Estep

    endif
 
    ! Energy conversion only affects output units. 
    ! The library writes energies as (E * negf%eneconv) 
    params%eneconv = HAR 

    if (allocated(sizes)) call log_deallocate(sizes)

    call set_params(negf,params)
 
    call tp2negf(transpar)

    if(negf%tSpinDegeneracy) then
       params%g_spin=2.
       negf%g_spin=2.
    else
       params%g_spin=1.
       negf%g_spin=1.
    end if

    negf%spin = 1
    
  end subroutine negf_init_nogeom

  !------------------------------------------------------------------------------------------------!  
  
  subroutine tp2negf(transpar)

    Type(TTranspar), intent(IN) :: transpar

    negf%tNoGeometry=transpar%tNoGeometry
    negf%tReadOverlap = transpar%tReadOverlap
    negf%tWriteDFTB=transpar%tWriteDFTB
    negf%tReadDFTB=transpar%tReadDFTB
    negf%tOrthonormal=transpar%tOrthonormal
    negf%tOrthonormalDevice=transpar%tOrthonormalDevice
    negf%tWriteOrthonormal=transpar%tWriteOrthonormal
    negf%tModel=transpar%tModel
    negf%tSpinDegeneracy=transpar%tSpinDegeneracy
    negf%tReadU=transpar%tReadU
    negf%tRead_negf_in=transpar%tRead_negf_in
    negf%NumStates=transpar%NumStates
    negf%FileName=transpar%HamiltonianFile
    negf%OverlapFile=transpar%OverlapFile
    negf%units_energy=transpar%units_energy
    negf%tManyBody=transpar%tManyBody
    negf%tElastic=transpar%tElastic
    negf%tDephasingVE=transpar%tDephasingVE
    negf%tDephasingBP=transpar%tDephasingBP
    negf%tZeroCurrent=transpar%tZeroCurrent
    negf%tMBNGF=transpar%tMBNGF
    negf%Mixing=transpar%Mixing
    negf%Tolerance=transpar%Tolerance
    negf%MaxIter=transpar%MaxIter
    negf%tWrite_negf_params=transpar%tWrite_negf_params
    negf%tWriteTagged=transpar%tWriteTagged
    allocate(negf%tranas%cont(transpar%ncont))
    negf%tranas%cont(:)%name=transpar%contacts(:)%name
    negf%tranas%cont(:)%tWriteSelfEnergy=transpar%contacts(:)%tWriteSelfEnergy
    negf%tranas%cont(:)%tReadSelfEnergy=transpar%contacts(:)%tReadSelfEnergy
    negf%tranas%cont(:)%tWriteSurfaceGF=transpar%contacts(:)%tWriteSurfaceGF
    negf%tranas%cont(:)%tReadSurfaceGF=transpar%contacts(:)%tReadSurfaceGF
    negf%tranas%cont(:)%tUnformatted=transpar%contacts(:)%tUnformatted
    negf%tranas%cont(:)%tWriteSeparatedSGF=transpar%contacts(:)%tWriteSeparatedSGF
    negf%tranas%cont(:)%tReadSeparatedSGF=transpar%contacts(:)%tReadSeparatedSGF 

  end subroutine tp2negf

  !------------------------------------------------------------------------------------------------!  
    
  subroutine negf_init_elph(elph)
  
    type(TElPh), intent(in) :: elph
    real(dp), dimension(:), allocatable :: coupling
    
    if (id0) write(*,*)

    allocate(coupling(size(elph%coupling)))
    coupling=elph%coupling
    
    if (elph%model .eq. 1) then
       if (id0) write(*,*) 'Setting local fully diagonal (FD) elastic dephasing model'
       call set_elph_dephasing(negf, coupling, elph%scba_niter)
    else if (elph%model .eq. 2) then
       if (id0) write(*,*) 'Setting local block diagonal (BD) elastic dephasing model'
       call set_elph_block_dephasing(negf, coupling, elph%orbsperatm, &
            elph%scba_niter)
    else if (elph%model .eq. 3) then
       if (id0) write(*,*) 'Setting overlap mask (OM) block diagonal elastic dephasing model'
       call set_elph_s_dephasing(negf, coupling, elph%orbsperatm, &
            elph%scba_niter)
    else
       write(*,*) "ERROR: el-ph model is not supported"
    endif

    negf%inter%Mixing=elph%Mixing
    negf%inter%scba_tol=elph%scba_tol
    
  end subroutine negf_init_elph
  !!$ !----------------------------------------------------------------------------
  !!$ subroutine negf_destroy_elph()
  !!$
  !!$    call destroy_elph_model(negf)
  !!$
  !!$ end subroutine negf_destroy_elph
  !!$ !----------------------------------------------------------------------------
  !!$

  !------------------------------------------------------------------------------------------------!  
  
  subroutine negf_init_bp(elph)
  
    type(TElPh), intent(in) :: elph
!    real(dp), dimension(:), allocatable :: coupling
    
    if (id0.and.negf%verbose.gt.30) write(*,*)
       
    if (elph%model .eq. 1) then
       if (id0.and.negf%verbose.gt.30) write(*,*) 'Setting local fully diagonal (FD) BP dephasing model'
       if(.not.allocated(negf%deph%bp%coupling)) &
            call log_allocate(negf%deph%bp%coupling,size(elph%coupling))      
       negf%deph%bp%coupling=elph%coupling
    else if (elph%model .eq. 2) then
       if (id0.and.negf%verbose.gt.30) write(*,*) 'Setting local block diagonal (BD) BP dephasing model'
       if (id0.and.negf%verbose.gt.30) write(*,*) 'NOT IMPLEMENTED! INTERRUPTED!'
       !call set_bp_block_dephasing(negf, elph%coupling, elph%orbsperatm, &
       !     elph%scba_niter)
       stop
    else if (elph%model .eq. 3) then
       if (id0.and.negf%verbose.gt.30) write(*,*) 'Setting overlap mask (OM) block diagonal BP dephasing model'
       if (id0.and.negf%verbose.gt.30) write(*,*) 'NOT IMPLEMENTED! INTERRUPTED!'
       !call set_bp_s_dephasing(negf, elph%coupling, elph%orbsperatm, &
       !     elph%scba_niter)
       stop
    else
       write(*,*) "ERROR: BP model is not supported"
    endif

    if (id0.and.negf%verbose.gt.30) write(*,*)'BP dephasing initialization is finished'

  end subroutine negf_init_bp
  
  !------------------------------------------------------------------------------------------------!
  
  subroutine negf_init_csr(iAtomStart, iNeighbor, nNeighbor, img2CentCell, orb)
    integer, intent(in) :: iAtomStart(:)    
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    type(TOrbitals), intent(in) :: orb

    if (allocated(csrHam%nzval)) call destroy(csrHam)
    call init(csrHam, iAtomStart, iNeighbor, nNeighbor, &
        &img2CentCell, orb) 
    if (allocated(csrOver%nzval)) call destroy(csrOver)    

    call init(csrOver, csrHam)

  end subroutine negf_init_csr
  
  !------------------------------------------------------------------------------------------------!
  
  subroutine negf_destroy()

    !write(*,'(A)') 'Release Negf Memory:'                                 !DAR
    call destruct(csrHam)
    call destruct(csrOver)
    call destroy_negf(negf)
    !call writePeakInfo(6)                                                 !DAR
    !call writeMemInfo(6)                                                  !DAR

  end subroutine negf_destroy

  !------------------------------------------------------------------------------------------------!
  
  subroutine negf_init_str(denseDescr, transpar, greendens, iNeigh, nNeigh, img2CentCell)

    use tranas_types_mbngf, only : Tmbngf
    
    Type(TDenseDescr), intent(in) :: denseDescr
    Type(TTranspar), intent(in) :: transpar
    Type(TGDFTBGreenDensInfo) :: greendens
    type(Tmbngf) :: mbngf_tmp
    Integer, intent(in) :: nNeigh(:)
    Integer, intent(in) :: img2CentCell(:)
    Integer, intent(in) :: iNeigh(0:,:)

    Integer, allocatable :: PL_end(:), cont_end(:), surf_end(:), cblk(:), ind(:)
    Integer, allocatable :: atomst(:), plcont(:)
    integer, allocatable :: minv(:,:)
    Integer :: natoms, ncont, nbl, iatm1, iatm2, iatc1, iatc2
    integer :: i, m, i1, j1, info
    integer, allocatable :: inRegion(:)

    iatm1 = transpar%idxdevice(1)
    iatm2 = transpar%idxdevice(2)
    ncont = transpar%ncont
    nbl = 0

    if (transpar%defined) then
       nbl = transpar%nPLs
    else if (greendens%defined) then
       nbl = greendens%nPLs
    endif

    if (nbl.eq.0) then
      call error('Internal ERROR: nbl = 0 ?!')
    end if

    natoms = size(denseDescr%iatomstart) - 1

!NEW    call check_pls(transpar, greendens, natoms, iNeigh, nNeigh, img2CentCell, info)

    allocate(PL_end(nbl))
    allocate(atomst(nbl+1))
    allocate(plcont(nbl))
    allocate(cblk(ncont))
    allocate(cont_end(ncont))
    allocate(surf_end(ncont))
    allocate(ind(natoms+1))
    allocate(minv(nbl,ncont))

    ind(:) = DenseDescr%iatomstart(:) - 1
    minv = 0
    cblk = 0

    do i = 1, ncont
       cont_end(i) = ind(transpar%contacts(i)%idxrange(2)+1)
       surf_end(i) = ind(transpar%contacts(i)%idxrange(1))
    enddo

    if (transpar%defined) then
      do i = 1, nbl-1
        PL_end(i) = ind(transpar%PL(i+1))
      enddo
      atomst(1:nbl) = transpar%PL(1:nbl)
      PL_end(nbl) = ind(transpar%idxdevice(2)+1)
      atomst(nbl+1) = iatm2 + 1
    else if (greendens%defined) then
      do i = 1, nbl-1
        PL_end(i) = ind(greendens%PL(i+1))
      enddo
      atomst(1:nbl) = greendens%PL(1:nbl)
      PL_end(nbl) = ind(natoms+1)
      atomst(nbl+1) = natoms + 1
    endif

    if (transpar%defined .and. ncont.gt.0) then

      if(.not.transpar%tNoGeometry) then

       ! For each PL finds the min atom index among the atoms in each contact
       ! At the end the array minv(iPL,iCont) can have only one value != 0
       ! for each contact and this is the interacting PL
       ! NOTE: the algorithm works with the asymmetric neighbor-map of dftb+
       !       because atoms in contacts have larger indices than in the device
       do m = 1, transpar%nPLs
          ! Loop over all PL atoms
          do i = atomst(m), atomst(m+1)-1

             ! Loop over all contacts
             do j1 = 1, ncont

                iatc1 = transpar%contacts(j1)%idxrange(1)
                iatc2 = transpar%contacts(j1)%idxrange(2)

                i1 = minval(img2CentCell(iNeigh(1:nNeigh(i),i)), &
                    & mask = (img2CentCell(iNeigh(1:nNeigh(i),i)).ge.iatc1 .and. &
                    & img2CentCell(iNeigh(1:nNeigh(i),i)).le.iatc2) )

                if (i1.ge.iatc1 .and. i1.le.iatc2) then
                    minv(m,j1) = j1
                endif

             end do
          end do
       end do


       do j1 = 1, ncont

         if (all(minv(:,j1) == 0)) then
           write(stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(stdOut,*) 'WARNING: contact',j1,' does no interact with any PL '
           write(stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           minv(1,j1) = j1
         end if

         if (count(minv(:,j1).eq.j1).gt.1) then
           write(stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(stdOut,*) 'ERROR: contact',j1,' interacts with more than one PL'
           write(stdOut,*) '       check structure and increase PL size         '
           write(stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!NEW           call error("")
         end if

         do m = 1, transpar%nPLs
           if (minv(m,j1).eq.j1) then
             cblk(j1) = m
           end if
         end do

       end do

      else

         cblk=transpar%cblk

      end if

      if (transpar%verbose.gt.50) then
        write(stdOut,"(/,'Structure info:')")
        write(stdOut,*) '  Number of PLs:',nbl
        write(stdOut,*) '  PLs coupled to contacts:',cblk(1:ncont)
      end if
      
    end if

    call init_structure(negf, ncont, cont_end, surf_end, nbl, PL_end, cblk)

    if(negf%tMBNGF) then
      mbngf_tmp%str=negf%str
      mbngf_tmp%tHartreeFock=transpar%tHartreeFock
      mbngf_tmp%tRPA=transpar%tRPA
      if(.not.allocated(negf%mbngf)) allocate(negf%mbngf, source=mbngf_tmp)
    end if

    deallocate(PL_end)
    deallocate(plcont)
    deallocate(atomst)
    deallocate(cblk)
    deallocate(cont_end)
    deallocate(surf_end)
    deallocate(ind)
    deallocate(minv)

  end subroutine negf_init_str

  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  
  !------------------------------------------------------------------------------------------------!
  !> INTERFACE subroutine to call Density Matrix computation (called only from this module)
  !------------------------------------------------------------------------------------------------!
  subroutine negf_density(miter, spin, nkpoint, HH, SS, mu, tundos, DensMat, EnMat)
  !DAR - tundos is added, it is necessary for 'call negf_init_elph(tundos%elph)'

    type(TTraNaS) :: tranas
    Type(TGDFTBTunDos), intent(IN) :: tundos                                !DAR
    integer, intent (in) :: miter          ! SCC step (used in SGF)
    integer, intent (in) :: spin           ! spin component (SGF)
    integer, intent (in) :: nkpoint        ! nk point (used in SGF)
    type(z_CSR), intent(in) :: HH          ! Hamiltonian
    type(z_CSR), intent(in) :: SS          ! Overlap
    real(dp), intent(in) :: mu(:)       
    type(z_CSR), optional :: DensMat   ! Density matrix (See NOTE)
    type(z_CSR), optional :: EnMat     ! Energy weighted DM (See NOTE)

    type(lnParams) :: params
    integer :: ncont

    call get_params(negf, params)
      
    params%iteration = miter
    params%kpoint = nkpoint
    params%spin = spin
    params%DorE='N'
    ncont=negf%str%num_conts
    params%mu(1:ncont) = mu(1:ncont)
    
    if(present(DensMat)) then
       params%DorE = 'D'
       call set_params(negf,params)
       call pass_DM(negf,rho=DensMat)
       if (id0.and.params%verbose.gt.30) then
         write(*,'(80("-"))')
         write(*,*) '                        COMPUTATION OF THE DENSITY MATRIX             '
         write(*,'(80("-"))') 
       endif
    endif
    if(present(EnMat)) then
       params%DorE = 'E'
       call set_params(negf,params)
       call pass_DM(negf,rhoE=EnMat)
       if (id0.and.params%verbose.gt.30) then
         write(*,'(80("-"))')
         write(*,*) '                     COMPUTATION OF E-WEIGHTED DENSITY MATRIX         '
         write(*,'(80("-"))') 
       endif
    endif
    if (present(DensMat).and.present(EnMat)) then
       params%DorE  = 'B'
       call set_params(negf,params)
       stop 'UNSUPPORTED CASE in negf_density'
    endif

    if (params%DorE.eq.'N') return 
    ! ---------------------------------------------

    call pass_HS(negf,HH,SS)
    if(negf%tDephasingVE) call negf_init_elph(tundos%elph)                  !DAR
    if(negf%tDephasingBP) call negf_init_bp(tundos%bp)                      !DAR
    tranas%negf = negf
    call compute_density_dft(tranas)
    negf = tranas%negf    
    
    call destroy_matrices(negf)

  end subroutine negf_density

  !------------------------------------------------------------------------------------------------!
  !> INTERFACE subroutine to call LDOS computation (called only from this module)
  !------------------------------------------------------------------------------------------------!
  subroutine negf_ldos(HH, SS, spin, kpoint, wght, ledos)

    type(TTraNaS) :: tranas
    type(z_CSR), intent(in) :: HH, SS
    integer, intent(in) :: spin      ! spin index 
    integer, intent(in) :: kpoint        ! kp index 
    real(dp), intent(in) :: wght      ! kp weight 
    real(dp), dimension(:,:), pointer :: ledos
    type(lnParams) :: params 

    call get_params(negf, params)

    if (id0) then
      write(*,*)
      write(*,'(73("="))')
      write(*,*) '                   COMPUTING  LOCAL  DOS          '
      write(*,'(73("="))') 
      write(*,*)
    end if

    params%spin = spin
    params%kpoint = kpoint
    params%wght = wght

    if (associated(negf%ldos_mat)) then
      call log_deallocatep(negf%ldos_mat)
      negf%ldos_mat=> null()
    endif
    ledos=>null()
    
    call pass_HS(negf,HH,SS)

    tranas%negf = negf
    call calcLDOS(tranas)
    negf = tranas%negf
    
    call destroy_matrices(negf)

    call associate_ldos(negf, ledos)

  end subroutine negf_ldos
  

  !----------------------------------------------------------------------------
  ! DEBUG routine dumping H and S on file in Matlab format
  !----------------------------------------------------------------------------
  subroutine negf_dumpHS(HH,SS)
    type(z_CSR), intent(in) :: HH, SS

    write(*,*) 'Dumping H and S on files...'    
    open(1121, file='HH.dat')
    write(1121,*) '% Size =',HH%nrow, HH%ncol
    write(1121,*) '% Nonzeros =',HH%nnz
    write(1121,*) '% '
    write(1121,*) 'zzz = ['
    call printcsr(1121,HH)
    write(1121,*) ']'
    close(1121) 
    open(1121, file='SS.dat')
    write(1121,*) '% Size =',SS%nrow, SS%ncol
    write(1121,*) '% Nonzeros =',SS%nnz
    write(1121,*) '% '
    write(1121,*) 'zzz = ['
    call printcsr(1121,SS)
    write(1121,*) ']'
    close(1121) 
  end subroutine negf_dumpHS
  
  !------------------------------------------------------------------------------------------------!
  !> INTERFACE subroutine to call current computation (called only from this module)
  !------------------------------------------------------------------------------------------------!   
  subroutine negf_current(HH, SS, spin, kpoint, wght, tunn, ledos, currents, tundos)
  !DAR - tundos is added, it is necessary for 'call negf_init_elph(tundos%elph)'
                                                               
    Type(TGDFTBTunDos), intent(IN) :: tundos                    !DAR
    type(z_CSR), intent(inout) :: HH, SS                        !DAR in -> inout
    integer, intent(in) :: spin      ! spin index 
    integer, intent(in) :: kpoint        ! kp index 
    real(dp), intent(in) :: wght      ! kp weight 
    real(dp), dimension(:,:), pointer :: tunn
    real(dp), dimension(:,:), pointer :: ledos
    real(dp), dimension(:), pointer :: currents 

    integer :: i
    type(lnParams) :: params

    type(TTraNaS) :: tranas

    !DAR begin - Preparation of the Hamiltonian
    !--------------------------------------------------------------------------!
    
    integer :: j,k,l,m,n,NumStates,icont
    real(dp), allocatable :: H(:,:),S(:,:)

    if(negf%tReadDFTB) call ReadDFTB
      
    if(negf%tModel) call ReadModel

    if(negf%tReadU) call ReadU

    if(negf%tOrthonormal) call Orthogonalization

    if(negf%tOrthonormalDevice) call Orthogonalization_dev

    if(negf%tElastic.and.(negf%tReadDFTB.or.negf%tModel.or.negf%tOrthonormal.or.negf%tOrthonormalDevice)) &
         call MakeHHSS(HH,SS)

    if(negf%tManyBody) call MakeHS_dev

    if(negf%tRead_negf_in) call ReadLibNEGF 
         
    !--------------------------------------------------------------------------!
    !DAR end

    call get_params(negf, params)

    params%spin = spin
    params%kpoint = kpoint
    params%wght = wght
    
    call set_params(negf, params)

    call pass_HS(negf,HH,SS)

    !DAR begin
    if((.not.negf%tElastic).and.(.not.negf%tManyBody)) then
         write(*,*)'Current is not calculated!'
         write(*,*)'Choose "Elastic = Yes" or "ManyBody = Yes"!'
         write(*,*)'Program is terminated!'
         stop
    end if 

    if(negf%tWrite_negf_params) call check_negf_params

    if(negf%tDephasingVE) call negf_init_elph(tundos%elph)                  
    if(negf%tDephasingBP) call negf_init_bp(tundos%bp)                      

    !--------------------------------------------------------------------------!
    ! Calculations with TraNaS
    !--------------------------------------------------------------------------!

    tranas%negf = negf 
    
    if (negf%tMBNGF) call calcMBNGF(tranas)

    if ((.not.allocated(negf%inter)).and.(.not.negf%tDephasingBP).and.(.not.negf%tMBNGF)) then
      if (tundos%writeLDOS.and.(.not.tundos%writeTunn)) call calcLDOS(tranas)
    else
      if (tundos%writeLDOS) call calcLDOS(tranas) 
    end if
       
    if (tundos%writeTunn) then
      if ((.not.allocated(negf%inter)).and.(.not.negf%tDephasingBP).and.(.not.negf%tMBNGF)) then
        call calcLandauer(tranas)
      else
        call calcMeirWingreen(tranas)
      end if
    end if

    negf = tranas%negf
    
    if (tundos%writeTunn) call associate_current(negf, currents)
    call associate_ldos(negf, ledos)  
    if (tundos%writeTunn) call associate_transmission(negf, tunn)
    
  end subroutine negf_current

  !------------------------------------------------------------------------------------------------!

  !> Subroutine to call transport computation without geometry (public)
  subroutine tranasNoGeom(mpicomm, tundos, transpar) 

    Type(TGDFTBTunDos), intent(IN) :: tundos
    Type(TTransPar), intent(in) :: transpar
    type(z_CSR) :: HH, SS
    type(mpifx_comm), intent(in) :: mpicomm
    real(dp), dimension(:,:), pointer :: tunn
    real(dp), dimension(:,:), pointer :: ledos
    real(dp), dimension(:), pointer :: currents
    type(unit) :: unitOfEnergy        ! Set the units of H
    type(unit) :: unitOfCurrent       ! Set desired units for Jel

    integer :: i,j,k,l,m,n,NumStates,icont
    real(8), dimension(:), allocatable :: coupling   
    real(dp), allocatable :: H(:,:),S(:,:)

    type(TTraNaS) :: tranas
    
    ! Tagged writer
    type(TTaggedWriter) :: taggedWriter

    unitOfEnergy%name = "H"
    unitOfCurrent%name = "A"
    
    if (id0.and.negf%verbose.gt.0) then
      write(*,*)
      write(*,'(80("="))')
      write(*,*) '                            COMPUTATION OF TRANSPORT         '
      write(*,'(80("="))') 
      write(*,*)
      write(*,"('Transport is started.')")
      write(*,"('Number of states (with electrodes)      = ',I5)")negf%NumStates
      write(*,"('Number of states in the central region  = ',I5)")negf%str%central_dim
    endif   

    !--------------------------------------------------------------------------!
    ! Preparation of the Hamiltonian
    !--------------------------------------------------------------------------!
    
    if(negf%tReadDFTB) call ReadDFTB
      
    if(.not.allocated(negf%H_all)) allocate(negf%H_all(negf%NumStates,negf%NumStates))
    negf%H_all=transpar%H_all
    if (id0.and.negf%verbose.gt.100) then
    write(*,"('Hamiltonian:')")   
    do i=1,negf%NumStates
       write(*,"(10000ES16.8)")negf%H_all(i,1:negf%NumStates)
    end do
    end if    
    if (negf%tReadOverlap) then
      if (id0.and.negf%verbose.gt.50) write(*,"(' Overlap is red from the file ',A)")trim(negf%OverlapFile)
      open(11,file=negf%OverlapFile,action="read")
      if(.not.allocated(negf%S_all)) allocate(negf%S_all(negf%NumStates,negf%NumStates))
      negf%S_all=0.0_dp
      do i=1,negf%NumStates
        read(11,*)negf%S_all(i,1:negf%NumStates)
      end do
      close(11)
    else   
      if (.not.allocated(negf%S_all)) allocate(negf%S_all(negf%NumStates,negf%NumStates))
      negf%S_all=transpar%S_all
    end if
    if (id0.and.negf%verbose.gt.100) then
    write(*,"('Overlap:')")   
    do i=1,negf%NumStates
       write(*,"(10000ES16.8)")negf%S_all(i,1:negf%NumStates)
    end do
    end if

    if(negf%tReadU) call ReadU

    if(negf%tOrthonormal) call Orthogonalization

    if(negf%tOrthonormalDevice) call Orthogonalization_dev

    if(negf%tElastic.and.(negf%tReadDFTB.or.negf%tModel.or.negf%tOrthonormal.or.negf%tOrthonormalDevice)) &
         call MakeHHSS(HH,SS)

    if(negf%tManyBody) call MakeHS_dev

    if(negf%tRead_negf_in) call ReadLibNEGF

    !--------------------------------------------------------------------------!
    
    call pass_HS(negf,HH,SS)
   
    if(negf%tWrite_negf_params) call check_negf_params

    if(negf%tDephasingVE) call negf_init_elph(tundos%elph)                  
    if(negf%tDephasingBP) call negf_init_bp(tundos%bp)

    negf%mpicomm = mpicomm !DAR! for compute_mbngf temporary solution

    !--------------------------------------------------------------------------!
    ! Calculations with TraNaS
    !--------------------------------------------------------------------------!

    tranas%negf = negf
    tranas%input = transpar%tranas_input
 
    if (negf%tMBNGF) call calcMBNGF(tranas)

    if ((.not.allocated(negf%inter)).and.(.not.negf%tDephasingBP).and.(.not.negf%tMBNGF)) then
      if (tundos%writeLDOS.and.(.not.tundos%writeTunn)) call calcLDOS(tranas)
    else
      if (tundos%writeLDOS) call calcLDOS(tranas) 
    end if

    if (tundos%writeTunn) then
      if ((.not.allocated(negf%inter)).and.(.not.negf%tDephasingBP).and.(.not.negf%tMBNGF)) then
        call calcLandauer(tranas)
      else
        call calcMeirWingreen(tranas)
      end if
    end if

    negf = tranas%negf

    if (tundos%writeTunn) call associate_current(negf, currents)
    call associate_ldos(negf, ledos)
    if (tundos%writeTunn) call associate_transmission(negf, tunn)

    call mpifx_barrier(mpicomm)
    
    !-------------------------------------------------------------------------   
    !WriteSelfEnergy / WriteSurfaceGF
    !-------------------------------------------------------------------------
    
    do icont=1,negf%str%num_conts
      if(negf%tranas%cont(icont)%tWriteSelfEnergy) &
         call mpifx_allreduceip(mpicomm, negf%tranas%cont(icont)%SelfEnergy, MPI_SUM)
    end do
  
    if (id0) then
       do icont=1,negf%str%num_conts
          if(negf%tranas%cont(icont)%tWriteSelfEnergy) then
            if (negf%tranas%cont(icont)%tUnformatted) then
              open(14,form="unformatted",file=trim(negf%tranas%cont(icont)%name)//'-SelfEnergy.mgf' &
                   ,action="write")
              do i = 1, size(negf%en_grid)
                write(14)real(negf%en_grid(i)%Ec),negf%tranas%cont(icont)%SelfEnergy(:,:,i)
              end do
            else
              open(14,file=trim(negf%tranas%cont(icont)%name)//'-SelfEnergy.mgf' &
                   ,action="write")
              do i = 1, size(negf%en_grid)
                write(14,*)real(negf%en_grid(i)%Ec),negf%tranas%cont(icont)%SelfEnergy(:,:,i)
              end do
            end if 
            close(14)
            if (negf%verbose.gt.30) write(*,"('The retarded contact self-energy is written into the file ',A)") &
                  trim(negf%tranas%cont(icont)%name)//'-SelfEnergy.mgf'
          end if
       end do
    end if

    do icont=1,negf%str%num_conts
       if(negf%tranas%cont(icont)%tWriteSurfaceGF) &
            call mpifx_allreduceip(mpicomm, negf%tranas%cont(icont)%SurfaceGF, MPI_SUM)
    end do

    !debug begin
    !   if (id0) then
    !   do i = 1, size(negf%en_grid)
    !     do icont=1,negf%str%num_conts  
    !        !print *, negf%tranas%cont(icont)%tWriteSurfaceGF
    !        if(negf%tranas%cont(icont)%tWriteSurfaceGF) then
    !          print *, 'SurfaceGF',icont,'IndexEnergy=',i
    !          do j=1,size(negf%tranas%cont(icont)%SurfaceGF,1)
    !            print *, negf%tranas%cont(icont)%SurfaceGF(j,1:size(negf%tranas%cont(icont)%SurfaceGF,1),i)
    !          end do
    !        end if 
    !     end do
    !   end do
    !   end if
    !debug end
    
    if (id0) then
       do icont=1,negf%str%num_conts
          if(negf%tranas%cont(icont)%tWriteSurfaceGF) then
            if (negf%tranas%cont(icont)%tUnformatted) then
              open(14,form="unformatted",file=trim(negf%tranas%cont(icont)%name)//'-SurfaceGF.mgf' &
                   ,action="write")
              do i = 1, size(negf%en_grid)
                write(14)real(negf%en_grid(i)%Ec),negf%tranas%cont(icont)%SurfaceGF(:,:,i)
              end do
            else
              open(14,file=trim(negf%tranas%cont(icont)%name)//'-SurfaceGF.mgf' &
                   ,action="write")
              do i = 1, size(negf%en_grid)
                write(14,*)real(negf%en_grid(i)%Ec),negf%tranas%cont(icont)%SurfaceGF(:,:,i)
              end do
            end if 
            close(14)
            if (negf%verbose.gt.30) write(*,"('The retarded surface Green function is written into the file ',A)") &
                  trim(negf%tranas%cont(icont)%name)//'-SurfaceGF.mgf'
          end if
       end do
    end if

    !-------------------------------------------------------------------------

    if (tundos%writeTunn) then   
    
      call mpifx_allreduceip(mpicomm, currents, MPI_SUM)

      currents = currents * convertCurrent(unitOfEnergy, unitOfCurrent) 
    
      if (id0.and.negf%verbose.gt.0) then
        write(*,"(/,'Current is calculated')") 
        do i=1, size(currents)
          write(*,'(1x,a,i3,i3,a,ES14.5,a,a)') &
             & ' contacts: ',negf%ni(i),negf%nf(i), &
             & '   current: ', currents(i),' ',unitOfCurrent%name
        enddo
      endif

      call mpifx_allreduceip(mpicomm, tunn, MPI_SUM)

      if (id0 .and. tundos%writeTunn) then  
         open(65000,file='transmission.dat')
         do i=1,size(tunn,1)
            write(65000,'(E18.8E3)',ADVANCE='NO') (negf%Emin+(i-1)*negf%Estep)*27.21138469
            do j=1,size(tunn,2)
               write(65000,'(E18.8E3)',ADVANCE='NO') tunn(i,j)
            enddo
            write(65000,*)
         enddo
         close(65000)
      endif

      if(negf%tZeroCurrent) then                                                
         call mpifx_allreduceip(mpicomm, negf%tunn_mat_bp, MPI_SUM)        
        if (id0 .and. tundos%writeTunn) then  
           open(65000,file='transmission_bp.dat')
           do i=1,size(negf%tunn_mat_bp,1)
              write(65000,'(E18.8E3)',ADVANCE='NO') (negf%Emin+(i-1)*negf%Estep)*27.21138469
              do j=1,size(negf%tunn_mat_bp,2)
                 write(65000,'(E18.8E3)',ADVANCE='NO') negf%tunn_mat_bp(i,j)
              enddo
              write(65000,*)
           enddo
           close(65000)
        endif 
      endif

    end if

    if (tundos%writeLDOS) call mpifx_allreduceip(mpicomm, ledos, MPI_SUM)

    if (id0 .and. tundos%writeLDOS) then
       open(65000,file='localDOS.dat')
       do i=1,size(ledos,1)
          write(65000,'(E18.8E3)',ADVANCE='NO') (negf%Emin+(i-1)*negf%Estep)*27.21138469
          do j=1,size(ledos,2)
             write(65000,'(E18.8E3)',ADVANCE='NO') ledos(i,j)
          enddo
          write(65000,*)
       enddo
       close(65000)
    end if

    if (negf%tWriteTagged) then
    
      call TTaggedWriter_init(taggedWriter)
    
      open(unit=10, file='autotest.tag', action="write", status="replace")

      if (id0 .and. tundos%writeTunn) call taggedWriter%write(10, tagLabels%tag_transmission, tunn)
    
!    if (tPeriodic) then
!      call writeTagged(10, tag_volume, cellVol)
!    end if
!    if (allocated(derivs)) then
!      call writeTagged(10, tag_forceTot, -derivs)
!    end if

      close(10)

    end if
        
  end subroutine tranasNoGeom

  !------------------------------------------------------------------------------------------------!
  !> Calculates density matrix with Green functions (public)
  !------------------------------------------------------------------------------------------------!
  subroutine calcdensity_green(iSCCIter, mpicomm, groupKS, ham, over, &
      & desc, iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, &
      & cellVec, orb, nEl, tempElec, kPoints, kWeights, &
      & rho, Eband, Ef, E0, TS, mu, tundos)
      !DAR - tundos is added,
      !      it is necessary for 'call negf_init_elph(tundos%elph)' in 'call negf_current'
    
    Type(TGDFTBTunDos), intent(IN) :: tundos                                !DAR
    integer, intent(in) :: iSCCIter
    type(mpifx_comm), intent(in) :: mpicomm
    integer, intent(in) :: groupKS(:,:)
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: desc(1) !Descriptor is not needed for now
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb 
    real(dp), intent(in) :: nEl(:), tempElec, kPoints(:,:), kWeights(:)
    real(dp), intent(out) :: rho(:,:)
    real(dp), intent(out) :: Eband(:), Ef(:), E0(:), TS(:)

    integer :: nSpin, nKPoint, nKS, iK, iS, iKS, ncont
    real(dp) :: prefac, EBandTmp, EfTmp, TSTmp, E0Tmp
    type(z_CSR) :: csrDens
    !! We need this now for different fermi levels in colinear spin
    !! Note: the spin polirized does not work with
    !! built-int potentials (the unpolarized does) in the poisson
    real(dp), intent(in) :: mu(:,:)
    type(lnParams) :: params

    call get_params(negf, params)

    nKS = size(groupKS, dim=2)
    nKPoint = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)
    ncont = size(mu,1)
    rho = 0.0_dp

    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)
      
      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)
      call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)

      call negf_density(iSCCIter, iS, iKS, csrHam, csrOver, mu(:,iS), tundos, DensMat=csrDens)
      !DAR - tundos is added,
      !      it is necessary for 'call negf_init_elph(tundos%elph)' in 'negf_current'

 
      ! NOTE:
      ! unfold adds up to rho the csrDens(k) contribution
      !
      call unfoldFromCSR(rho(:,iS), csrDens, kPoints(:,iK), &
          kWeights(iK), iAtomStart, iPair, iNeighbor, nNeighbor, &
          img2CentCell, iCellVec, cellVec, orb)

      call destruct(csrDens)

      !! Set some fake energies:
      Eband(iS) = 0.0_dp
      Ef(iS) = 0.0_dp
      TS(iS) = 0.0_dp
      E0(iS) = 0.0_dp

    end do

    do iS = 1, nSpin
      ! In place all-reduce of the density matrix
      call mpifx_allreduceip(mpicomm, rho(:,iS), MPI_SUM)
    end do
  
  end subroutine calcdensity_green

  !------------------------------------------------------------------------------------------------!
  !> Calculates E-density matrix with Green functions (public)
  !------------------------------------------------------------------------------------------------!
  subroutine calcEdensity_green(iSCCIter, mpicomm, groupKS, ham, over, &
      & iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, &
      & cellVec, orb, kPoints, kWeights, rhoE, mu, tundos)
      !DAR - tundos is added,
      !      it is necessary for 'call negf_init_elph(tundos%elph)' in 'call negf_current'
    
    Type(TGDFTBTunDos), intent(IN) :: tundos                                !DAR
    integer, intent(in) :: iSCCIter
    type(mpifx_comm), intent(in) :: mpicomm
    integer, intent(in) :: groupKS(:,:)
    real(dp), intent(in) :: ham(:,:), over(:)
    !integer, intent(in) :: desc(1) !Descriptor is not needed for now
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb   !Needs only orb%nOrbAtom, orb%mOrb    
    real(dp), intent(in) :: kPoints(:,:), kWeights(:)
    real(dp), intent(out) :: rhoE(:)

    integer :: nSpin, nKPoint, nKS, iK, iS, iKS, ncont
    real(dp) :: prefac, EBandTmp, EfTmp, TSTmp, E0Tmp
    type(z_CSR) :: csrEDens
    !! We need this now for different fermi levels in colinear spin
    !! Note: the spin polirized does not work with
    !! built-int potentials (the unpolarized does) in the poisson
    !! I do not set the fermi because it seems that in libnegf it is 
    !! not really needed
    real(dp), intent(in) :: mu(:,:)
    type(lnParams) :: params

    call get_params(negf, params)

    nKS = size(groupKS, dim=2)
    nKPoint = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)
    ncont = size(mu,1)
    rhoE = 0.0_dp

    
    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)
      
      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)
      call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)
      
      call negf_density(iSCCIter, iS, iKS, csrHam, csrOver, mu(:,iS), tundos, EnMat=csrEDens)
      !DAR - tundos is added,
      !      it is necessary for 'call negf_init_elph(tundos%elph)' in 'negf_current')

      ! NOTE:
      ! unfold adds up to rhoEPrim the csrEDens(k) contribution
      !
      call unfoldFromCSR(rhoE, csrEDens, kPoints(:,iK), &
          kWeights(iK), iAtomStart, iPair, iNeighbor, nNeighbor, &
          img2CentCell, iCellVec, cellVec, orb)

      call destruct(csrEDens)
      
    end do

    ! In place all-reduce of the density matrix
    call mpifx_allreduceip(mpicomm, rhoE, MPI_SUM)
    
  end subroutine calcEdensity_green

  !------------------------------------------------------------------------------------------------!
  !>
  !------------------------------------------------------------------------------------------------!
  subroutine calcPDOS_green(mpicomm, groupKS, ham, over, &
      & desc, iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, &
      & cellVec, orb, nEl, tempElec, kPoints, kWeights, ldosTot, writeLDOS)
    integer, intent(in) :: groupKS(:,:)
    type(mpifx_comm), intent(in) :: mpicomm
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: desc(1) !Descriptor is not needed for now
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: nEl(:), tempElec, kPoints(:,:), kWeights(:)
    real(dp), allocatable, intent(inout) :: ldosTot(:,:)
    logical, intent(in) :: writeLDOS

    real(dp), pointer    :: ldosMat(:,:)=>null()
    integer :: iKS, iK, iS, nKS, nKPoint, nSpin, ii, jj, err
    type(lnParams) :: params

    call get_params(negf, params)
    
    nKS = size(groupKS, dim=2)
    nKPoint = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)

    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)
      
      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)
     
      call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)


      call negf_ldos(csrHam, csrOver, iS, iK, kWeights(iK), ldosMat)
      
      if (associated(ldosMat)) then
        if(.not.allocated(ldosTot)) then
          allocate(ldosTot(size(ldosMat,1),size(ldosMat,2)), stat=err)
          if (err/=0) STOP 'Allocation error (ldosTot)'
          ldosTot = 0.0_dp
        endif
        ldosTot = ldosTot + ldosMat
      endif
      
    end do

    if (allocated(ldosTot)) then
            
      call mpifx_allreduceip(mpicomm, ldosTot, MPI_SUM)
      
      ! Write Total localDOS on a separate file (optional)
      if (id0 .and. writeLDOS) then
        open(65000,file='localDOS.dat')
        do ii=1,size(ldosTot,1)
          write(65000,*) (params%Emin+(ii-1)*params%Estep) * HAR, &
              (ldosTot(ii,jj), jj=1,size(ldosTot,2))
        enddo
        close(65000)
      endif
    else
      allocate(ldosTot(0,0))
    endif
    
    
  end subroutine calcPDOS_green
     
  !------------------------------------------------------------------------------------------------!
  !> Subroutine to call current computation (public)
  !------------------------------------------------------------------------------------------------!
  subroutine calc_current(mpicomm, groupKS, ham, over, &
      & desc, iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, &
      & cellVec, orb, nEl, tempElec, kPoints, kWeights, tunnTot,&
      & ldosTot, currTot, writeTunn, writeLDOS, mu, tundos)
      !DAR - tundos is added,
      !      it is necessary for 'call negf_init_elph(tundos%elph)' in 'call negf_current'
    
    Type(TGDFTBTunDos), intent(IN) :: tundos                                !DAR
    integer, intent(in) :: groupKS(:,:)
    type(mpifx_comm), intent(in) :: mpicomm
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: desc(1) !Descriptor is not needed for now
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: nEl(:), tempElec, kPoints(:,:), kWeights(:)
    real(dp), allocatable, intent(inout) :: tunnTot(:,:), ldosTot(:,:)
    real(dp), allocatable :: tunnSKRes(:,:,:), ldosSKRes(:,:,:)
    real(dp), allocatable, intent(inout) :: currTot(:)
    logical, intent(in) :: writeLDOS
    logical, intent(in) :: writeTunn
    !! We need this now for different fermi levels in colinear spin
    !! Note: the spin polirized does not work with
    !! built-int potentials (the unpolarized does) in the poisson
    !! I do not set the fermi because it seems that in libnegf it is 
    !! not really needed
    real(dp), intent(in) :: mu(:,:)
    
    real(dp), pointer    :: tunnMat(:,:)=>null()
    real(dp), pointer    :: ldosMat(:,:)=>null()
    real(dp), pointer    :: currVec(:)=>null()    
    integer :: iKS, iK, iS, nKS, nKPoint, nSpin, ii, jj, err, ncont
    type(unit) :: unitOfEnergy        ! Set the units of H
    type(unit) :: unitOfCurrent       ! Set desired units for Jel
    type(lnParams) :: params

    integer :: i,j,k,NumStates,icont                                       !DAR

    call get_params(negf, params)

    unitOfEnergy%name = "H"
    unitOfCurrent%name = "A"

       nKS = size(groupKS, dim=2)
       nKPoint = size(kPoints, dim=2)
       nSpin = size(ham, dim=2)
       ncont = size(mu,1)

  ! Start iKS     
  do iKS = 1, nKS 
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS) 

      params%mu(:ncont) = mu(:ncont,iS) 
   
      call set_params(negf, params)

      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, &
           &iPair, iNeighbor, nNeighbor, img2CentCell, &
           &iCellVec, cellVec, orb)
    
      call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, &
           &iPair, iNeighbor, nNeighbor, img2CentCell, &
           &iCellVec, cellVec, orb)

      !DAR begin - calc_current
      !------------------------------------------------------------------------

      if(negf%NumStates.eq.0) negf%NumStates=csrHam%ncol
      
      if (id0.and.params%verbose.gt.0) then
      write(*,*)
      write(*,'(80("="))')
      write(*,*) '                            COMPUTATION OF TRANSPORT         '
      write(*,'(80("="))') 
      write(*,*)
      write(*,"('Transport is started.')")
      write(*,"('Number of states (with electrodes)      = ',I5)")negf%NumStates
      write(*,"('Number of states in the central region  = ',I5)")negf%str%central_dim
      endif
      
      if( negf%tWriteDFTB &
      .or.(negf%tManyBody.and.(.not.negf%tReadDFTB).and.(.not.negf%tModel).and.(.not.negf%tRead_negf_in)) &
      .or.(negf%tElastic.and.(negf%tOrthonormal.or.negf%tOrthonormalDevice).and.(.not.negf%tReadDFTB)) &  
        ) then 
      
         if (id0) then
         !Writing H to file H_dftb.mtr      
         call writeSparseAsSquare_old('H_dftb.mtr', ham(:,iS)*27.21138469, iNeighbor, nNeighbor, &
                                       &iAtomStart, iPair, img2CentCell)
         if (negf%verbose.gt.30) write(*,"(' Hamiltonian is written to the file ',A)")trim('H_dftb.mtr')
         
         !Writing S to file S_dftb.mtr
         call writeSparseAsSquare_old('S_dftb.mtr', over, iNeighbor, nNeighbor, &
                                       &iAtomStart, iPair, img2CentCell)
         if (negf%verbose.gt.30) write(*,"(' Overlap is written to the file ',A)")trim('S_dftb.mtr')
         end if
         call mpifx_barrier(mpicomm) 
      end if

      if( (negf%tManyBody.and.(.not.negf%tReadDFTB).and.(.not.negf%tModel).and.(.not.negf%tRead_negf_in)) &
         .or. (negf%tElastic.and.(.not.negf%tReadDFTB).and.(negf%tOrthonormal.or.negf%tOrthonormalDevice)) ) then

         NumStates = negf%NumStates

         if (id0.and.negf%verbose.gt.50) write(*,"(' Hamiltonian is red from the file ',A)")trim('H_dftb.mtr')
     
         open(11,file='H_dftb.mtr',action="read")
         read(11,*);read(11,*);read(11,*);read(11,*);read(11,*)
         if(.not.allocated(negf%H_all)) allocate(negf%H_all(NumStates,NumStates))
         negf%H_all=0.0_dp
         do i=1,NumStates
            read(11,*)negf%H_all(i,1:NumStates)
         end do
         close(11)
       
         negf%H_all=negf%H_all/27.21138469
        
         if (id0.and.negf%verbose.gt.50) write(*,"(' Overlap is red from the file ',A)")trim('S_dftb.mtr')
         open(12,file='S_dftb.mtr',action="read")
         read(12,*);read(12,*);read(12,*);read(12,*);read(12,*)
         if(.not.allocated(negf%S_all)) allocate(negf%S_all(NumStates,NumStates))
         negf%S_all=0.0_dp
         do i=1,NumStates
            read(12,*)negf%S_all(i,1:NumStates)
         end do
         close(12)
         !debug begin
         !do i=1,1 !NumStates
         !   write(*,"(10000F9.3)")negf%H_all(i,1:10) !NumStates)
         !end do
         !debug end

      end if

      negf%mpicomm = mpicomm !DAR! for Orthogonalization_dev temporary solution
      
      !-------------------------------------------------------------------------
      !DAR end

      call negf_current(csrHam, csrOver, iS, iK, kWeights(iK), &
           tunnMat, ldosMat, currVec, tundos)
       
      !DAR - tundos is added,
      !      it is necessary for 'call negf_init_elph(tundos%elph)' in 'negf_current'

      !-------------------------------------------------------------------------   
      !DAR begin - WriteSelfEnergy / WriteSurfaceGF
      !-------------------------------------------------------------------------

      do icont=1,negf%str%num_conts
         if(negf%tranas%cont(icont)%tWriteSelfEnergy) &
              call mpifx_allreduceip(mpicomm, negf%tranas%cont(icont)%SelfEnergy, MPI_SUM)
      end do
       
    if (id0) then
       do icont=1,negf%str%num_conts
          if(negf%tranas%cont(icont)%tWriteSelfEnergy) then
            if (negf%tranas%cont(icont)%tUnformatted) then
              open(14,form="unformatted",file=trim(negf%tranas%cont(icont)%name)//'-SelfEnergy.mgf' &
                   ,action="write")
              do i = 1, size(negf%en_grid)
                write(14)real(negf%en_grid(i)%Ec),negf%tranas%cont(icont)%SelfEnergy(:,:,i)
              end do
            else
              open(14,file=trim(negf%tranas%cont(icont)%name)//'-SelfEnergy.mgf' &
                   ,action="write")
              do i = 1, size(negf%en_grid)
                write(14,*)real(negf%en_grid(i)%Ec),negf%tranas%cont(icont)%SelfEnergy(:,:,i)
              end do
            end if 
            close(14)
            if (negf%verbose.gt.30) write(*,"('The retarded contact self-energy is written into the file ',A)") &
                  trim(negf%tranas%cont(icont)%name)//'-SelfEnergy.mgf'
          end if
       end do
    end if
      
      do icont=1,negf%str%num_conts
         if(negf%tranas%cont(icont)%tWriteSurfaceGF) &
              call mpifx_allreduceip(mpicomm, negf%tranas%cont(icont)%SurfaceGF, MPI_SUM)
      end do
       
    if (id0) then
       do icont=1,negf%str%num_conts
          if(negf%tranas%cont(icont)%tWriteSurfaceGF) then
            if (negf%tranas%cont(icont)%tUnformatted) then
              open(14,form="unformatted",file=trim(negf%tranas%cont(icont)%name)//'-SurfaceGF.mgf' &
                   ,action="write")
              do i = 1, size(negf%en_grid)
                write(14)real(negf%en_grid(i)%Ec),negf%tranas%cont(icont)%SurfaceGF(:,:,i)
              end do
            else
              open(14,file=trim(negf%tranas%cont(icont)%name)//'-SurfaceGF.mgf' &
                   ,action="write")
              do i = 1, size(negf%en_grid)
                write(14,*)real(negf%en_grid(i)%Ec),negf%tranas%cont(icont)%SurfaceGF(:,:,i)
              end do
            end if 
            close(14)
            if (negf%verbose.gt.30) write(*,"('The retarded surface Green function is written into the file ',A)") &
                  trim(negf%tranas%cont(icont)%name)//'-SurfaceGF.mgf'
          end if
       end do
    end if
     
      !-------------------------------------------------------------------------
      !DAR end
      !-------------------------------------------------------------------------
      
      ! Adding up k and spin components:
      ! on first iteration allocate not allocated arrays

    if (tundos%writeTunn) then
      
      if(.not.allocated(currTot)) then
        allocate(currTot(size(currVec)), stat=err)
        if (err/=0) STOP 'Allocation error (currTot)'
        currTot = 0.0_dp
      endif
      currTot = currTot + currVec

      call add_partial_results(tunnMat, tunnTot, tunnSKRes, iKS, nKS) 
      
      if (associated(tunnMat).and.(nKS.gt.1)) then
        call mpifx_allreduceip(mpicomm, tunnSKRes(:,:,iKS), MPI_SUM)
      end if

    end if  

    if (tundos%writeLDOS) then
    
      call add_partial_results(ldosMat, ldosTot, ldosSKRes, iKS, nKS) 

      if (associated(ldosMat).and.(nKS.gt.1)) then
        call mpifx_allreduceip(mpicomm, ldosSKRes(:,:,iKS), MPI_SUM)
      end if

    end if  
   
  end do
  ! End iKS 
    
  if (tundos%writeTunn) then  
    
    call mpifx_allreduceip(mpicomm, currTot, MPI_SUM)
    
    ! converts from internal atomic units into A
    currTot = currTot * convertCurrent(unitOfEnergy, unitOfCurrent) 
    
    if (id0.and.negf%verbose.gt.0) then
      write(*,"(/,'Current is calculated')") 
      do ii=1, size(currTot)
        write(*,'(1x,a,i3,i3,a,ES14.5,a,a)') &
             & ' contacts: ',params%ni(ii),params%nf(ii), &
             & '   current: ', currTot(ii),' ',unitOfCurrent%name
      enddo
    endif
     
    if (allocated(tunnTot)) then
      !print*,'reduce tunnTot cpu #',id                                     !DAR
      call mpifx_allreduceip(mpicomm, tunnTot, MPI_SUM)
      ! Write Total transmission on a separate file (optional)
      if (id0 .and. writeTunn) then
        call write_file(negf, tunnTot, tunnSKRes, 'transmission', &
                         & groupKS, kpoints, kWeights)
        if (allocated(tunnSKRes)) deallocate(tunnSKRes)
      endif 
    else
      allocate(tunnTot(0,0))
    endif

    !DAR begin - print tunn_mat_bp
    if(negf%tZeroCurrent) then                                                
       call mpifx_allreduceip(mpicomm, negf%tunn_mat_bp, MPI_SUM)        
       if (id0) call write_file(negf, negf%tunn_mat_bp, tunnSKRes, 'transmission_bp', &
            & groupKS, kpoints, kWeights)
        if (allocated(tunnSKRes)) deallocate(tunnSKRes)
    endif
    !DAR end

  end if   

  if (tundos%writeLDOS) then
     
    if (allocated(ldosTot)) then
      !print*,'reduce ldosTot cpu #',id                                     
      call mpifx_allreduceip(mpicomm, ldosTot, MPI_SUM)
      ! Write Total localDOS on a separate file (optional)
      if (id0 .and. writeLDOS) then
        call write_file(negf, ldosTot, ldosSKRes, 'localDOS', &
                         & groupKS, kpoints, kWeights)
        if (allocated(ldosSKRes)) deallocate(ldosSKRes)
      end if
    else
      allocate(ldosTot(0,0))
    end if

  end if
    
  end subroutine calc_current

  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------
  !   utility to allocate and sum partial results
  !----------------------------------------------------------------------------
  subroutine add_partial_results(pMat, pTot, pSKRes, iK, nK)
    real(dp), pointer :: pMat(:,:)
    real(dp), allocatable :: pTot(:,:)
    real(dp), allocatable :: pSKRes(:,:,:)
    integer, intent(in) :: iK, nK
      
    integer :: err 

    if (associated(pMat)) then
      if(.not.allocated(pTot)) then 
        allocate(pTot(size(pMat,1), size(pMat,2)), stat=err)
        if (err/=0) STOP 'Allocation error (tunnTot)'
        pTot = 0.0_dp
      endif
      pTot = pTot + pMat
   
      if (nK.gt.1) then
        if(.not.allocated(pSKRes)) then 
          allocate(pSKRes(size(pMat,1), size(pMat,2), nK), stat=err)
          if (err/=0) STOP 'Allocation error (tunnSKRes)'
          pSKRes = 0.0_dp
        endif
        pSKRes(:,:,iK) = pMat(:,:)
      end if
    endif
  
  end subroutine add_partial_results
  ! ----------------------------------------------------------------------------
  !    utility to write transmission or ldos on files 
  ! ----------------------------------------------------------------------------
  subroutine write_file(negf, pTot, pSKRes, filename, groupKS, kpoints, kWeights)
    type(TNegf) :: negf
    real(dp), intent(in) :: pTot(:,:)
    real(dp), intent(in) :: pSKRes(:,:,:)
    character(*), intent(in) :: filename
    integer, intent(in) :: groupKS(:,:)
    real(dp), intent(in) :: kPoints(:,:)
    real(dp), intent(in) :: kWeights(:)

    integer :: ii, jj, nKS, iKS, nK, nS, iK, iS
    type(lnParams) :: params

    call get_params(negf, params)

    nKS = size(groupKS, dim=2) 
    nK = size(kpoints, dim=2)
    nS = nKS/nK 

    open(65000,file=trim(filename)//'.dat')
    do ii=1,size(pTot,1)
       !write(65000,'(f20.8)',ADVANCE='NO') (params%Emin+ii*params%Estep) * HAR 
       write(65000,'(E17.8)',ADVANCE='NO') (params%Emin+ii*params%Estep) * HAR !DAR
      do jj=1,size(pTot,2)
        !write(65000,'(es20.8)',ADVANCE='NO') pTot(ii,jj)
        write(65000,'(E17.8)',ADVANCE='NO') pTot(ii,jj)                     !DAR
      enddo
      write(65000,*)
    enddo
    close(65000)

    if (nKS.gt.1) then
      open(65000,file=trim(filename)//'_kpoints.dat')
      write(65000,*)  '# NKpoints = ', nK
      write(65000,*)  '# NSpin = ', nS
      write(65000,*)  '# Energy [eV], <spin k1 k2 k3 weight> '
      write(65000,'(A1)', ADVANCE='NO') '# '
      do iKS = 1,nKS
        iK = groupKS(1,iKS)
        iS = groupKS(2,iKS)    
        write(65000,'(i5.2)', ADVANCE='NO') iS
        write(65000,'(es15.5, es15.5, es15.5, es15.5)', ADVANCE='NO') kpoints(:,iK),&
                                                                      & kWeights(iK) 
      end do
      write(65000,*)
      do ii=1,size(pSKRes(:,:,1),1)
        write(65000,'(f20.8)',ADVANCE='NO') (params%Emin+ii*params%Estep) * HAR
        do jj=1,size(pSKRes(:,:,1),2)
          do iKS = 1,nKS
            write(65000,'(es20.8)',ADVANCE='NO') pSKRes(ii,jj, iKS)
          enddo
          write(65000,*)
        enddo
      enddo
      close(65000)
    end if

  end subroutine write_file

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! THIS is a first version of local current computation. 
  ! It has been placed here since it depends on internal representations of DFTB
  !       
  ! NOTE: Limited to non-periodic systems             s !!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine local_currents(mpicomm, groupKS, ham, over, &
      & iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, &
      & cellVec, orb, kPoints, kWeights, coord, dumpDens, chempot, tundos)
      !DAR - tundos is added,
      !      it is necessary for 'call negf_init_elph(tundos%elph)' in 'call negf_current'
    
    Type(TGDFTBTunDos), intent(IN) :: tundos                                !DAR
    type(mpifx_comm), intent(in) :: mpicomm
    integer, intent(in) :: groupKS(:,:)
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: kPoints(:,:), kWeights(:)
    real(dp), intent(in) :: coord(:,:)
    logical, intent(in) :: dumpDens
    !! We need this now for different fermi levels in colinear spin
    !! Note: the spin polirized does not work with
    !! built-int potentials (the unpolarized does) in the poisson
    !! I do not set the fermi because it seems that in libnegf it is 
    !! not really needed
    real(dp), intent(in) :: chempot(:,:)
    
    integer :: n,m, mu, nu, nAtom, irow, nrow, ncont
    integer :: nKS, nKPoint, nSpin, iK, iKS, iS, inn, startn, endn, morb
    real(dp), dimension(:), allocatable :: Im 
    integer, dimension(:,:), allocatable :: nneig
    integer, dimension(:), allocatable :: nn
    integer, parameter :: NMAX=40
    complex(dp) :: c1,c2
    character(1) :: sp
    integer :: iSCCiter=2
    type(z_CSR) :: csrDens, csrEDens
    type(lnParams) :: params

    call get_params(negf, params)

    if (id0.and.params%verbose.gt.30) then
      write(*,*)
      write(*,'(73("="))')
      write(*,*) '                   COMPUTING LOCAL CURRENTS          '
      write(*,'(73("="))') 
      write(*,*)
    endif

    nKS = size(groupKS, dim=2)
    nKPoint = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)
    nAtom = size(orb%nOrbAtom)
    
    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)
      
      ! We need to recompute Rho and RhoE .....
      call foldToCSR(csrHam, ham(:,iS), kPoints(:,1), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)
      call foldToCSR(csrOver, over, kPoints(:,1), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)

      call negf_density(iSCCIter, iS, iKS, csrHam, csrOver, chempot(:,iS), tundos, DensMat=csrDens)
      !DAR - tundos is added,
      !      it is necessary for 'call negf_init_elph(tundos%elph)' in 'negf_current')
      call negf_density(iSCCIter, iS, iKS, csrHam, csrOver, chempot(:,iS), tundos, EnMat=csrEDens)
      !DAR - tundos is added,
      !      it is necessary for 'call negf_init_elph(tundos%elph)' in 'negf_current')

      call mpifx_allreduceip(mpicomm, csrDens%nzval, MPI_SUM)
      
      call mpifx_allreduceip(mpicomm, csrEDens%nzval, MPI_SUM)
      ! Gather is done here 
      if (dumpDens) then
        open(65000, file='dens.csr',form='formatted')
        call print_mat(65000, csrDens, .true.)
        close(65000)
      end if
      
      call log_allocate(nneig,nAtom,NMAX)
      call log_allocate(nn,nAtom)
      call symmetrize_neiglist(nAtom,img2CentCell,iNeighbor,nNeighbor,coord&
          &,nneig,nn)
      
      
      if (iS .eq. 1) sp = 'u'
      if (iS .eq. 2) sp = 'd'
      open(207,file='lcurrents_'//sp//'.dat')

      do m = 1, nAtom

        call log_allocate(Im,nn(m))
        Im(:) = 0.0_dp

        mOrb = orb%nOrbAtom(m)
        iRow = iAtomStart(m)

        write(207,'(I5,3(F12.6),I3)',advance='NO') m, coord(:,m), nn(m) 

        do inn = 1, nn(m)
          n = nneig(m,inn)
          startn = iAtomStart(n)
          endn = startn + orb%nOrbAtom(n) - 1

          ! tracing orbitals of atoms  n  m
          ! More efficient without getel ?    
          do mu = iRow, iRow+mOrb-1
            do nu = startn, endn
              c1=conjg(getel(csrDens,mu,nu))
              c2=conjg(getel(csrEDens,mu,nu))
              Im(inn)=Im(inn)+aimag(getel(csrHam,mu,nu)*c1 - &
                  &getel(csrOver,mu,nu)*c2)
            enddo
          enddo
          ! pi-factor  comes from  Gn = rho * pi 
          write(207,'(I5,ES17.8)',advance='NO') n, 2.d0*params%g_spin*pi*eovh*Im(inn)

        enddo

        write(207,*)

        call log_deallocate(Im) 

      enddo
    enddo

    close(207)
    call destruct(csrDens)
    call destruct(csrEDens)
    call log_deallocate(nneig)
    call log_deallocate(nn)

  end subroutine local_currents

  !--------------------------------------------------------------
  ! Neighbor is non-symmetric: i->j  j>i
  ! The following routine symmetrizes the neighbor list 
  !
  subroutine symmetrize_neiglist(nAtom,img2CentCell,iNeighbor,nNeighbor,coord&
      &,nneig,nn)
    integer, intent(in) :: nAtom
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    real(dp), intent(in) :: coord(:,:)
    integer, dimension(:,:) :: nneig
    integer, dimension(:) :: nn    
    
    real(dp) :: dist, dr(3)
    integer :: m, n, inn, ii, jj, morb
    integer, parameter :: NMAX=40
    
    nn=0
    nneig=0

    do m = 1, nAtom
      do inn = 1, nNeighbor(m)
        if (inn.gt.NMAX) exit
        n = img2CentCell(iNeighbor(inn,m))
        if(nn(m).le.NMAX-1) then
          nn(m)=nn(m)+1 
          nneig(m,nn(m))=n
          ! sort by distances
          jj = nn(m)
          dr(:) = coord(:,m)-coord(:,nneig(m,jj))
          dist= dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
          do ii = nn(m)-1,1,-1
            dr(:) = coord(:,m)-coord(:,nneig(m,ii))
            if (dist.lt.dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)) then
              morb=nneig(m,ii)     
              nneig(m,ii)=nneig(m,jj)
              nneig(m,jj)=morb
              jj = jj - 1 
            endif
          enddo
          ! -----------------
        endif
        if(nn(n).le.NMAX-1) then
          nn(n)=nn(n)+1 
          nneig(n,nn(n))=m
          ! sort by distances
          jj = nn(n)
          dr(:) = coord(:,n)-coord(:,nneig(n,jj))
          dist= dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3) 
          do ii = nn(n)-1,1,-1
            dr(:) = coord(:,n)-coord(:,nneig(n,ii))
            if (dist.lt.dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)) then
              morb=nneig(n,ii)     
              nneig(n,ii)=nneig(n,jj)
              nneig(n,jj)=morb
              jj=jj-1
            endif
          enddo
          ! -----------------
        endif
      enddo
    enddo
    !--------------------------------------------------------------
    !nn = 36 
    !print*,'nneig:'
    !do m =1,nAtom
    !   write(*,'(8(i6))',advance='NO') m
    !   do jj=1,nn(m)
    !     dr(:) = coord(:,m)-coord(:,nneig(m,jj))
    !     dist= sqrt(dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)) 
    !     write(*,'(i5,ES15.6)',advance='NO') nneig(m,nn(jj)),dist
    !   enddo 
    !   write(*,*) 
    !enddo
    !print*,'done'
    !--------------------------------------------------------------

  end subroutine symmetrize_neiglist

  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!

  subroutine ReadDFTB

    use dftbp_unitconversion
    
    integer :: i,j,k,l,m,n,NumStates

    NumStates = negf%NumStates

    if (id0.and.negf%verbose.gt.50) write(*,"(' Hamiltonian is red from the file ',A)")trim('H_dftb.mtr')
    if(.not.allocated(negf%H_all)) allocate(negf%H_all(NumStates,NumStates))
    open(11,file='H_dftb.mtr',action="read")
    read(11,*);read(11,*);read(11,*);read(11,*);read(11,*)
    negf%H_all=0.0_dp 
    do i=1,NumStates
       read(11,*)negf%H_all(i,1:NumStates)
    end do 
    close(11)
    if (id0.and.negf%verbose.gt.50) write(*,"(' Overlap is red from the file ',A)")trim('S_dftb.mtr')
    open(12,file='S_dftb.mtr',action="read")
    read(12,*);read(12,*);read(12,*);read(12,*);read(12,*)
    if(.not.allocated(negf%S_all)) allocate(negf%S_all(NumStates,NumStates))
    negf%S_all=0.0_dp
    do i=1,NumStates
      read(12,*)negf%S_all(i,1:NumStates)
    end do
    close(12)
    !debug begin
    if (id0.and.negf%verbose.gt.100) then
    do i=1,NumStates
       write(*,"(10000F9.3)")negf%H_all(i,1:NumStates)
    end do
    end if   
    !debug end
    call convertByMul_NoNode(negf%units_energy, energyUnits, negf%H_all)

  end subroutine ReadDFTB

  !------------------------------------------------------------------------------------------------!!

  subroutine ReadModel
                                                     
    use dftbp_unitconversion
    
    integer :: i,j,k,l,m,n,NumStates

    NumStates = negf%NumStates

    if (id0.and.negf%verbose.gt.50) write(*,"(' Hamiltonian is red from the file ',A)")trim(negf%FileName)
    open(11,file=negf%FileName,action="read")
    if(.not.allocated(negf%H_all)) allocate(negf%H_all(NumStates,NumStates))
    negf%H_all=0.0_dp
    do i=1,NumStates
       read(11,*)negf%H_all(i,1:NumStates)
    end do
    close(11)
    !extended output
    if (id0.and.negf%verbose.gt.100) then
    do i=1,NumStates
       write(*,"(10000F9.3)")negf%H_all(i,1:NumStates)
    end do
    end if
    !end
    call convertByMul_NoNode(negf%units_energy, energyUnits, negf%H_all)

    if (negf%tReadOverlap) then
      if (id0.and.negf%verbose.gt.50) write(*,"(' Overlap is red from the file ',A)")trim(negf%OverlapFile)
      open(11,file=negf%OverlapFile,action="read")
      if(.not.allocated(negf%S_all)) allocate(negf%S_all(NumStates,NumStates))
      negf%S_all=0.0_dp
      do i=1,NumStates
        read(11,*)negf%S_all(i,1:NumStates)
      end do
      close(11)
    else   
      if (.not.allocated(negf%S_all)) allocate(negf%S_all(NumStates,NumStates))
      negf%S_all=0.0_dp
      do i=1,NumStates
        negf%S_all(i,i)=1._dp
      end do
    end if
      
  end subroutine ReadModel

  !------------------------------------------------------------------------------------------------!
  !> Implementation of convertByMul for real scalar.
  subroutine convertByMul_NoNode(modifier, units, convertValue)
                                                           
    use dftbp_unitconversion

    !> Modifier (name of the unit to use)
    character(len=*), intent(in) :: modifier

    !> Array of the possible units
    type(unit), intent(in) :: units(:)

    !> Value to convert, converted value on return.
    real(dp), intent(inout) :: convertValue(:,:)
    
    integer :: ind   ,i

    if (len(modifier) > 0) then
      ind = getModifierIndex_NoNode(modifier, units)    
      convertValue = convertValue * units(ind)%convertValue
    end if

  end subroutine convertByMul_NoNode

  !> Gets the index of a modifier from an array of possible modifier names.
  function getModifierIndex_NoNode(modifier, modifiers) result(ind)
                                                       
    use dftbp_charmanip                                                            
    use dftbp_unitconversion
    
    !> String containing the parsed modifier
    character(len=*), intent(in) :: modifier

    !> Array containing the names of the possible modifiers
    type(unit), intent(in) :: modifiers(:)

    !> Index of the modifer (zero if not found)
    integer :: ind

    character(len=len(modifier)) :: modifierLo
    logical :: mandatory
    integer :: ii

    modifierLo = tolower(modifier)

    ind = 0
    do ii = 1, size(modifiers)
      if (trim(modifiers(ii)%name) == modifierLo) then
        ind = ii
        exit
      end if
    end do

  end function getModifierIndex_NoNode

  !------------------------------------------------------------------------------------------------!

  subroutine ReadU
      
    integer :: i,j,k,l,m,n,NumStates

    NumStates = negf%str%central_dim

    if (id0.and.negf%verbose.gt.50) write(*,"(' Hubbard Us are red from the file ',A)")trim('U.mtr')
    open(11,file='U.mtr',action="read")
    if(.not.allocated(negf%U)) allocate(negf%U(NumStates,NumStates))
    negf%U=0.0_dp
    do i=1,NumStates
       read(11,*)negf%U(i,1:NumStates)
    end do
    close(11)
    !extended output
    if (id0.and.negf%verbose.gt.100) then
    do i=1,NumStates
       write(*,"(10000F9.3)")negf%U(i,1:NumStates)
    end do
    end if
    !end
    !call convertByMul_NoNode(negf%units_energy, energyUnits, negf%U)
    negf%U=negf%U/27.21138469
    
  end subroutine ReadU

  !------------------------------------------------------------------------------------------------!

  subroutine ReadLibNEGF
   
    write(*,*) 'Import from the negf.in file'
    call read_negf_in(negf)
    call negf_partition_info(negf)

  end subroutine ReadLibNEGF

  !------------------------------------------------------------------------------------------------!

  subroutine MakeHHSS(HH,SS)

    type(z_CSR), intent(inout) :: HH, SS
    integer :: i,j,k,l,m,n,NumStates

    NumStates = negf%NumStates
    
    !open(12,file='H_real_0.dat',action="write")
    !open(13,file='H_imm_0.dat',action="write")
    HH%nnz=0
    HH%nrow=NumStates
    HH%ncol=NumStates
    do i=1,NumStates
    do j=1,NumStates   
       if((i.eq.j).or.(abs(negf%H_all(i,j)).gt.0.00001)) then
          HH%nnz=HH%nnz+1
    !      write(12,*)i,j,H(i,j)
    !      write(13,*)i,j,0.
       end if
    end do
    end do
    !close(12)
    !close(13)

    if(allocated(HH%nzval)) deallocate(HH%nzval)
    if(allocated(HH%colind)) deallocate(HH%colind)
    if(allocated(HH%rowpnt)) deallocate(HH%rowpnt)
    allocate(HH%nzval(HH%nnz))
    allocate(HH%colind(HH%nnz))
    allocate(HH%rowpnt(HH%nrow+1))
    if(allocated(SS%nzval)) deallocate(SS%nzval)
    if(allocated(SS%colind)) deallocate(SS%colind)
    if(allocated(SS%rowpnt)) deallocate(SS%rowpnt)
    allocate(SS%nzval(HH%nnz))
    allocate(SS%colind(HH%nnz))
    allocate(SS%rowpnt(HH%nrow+1))

    HH%nnz=0
    HH%rowpnt(1)=1
    SS%nnz=0
    do i=1,NumStates
       k=0
       do j=1,NumStates   
          if((i.eq.j).or.(abs(negf%H_all(i,j)).gt.0.00001)) then
             k=k+1
             HH%nnz=HH%nnz+1            
             if(i.eq.j) then
                HH%nzval(HH%nnz)=negf%H_all(i,j)
                SS%nzval(HH%nnz)=negf%S_all(i,j)
             else
                HH%nzval(HH%nnz)=negf%H_all(i,j)
                SS%nzval(HH%nnz)=negf%S_all(i,j)
             end if
             HH%colind(HH%nnz)=j
             !print *, HH%nnz, negf%H_all(i,j),  HH%nzval(HH%nnz), SS%nzval(HH%nnz)   !debug
          end if
       end do
       HH%rowpnt(i+1)=HH%rowpnt(i)+k
    end do

    SS%nnz=HH%nnz
    SS%ncol=HH%ncol
    SS%nrow=HH%nrow
    SS%colind=HH%colind
    SS%rowpnt=HH%rowpnt

  end subroutine MakeHHSS
  
  !------------------------------------------------------------------------------------------------!
  
  subroutine MakeHS_dev

    allocate(negf%H_dev(negf%str%central_dim,negf%str%central_dim))
    negf%H_dev=0.0_dp
    negf%H_dev=negf%H_all(1:negf%str%central_dim,1:negf%str%central_dim)
    allocate(negf%S_dev(negf%str%central_dim,negf%str%central_dim))
    negf%S_dev=0.0_dp
    negf%S_dev=negf%S_all(1:negf%str%central_dim,1:negf%str%central_dim)

  end subroutine MakeHS_dev
 
  !------------------------------------------------------------------------------------------------!
  
  subroutine check_negf_params

    use globals, only : LST
    
    Integer :: ncont, nbl, ii, jj, ist, iend
    Integer, dimension(:), allocatable :: PL_end, cont_end, surf_end
    character(32) :: tmp
    character(LST) :: file_re_H, file_im_H, file_re_S, file_im_S
    integer, dimension(:), allocatable :: cblk

    open(10,file='log',action="write")
    
    write(10,"('negf%H%nnz               = ',I6)")negf%H%nnz
    write(10,"('negf%H%nrow              = ',I6)")negf%H%nrow
    write(10,"('negf%H%ncol              = ',I6)")negf%H%ncol
    write(10,"('negf%H%nzval             = ',10000000F8.4)")negf%H%nzval
    write(10,"('negf%H%colind            = ',10000000I6)")negf%H%colind
    write(10,"('negf%H%rowpnt            = ',10000000I6)")negf%H%rowpnt
    write(10,"('negf%isSid               = ',L6)")negf%isSid
    write(10,"('negf%S%nnz               = ',I6)")negf%S%nnz
    write(10,"('negf%S%nrow              = ',I6)")negf%S%nrow
    write(10,"('negf%S%ncol              = ',I6)")negf%S%ncol
    write(10,"('negf%S%nzval             = ',10000000F8.4)")negf%S%nzval
    write(10,"('negf%S%colind            = ',10000000I6)")negf%S%colind
    write(10,"('negf%S%rowpnt            = ',10000000I6)")negf%S%rowpnt
    write(10,"('negf%str%num_conts       = ',I6)")negf%str%num_conts
    write(10,"('negf%str%num_PLs         = ',I6)")negf%str%num_PLs
    write(10,"('negf%str%active_cont     = ',I6)")negf%str%active_cont
    write(10,"('negf%str%mat_B_start     = ',I6,I6)")negf%str%mat_B_start
    write(10,"('negf%str%mat_C_start     = ',I6,I6)")negf%str%mat_C_start
    write(10,"('negf%str%mat_C_end       = ',I6,I6)")negf%str%mat_C_end
    write(10,"('negf%str%cblk            = ',I6,I6)")negf%str%cblk
    write(10,"('negf%str%cont_dim        = ',I6,I6)")negf%str%cont_dim
    write(10,"('negf%str%mat_PL_start    = ',10000000I6)")negf%str%mat_PL_start
    write(10,"('negf%str%mat_PL_end      = ',10000000I6)")negf%str%mat_PL_end
    write(10,"('negf%str%central_dim     = ',I6)")negf%str%central_dim
    write(10,"('negf%str%total_dim       = ',I6)")negf%str%total_dim
    write(10,"('negf%Ec, negf%Ev         = ',3F8.4)")negf%Ec, negf%Ev
    write(10,"('negf%DeltaEc, %DeltaEv   = ',3F8.4)")negf%DeltaEc, negf%DeltaEv
    write(10,"('negf%Emin, %Emax, %Estep = ',3F8.4)")negf%Emin, negf%Emax, negf%Estep
    write(10,"('negf%kbT                 = ',10000000E12.4)")negf%kbT
    write(10,"('negf%mu_n                = ',10000000F8.4)")negf%mu_n
    write(10,"('negf%mu_p                = ',10000000F8.4)")negf%mu_p
    write(10,"('negf%mu                  = ',10000000F8.4)")negf%mu
    write(10,"('negf%delta               = ',10000000E12.4)")negf%delta

    write(10,"('negf%wght                = ',10000000F8.4)")negf%wght
    write(10,"('negf%Np_n(1:2)           = ',I6,I6)")negf%Np_n(1:2)
    write(10,"('negf%Np_p(1:2)           = ',I6,I6)")negf%Np_p(1:2)
    write(10,"('negf%n_poles              = ',I6)")negf%n_poles
    write(10,"('negf%Np_real(1)           = ',I6)")negf%Np_real(1)
    write(10,"('negf%n_kt                 = ',I6)")negf%n_kt
    write(10,"('negf%g_spin               = ',F8.4)")negf%g_spin
    write(10,"('negf%nLDOS                = ',I6)")negf%nLDOS
    write(10,"('negf%LDOS(1)%indexes      = ',10000000I6)")negf%LDOS(1)%indexes
    write(10,"('negf%verbose              = ',I6)")negf%verbose

    close (10)

    write(*,*)
    write(*,"(' negf parameters are written to the log file')")
        
  end subroutine check_negf_params
  
  !------------------------------------------------------------------------------------------------!

  subroutine orthogonalization

    integer :: i,m,n1_first,n1_last,n2_first,n2_last
    integer :: INFO, N
    real(dp), allocatable :: A(:,:), WORK(:), W(:)
    real(dp), allocatable :: B(:,:),C(:,:)
    character(mc) :: strForm

    if (id0.and.negf%verbose.gt.50) write(*,"(' Löwdin orthogonalization is started')")

    N=negf%NumStates

    allocate(A(N,N),WORK(3*N),W(N))
    allocate(B(N,N),C(N,N))
    W=0.0_dp
    WORK=0.0_dp
    A=0.0_dp
    B=0.0_dp
    C=0.0_dp
    
    A=negf%S_all

    call DSYEV('V','U',N,A,N,W,WORK,3*N,INFO )

    !print  *,'U matrix, Eigenvectors for S diagonalization'
    !do i=1,N
    !   write(*,*)A(i,1:N)
    !end do

    B=matmul(transpose(A),A)
    !print *,'U matrix unitarity check'
    !do i=1,N
    !   write(*,*)B(i,1:N)
    !end do

    B=matmul(transpose(A),matmul(negf%S_all,A))
    !print *,'S diagonal (transformed)'
    !do i=1,N
    !   write(*,*)B(i,1:N)
    !end do

    !B=sqrt(Sdiag)inverted
    do i=1,N
       B(i,i)=1./sqrt(B(i,i))
    end do

    !C=sqrt(S)
    C=matmul(A,matmul(B,transpose(A)))
    !print *,'sqrt(S) inverted'
    !do i=1,N
    !   write(*,*)C(i,1:N)
    !end do

    B=matmul(transpose(C),matmul(negf%S_all,C))
    !print *,'S unity check'
    !do i=1,N
    !   write(*,*)B(i,1:N)
    !end do

    !print *,'H_dftb before orthogonalization'
    !do i=1,N
    !   write(*,*)negf%H_all(i,1:N)
    !end do

    negf%H_all=matmul(transpose(C),matmul(negf%H_all,C))

    !print *,'H_dftb_orth before replacement'
    !do i=1,N
    !   write(*,*)negf%H_all(i,1:N)
    !end do

    do m=1,negf%str%num_conts
       n1_first=negf%str%mat_B_start(m)
       n1_last=(negf%str%mat_C_end(m)+negf%str%mat_B_start(m))/2
       n2_first=n1_last+1
       n2_last=negf%str%mat_C_end(m)
       !print *,n1_first,n1_last,n2_first,n2_last
       negf%H_all(n2_first:n2_last,n2_first:n2_last)=negf%H_all(n1_first:n1_last,n1_first:n1_last)
    end do

    !print *,'H_dftb_orth after replacement'
    !do i=1,N
    !   write(*,*)negf%H_all(i,1:N)
    !end do

    negf%S_all=0.0_dp
    do i=1,N
       negf%S_all(i,i)=1.0_dp
    end do

    write (strForm, "(A,I0,A)") "(", N, "ES24.15)"
    if (id0) then
      !Save H_dftb_orth.mtr to file
      open(12,file='H_dftb_orth.mtr',action="write")
      do i=1,N
        write(12, strForm)negf%H_all(i,1:N)*27.21138469
      end do
      close(12)
    end if  

    negf%delta=0.0001
    if (id0.and.negf%verbose.gt.50) write(*,"(' Hamiltonian is written to the file ',A)")trim('H_dftb_orth.mtr')
    if (id0.and.negf%verbose.gt.50) write(*,"(' Löwdin orthogonalization is done! ')")
    
  end subroutine orthogonalization

  !------------------------------------------------------------------------------------------------!

  subroutine orthogonalization_dev

    integer :: i,m,n1_first,n1_last,n2_first,n2_last
    integer :: INFO, N, N2
    real(dp), allocatable :: A(:,:), WORK(:), W(:)
    real(dp), allocatable :: B(:,:),C(:,:)
    character(mc) :: strForm

    if (id0.and.negf%verbose.gt.50) write(*,"(' Löwdin orthogonalization for device only is started')")

    N=negf%NumStates
    N2=negf%str%central_dim

    allocate(A(N2,N2),WORK(3*N2),W(N2))
    allocate(B(N2,N2),C(N2,N2))
    W=0.0_dp
    WORK=0.0_dp
    A=0.0_dp
    B=0.0_dp
    C=0.0_dp
   
    A=negf%S_all(1:N2,1:N2)

    call DSYEV('V','U',N2,A,N2,W,WORK,3*N2,INFO )

    !print  *,'U matrix, Eigenvectors for S diagonalization'
    !do i=1,N2
    !   write(*,*)A(i,1:N2)
    !end do

    B=matmul(transpose(A),A)
    !print *,'U matrix unitarity check'
    !do i=1,N2
    !   write(*,*)B(i,1:N2)
    !end do

    B=matmul(transpose(A),matmul(negf%S_all(1:N2,1:N2),A))
    !print *,'S diagonal (transformed)'
    !do i=1,N2
    !   write(*,*)B(i,1:N2)
    !end do

    !B=sqrt(Sdiag)inverted
    do i=1,N2
       B(i,i)=1./sqrt(B(i,i))
    end do

    !C=sqrt(S)
    C=matmul(A,matmul(B,transpose(A)))
    !print *,'sqrt(S) inverted'
    !do i=1,N2
    !   write(*,*)C(i,1:N2)
    !end do

    !S unity check
    B=matmul(transpose(C),matmul(negf%S_all(1:N2,1:N2),C))
    !print *,'S unity check'
    !do i=1,N2
    !   write(*,*)B(i,1:N2)
    !end do

    B=C
    deallocate(C)
    allocate(C(N,N))
    C=0._dp
    do i=N2+1,N
    C(i,i)=1._dp
    end do
    C(1:N2,1:N2)=B

    !C=sqrt(S) big matrix
    !print *,'C=sqrt(S) big matrix'
    !do i=1,N
    !   write(*,*)C(i,1:N)
    !end do
        
    !print *,'H_dftb before orthogonalization'
    !do i=1,1 !N
    !   write(*,*)negf%H_all(i,1:10)
    !end do

    negf%H_all=matmul(transpose(C),matmul(negf%H_all,C))
    negf%S_all=matmul(transpose(C),matmul(negf%S_all,C))

    !print *,'H_dftb_orth before replacement'
    !do i=1,1 !N
    !   write(*,*)negf%H_all(i,1:10)
    !end do

    do m=1,negf%str%num_conts
       n1_first=negf%str%mat_B_start(m)
       n1_last=(negf%str%mat_C_end(m)+negf%str%mat_B_start(m))/2
       n2_first=n1_last+1
       n2_last=negf%str%mat_C_end(m)
       !print *,n1_first,n1_last,n2_first,n2_last
       negf%H_all(n2_first:n2_last,n2_first:n2_last)=negf%H_all(n1_first:n1_last,n1_first:n1_last)
       negf%S_all(n2_first:n2_last,n2_first:n2_last)=negf%S_all(n1_first:n1_last,n1_first:n1_last)
    end do

    !print *,'H_dftb_orth after replacement'
    !do i=1,1 !N
    !   write(*,*)negf%H_all(i,1:10)
    !end do

    !print *,'S_dftb_orth after replacement'
    !do i=1,N
    !   write(*,*)negf%S_all(i,1:N)
    !end do

    write (strForm, "(A,I0,A)") "(", N, "ES24.15)"
    if (id0) then 
      !Save H_dftb_orth.mtr to file
      open(12,file='H_dftb_orth_dev.mtr',action="write")
      do i=1,N
        write(12, strForm)negf%H_all(i,1:N)*27.21138469
      end do
      close(12)
      !Save S_dftb_orth.mtr to file
      open(12,file='S_dftb_orth_dev.mtr',action="write")
      do i=1,N
        write(12, strForm)negf%S_all(i,1:N)
      end do
      close(12)
    end if
   
    call mpifx_barrier(negf%mpicomm)   

    negf%delta=0.000001
    
    if (id0.and.negf%verbose.gt.50) write(*,"(' Hamiltonian is written to the file ',A)")trim('H_dftb_orth_dev.mtr')
    if (id0.and.negf%verbose.gt.50) write(*,"(' Overlap is written to the file ',A)")trim('S_dftb_orth_dev.mtr')
    if (id0.and.negf%verbose.gt.50) write(*,"(' Löwdin orthogonalization for device only is done! ')")
    
  end subroutine orthogonalization_dev
  
  !------------------------------------------------------------------------------------------------!

  function inv(A) result(Ainv)
    
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), dimension(size(A,1),size(A,2)) :: Ainv

  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
  
  end function inv
  
end module tranas_interface

