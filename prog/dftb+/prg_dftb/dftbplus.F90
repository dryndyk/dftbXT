!--------------------------------------------------------------------------------------------------!
! DFTB+XT open software package for quantum nanoscale modeling (TraNaS OpenSuite)                  !
! Copyright (C) 2018-2019 Dmitry A. Ryndyk                                                         !
! DFTB+: general package for performing fast atomistic simulations                                 !
! Copyright (C) 2006-2019 DFTB+ developers group                                                   !
!--------------------------------------------------------------------------------------------------!
! GNU Lesser General Public License version 3 or (at your option) any later version.               !
! See the LICENSE file for terms of usage and distribution.                                        !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program dftbplus
  use dftbp_globalenv
  use dftbp_environment
  use dftbp_main, only : runDftbPlus
  use dftbp_inputdata_module, only : inputData
  use dftbp_formatout, only : printDftbHeader
  use dftbp_hsdhelpers, only : parseHsdInput
  use dftbp_initprogram, only : initProgramVariables, destructProgramVariables
#:if WITH_TRANSPORT
  use dftbp_initprogram, only : initProgramVariables, &
                          negf_init_nogeom, negf_init_str, tranasNoGeom !DAR
  use libmpifx_module !DAR
  use tranas_vars !DAR
  use dftbp_periodic !DAR
#:else
#:endif
  use dftbp_densedescr
  
  implicit none

  character(len=*), parameter :: releaseName = '${RELEASE}$'
  integer, parameter :: releaseYear = 2019

  type(TEnvironment) :: env
  type(inputData), allocatable :: input
#:if WITH_TRANSPORT  
  integer :: iAtom, nAtom !DAR
  logical :: tInitNEGF !DAR
  type(TDenseDescr) :: denseDescr !DAR
#:endif  
  integer, allocatable :: nNeigh(:)
  integer, allocatable :: img2CentCell(:)
  integer, allocatable :: iNeigh(:,:)

  call initGlobalEnv()
  call printDftbHeader(releaseName, releaseYear, stdOut)
  allocate(input)
  call parseHsdInput(input)
  call TEnvironment_init(env)
  
  !------------------------------------------------------------------------------------------------!
  !DAR begin - NoGeometry
  !------------------------------------------------------------------------------------------------!
#:if WITH_TRANSPORT
  if(input%transpar%tNoGeometry) then
    call env%initMpi(1)    
    call mpifx_barrier(env%mpi%globalComm)    
    if (input%ctrl%verbose.gt.0) then
      write(stdout, "(/,A)") repeat("-", 80)
      write(stdOut, "(A)") "-- Initialization is started (without geometry)                               --"
      write(stdout, "(A)") repeat("-", 80)
    end if
      
    input%transpar%tWriteTagged = input%ctrl%tWriteTagged
    
    call negf_init_nogeom(input%transpar,input%ginfo%greendens,input%ginfo%tundos,env%mpi%globalComm,tInitNEGF)
    if (.not.tInitNEGF) write(stdout, "('TraNaS Lib initialization error (negf_init_nogeom)')")
    nAtom=input%transpar%NumStates
    allocate(denseDescr%iAtomStart(nAtom+1))
    do iAtom=1,nAtom+1          
      denseDescr%iAtomStart(iAtom)=iAtom
    end do
    call negf_init_str(denseDescr, input%transpar, input%ginfo%greendens, &
         iNeigh, nNeigh, img2CentCell)
    if (input%ctrl%verbose.gt.0) then
      write(stdout, "(/,A)") repeat("-", 80)
      write(stdOut, "(A)") "-- Initialization is finished                                                 --"
      write(stdout, "(A)") repeat("-", 80)
    end if
    !call destroy(input)    
    call tranasNoGeom(env%mpi%globalComm,input%ginfo%tundos, input%transpar)
    !DAR - input%ginfo%tundos is added,
    !      it is necessary for 'call negf_init_elph(tundos%elph)' in tranas_interface_nogeom)
    deallocate(input)
  call destructProgramVariables()
#:if WITH_GPU  
  call magmaf_finalize()
#:endif
    call env%destruct()
    call destructGlobalEnv()
    
  else
 
  !------------------------------------------------------------------------------------------------!
  !DAR end
  !------------------------------------------------------------------------------------------------!
#:endif  
  call initProgramVariables(input, env)
  !!DAR!! deallocate(input)      !! Big temporary hack.
  call runDftbPlus(env, input)   !!DAR!! + input
  call env%destruct()
  call destructGlobalEnv()
#:if WITH_TRANSPORT
  end if
#:endif
  write(stdout, "(/,A)") repeat("=", 80)
  write(stdout,"('DFTB+XT is finished')")
  write(stdout, "(A,/)") repeat("=", 80)
end program dftbplus
