!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!--------------------------------------------------------------------------------------------------!
!  DFTB+XT open software package for quantum nanoscale modeling                                    !
!  Copyright (C) 2018 Dmitry A. Ryndyk.                                                            !
!--------------------------------------------------------------------------------------------------!
!  GNU Lesser General Public License version 3 or (at your option) any later version.              !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program dftbplus
  use globalenv
  use environment
  use main, only : runDftbPlus
  use inputdata_module, only : inputData
  use formatout, only : printDftbHeader
  use parser, only : parseHsdInput
  use initprogram, only : initProgramVariables, &
                          negf_init_nogeom, negf_init_str, negf_current_nogeom !DAR
  use libmpifx_module !DAR
  use libnegf_vars !DAR
  use periodic !DAR
  implicit none

  character(len=*), parameter :: releaseName = '${RELEASE}$'
  integer, parameter :: releaseYear = 2018

  type(TEnvironment) :: env
  type(inputData), allocatable :: input
  integer :: iAtom !DAR
  logical :: tInitNEGF !DAR
  Type(TGDFTBstructure) :: gdftbStr !DAR
  integer, allocatable :: nNeigh(:)
  integer, allocatable :: img2CentCell(:)
  integer, allocatable :: iNeigh(:,:)

  call initGlobalEnv()
  call printDftbHeader(releaseName, releaseYear, stdOut)
  allocate(input)
  call parseHsdInput(input)
  call TEnvironment_init(env)
 
  call env%initMpi(1) !!DAR
 
  !------------------------------------------------------------------------------------------------!
  !DAR begin - NoGeometry
  !------------------------------------------------------------------------------------------------!

  if(input%transpar%tNoGeometry) then
    call mpifx_barrier(env%mpi%globalComm) 
    write(stdout, "(A)") repeat("-", 80)
    write(stdOut, "(A)") "-- Initialization is started (without geometry)                               --"
    write(stdout, "(A)") repeat("-", 80)

    input%transpar%tWriteTagged = input%ctrl%tWriteTagged
    
    call negf_init_nogeom(input%transpar,input%ginfo%greendens,input%ginfo%tundos,env%mpi%globalComm,tInitNEGF)
    if (.not.tInitNEGF) write(stdout, "('libnegf initialization error (negf_init_nogeom)')")
    gdftbStr%nAtom=input%transpar%NumStates
    allocate(gdftbStr%iAtomStart(gdftbStr%nAtom+1))
    do iAtom=1,gdftbStr%nAtom+1          
      gdftbStr%iAtomStart(iAtom)=iAtom
    end do
    call negf_init_str(gdftbStr, input%transpar, input%ginfo%greendens, &
         iNeigh, nNeigh, img2CentCell)     
    write(stdout, "(A)") repeat("-", 80)
    write(stdOut, "(A)") "-- Initialization is finished                                                 --"
    write(stdout, "(A)") repeat("-", 80)
    !call destroy(input)
    call negf_current_nogeom(env%mpi%globalComm,input%ginfo%tundos)
    !DAR - input%ginfo%tundos is added,
    !      it is necessary for 'call negf_init_elph(tundos%elph)' in negf_int_nogeom)
    write(stdout,*)
      
    !call myclock_full%stop()
    !call myclock_full%print_times(msg="program total", fp=stdout)
    deallocate(input)
    call env%destruct()
    call destructGlobalEnv()
    
  else
 
  !------------------------------------------------------------------------------------------------!
  !DAR end
  !------------------------------------------------------------------------------------------------!
  
  call initProgramVariables(input, env)
  !!DAR deallocate(input)        !! Big temporary hack.
  call runDftbPlus(env, input)   !!DAR + input
  call env%destruct()
  call destructGlobalEnv()

  end if

end program dftbplus
