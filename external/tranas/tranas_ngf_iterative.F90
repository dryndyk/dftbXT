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

module tranas_ngf_iterative

  use ln_precision
  use ln_constants, only : pi
  use ln_allocation
  use mat_def
  use sparsekit_drv
  use inversions
  use elph
  use ln_structure, only : TStruct_Info
  use tranas_types_main
  use tranas_types_mbngf, only : TMBNGF
  use mpi_globals, only : id, numprocs, id0                                 
  use outmatrix, only : outmat_c, inmat_c, direct_out_c, direct_in_c 
  use clock
  
  implicit none
  private

  public :: transmission_dns
  public :: tun_and_dos

  public :: iterativeGreenRetarded
  public :: iterativeGreenLesser_Landauer
  public :: iterativeGreenLesser
  public :: iterativeMeirWingreen

  !public :: Make_gsmr_mem_dns
  !public :: Make_gsml_mem_dns
  !public :: Make_Gr_mem_dns
  !public :: Make_Grcol_mem_dns
  !public :: Make_Gn_mem_dns
  !public :: Outer_Gr_mem_dns
  !public :: Outer_Gn_mem_dns
  !public :: Outer_A_mem_dns
  
  public :: destroy_blk                                                     

  LOGICAL, PARAMETER :: debug=.false. 

  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsmr
  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsml
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Gr

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

  SUBROUTINE transmission_dns(H,S,Ec,SelfEneR,ni,nf,size_ni,str,TUN_MAT)

    implicit none

    Type(z_CSR) :: H
    Type(z_CSR) :: S           
    Complex(dp) :: Ec
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    !Type(z_DNS) :: SelfEner_d 
    Real(dp), Dimension(:) :: TUN_MAT
    Type(z_CSR) :: ESH_tot
    Type(TStruct_Info) :: str

    ! Local variables
    Type(z_DNS), Dimension(:,:), allocatable :: ESH
    Real(dp) :: tun
    Integer :: ni(MAXNCONT)
    Integer :: nf(MAXNCONT)
    Integer :: nbl,ncont,size_ni
    Integer :: i, ierr, icpl, nit, nft, nt, nt1

    nbl = str%num_PLs
    ncont = str%num_conts

    ! Create ESH : Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(H,S,(-1.d0, 0.d0),Ec,ESH_tot)    
    call allocate_blk_dns(ESH,nbl)
    call zcsr2blk_sod(ESH_tot,ESH,str%mat_PL_start)
    call destroy(ESH_tot)

    ! Add the retarded contact self-energies to the relevant blocks
    do i=1,ncont
      ESH(str%cblk(i),str%cblk(i))%val = ESH(str%cblk(i),str%cblk(i))%val-SelfEneR(i)%val
    enddo

    !--------------------------------------------------------------------------!
    ! Retarded Green function calculation
    !--------------------------------------------------------------------------!

    call allocate_gsm_dns(gsmr,nbl)
    call Make_gsmr_mem_dns(ESH,nbl,2) 

    !Computation of transmission(s) between contacts ni(:) -> nf(:)  
    nit=ni(1)
    nft=nf(1)

    !Arrange contacts in a way that the order between first and second is always the
    !same (always ct1 < ct2), so nt is the largest contact block 

    if (str%cblk(nit).gt.str%cblk(nft)) then
      nt = str%cblk(nit)
    else
      nt = str%cblk(nft)
    endif

    ! Iterative calculation of Gr down to nt1 
    call allocate_blk_dns(Gr,nbl)
    call Make_Gr_mem_dns(ESH,1)
    if (nt.gt.1) call Make_Gr_mem_dns(ESH,2,nt)

    if (ncont.eq.2) then
      call trasmission_dns(nit,nft,ESH,SelfEneR,str%cblk,tun) 
    else
      call trasmission_old(nit,nft,ESH,SelfEneR,str%cblk,tun) 
    endif

    TUN_MAT(1) = tun 

    ! When more contacts are present sometimes we can re-use previous Gr 
    do icpl = 2, size_ni

      nit=ni(icpl)
      nft=nf(icpl)

      if (str%cblk(nit).gt.str%cblk(nft)) then
        nt1 = str%cblk(nit)
      else
        nt1 = str%cblk(nft)
      endif

      ! if nt1 > nt extend the Gr calculation
      if (nt1 .gt. nt) then
        call Make_Gr_mem_dns(ESH,nt+1,nt1)
        nt = nt1
      endif

      call trasmission_old(nit,nft,ESH,SelfEneR,str%cblk,tun) 

      TUN_MAT(icpl) = tun 

    enddo

    !Distruzione delle Green
    do i=2, nt 
      call destroy(Gr(i,i))
      call destroy(Gr(i-1,i))
      call destroy(Gr(i,i-1))
    enddo
    call destroy(Gr(1,1))

    call deallocate_blk_dns(Gr)

    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)

    call destroy_ESH(ESH)
    call deallocate_blk_dns(ESH)

  end SUBROUTINE transmission_dns

  !------------------------------------------------------------------------------------------------!

  !************************************************************************
  !
  ! Subroutine for transmission calculation
  !
  !************************************************************************
  ! NOTE:
  !
  !  This subroutine was hacked quickly to obain effiecient transmission calcs
  !  Useful only when there are 2 contacts
  !                ===================
  !************************************************************************

  subroutine trasmission_dns(ni,nf,ESH,SelfEneR,cblk,TUN)

    implicit none

    !In/Out
    Integer :: ni,nf
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    Type(z_DNS), Dimension(:,:) :: ESH
    Integer, Dimension(:), pointer :: cblk
    Real(dp) :: TUN

    !Work variables
    Integer :: ct1, bl1
    Type(z_DNS) :: work1, work2, GAM1_dns, GA, TRS, AA
    Complex(dp), PARAMETER ::    j = (0.d0,1.d0)  ! CMPX unity

    if (size(cblk).gt.2) then
      write(*,*) "ERROR: transmission_dns is valid only for 2 contacts"
      return
    endif

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
      ct1=ni;
    else
      ct1=nf;
    endif

    bl1=cblk(ct1); 

    call zdagger(Gr(bl1,bl1),GA)

    ! Computes the Gamma matrices
    call zspectral(SelfEneR(ct1),SelfEneR(ct1),0,GAM1_dns)

    ! Work to compute transmission matrix (Gamma G Gamma G)
    call prealloc_mult(GAM1_dns,Gr(bl1,bl1),work1)

    call prealloc_mult(work1,GAM1_dns,work2)

    call destroy(work1)

    call prealloc_mult(work2,GA,work1)

    call destroy(work2)

    call create(AA,GA%nrow,GA%ncol)

    AA%val = j * (Gr(bl1,bl1)%val-GA%val)

    call destroy(GA) 

    call prealloc_mult(GAM1_dns,AA,work2) 

    call destroy(GAM1_dns,AA)

    call create(TRS,work1%nrow,work1%ncol)

    TRS%val = work2%val - work1%val

    TUN = abs( real(trace(TRS)) )  

    call destroy(TRS,work1,work2)

  end subroutine trasmission_dns

  !------------------------------------------------------------------------------------------------!
  
  !************************************************************************
  !
  ! Subroutine for transmission calculation (GENERIC FOR N CONTACTS)
  !
  !************************************************************************ 
  subroutine trasmission_old(ni,nf,ESH,SelfEneR,cblk,TUN)

    !In/Out
    Integer :: ni,nf
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    Type(z_DNS), Dimension(:,:) :: ESH
    Integer, Dimension(:), pointer :: cblk
    Real(kind=dp) :: TUN

    !Work variables
    Integer :: ct1, ct2, bl1, bl2, i, nbl
    Type(z_DNS) :: work1, work2, GAM1_dns, GAM2_dns, GA, TRS
    Real(kind=dp) :: max

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
      ct1=ni;ct2=nf;
    else
      ct1=nf;ct2=ni;
    endif

    bl1=cblk(ct1); bl2=cblk(ct2);
    nbl = size(cblk)
    ! in this way nt1 < nt2 by construction

    if ( nbl.gt.1 .and. (bl2-bl1).gt.1) then

      ! Compute column-blocks of Gr(i,bl1) up to i=bl2
      ! Gr(i,bl1) = -gr(i) T(i,i-1) Gr(i-1,bl1)
      do i = bl1+1, bl2
        !Checks whether previous block is non null. 
        !If so next block is also null => TUN = 0       
        max=maxval(abs(Gr(i-1,bl1)%val))

        if (max.lt.EPS) then
          TUN = EPS*EPS !for log plots 
          !Destroy also the block adjecent to diagonal since 
          !this is not deallocated anymore in calling subroutine
          if (i.gt.(bl1+1)) call destroy(Gr(i-1,bl1))
          return
        endif

        !Checks whether block has been created, if not do it 
        if (.not.allocated(Gr(i,bl1)%val)) then 

          call prealloc_mult(gsmr(i),ESH(i,i-1),(-1.d0, 0.d0),work1)

          call prealloc_mult(work1,Gr(i-1,bl1),Gr(i,bl1))

          call destroy(work1)

        endif

        ! avoid destroying blocks closer to diagonal
        if (i.gt.(bl1+2)) call destroy(Gr(i-1,bl1))

      enddo

    endif

    ! Computes the Gamma matrices
    call zspectral(SelfEneR(ct1),SelfEneR(ct1),0,GAM1_dns)
    call zspectral(SelfEneR(ct2),SelfEneR(ct2),0,GAM2_dns)

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call prealloc_mult(GAM2_dns,Gr(bl2,bl1),work1)

    call destroy(GAM2_dns)

    call prealloc_mult(work1,GAM1_dns,work2)

    call destroy(work1)

    call destroy(GAM1_dns)

    call zdagger(Gr(bl2,bl1),GA)

    if (bl2.gt.bl1+1) call destroy( Gr(bl2,bl1) )

    call prealloc_mult(work2,GA,TRS)

    call destroy(work2)

    call destroy(GA) 

    TUN = real(trace(TRS))

    call destroy(TRS)

  end subroutine trasmission_old

  !------------------------------------------------------------------------------------------------!
  !Subroutine for transmission and dos calculation    
  !------------------------------------------------------------------------------------------------!  

  subroutine tun_and_dos(H,S,Ec,SelfEneR,Gs,ni,nf,nLDOS,LDOS,size_ni,str,TUN_MAT,LEDOS)

    implicit none

    Type(z_CSR), intent(in) :: H
    Type(z_CSR), intent(in) :: S           
    Complex(dp), intent(in) :: Ec
    Type(z_DNS), Dimension(MAXNCONT), intent(in) :: SelfEneR, Gs
    Integer, intent(in) :: ni(MAXNCONT)
    Integer, intent(in) :: nf(MAXNCONT)
    Type(TStruct_Info), intent(in) :: str
    integer, intent(in)  :: nLdos, size_ni
    type(intarray), dimension(:), intent(in) :: LDOS      
    Real(dp), Dimension(:), intent(inout) :: TUN_MAT
    Real(dp), Dimension(:), intent(inout) :: LEDOS

    ! Local variables
    Type(z_CSR) :: ESH_tot, GrCSR
    Type(z_DNS), Dimension(:,:), allocatable :: ESH
    Type(r_CSR) :: Grm                          ! Green Retarded nella molecola           
    real(dp), dimension(:), allocatable :: diag
    Real(dp) :: tun
    Complex(dp) :: zc
    Integer :: nbl,ncont, ierr
    Integer :: nit, nft, icpl
    Integer :: iLDOS, i2, i
    Character(1) :: Im
   

    nbl = str%num_PLs
    ncont = str%num_conts
    Im = 'I'

    ! Create ESH : Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(H,S,(-1.d0, 0.d0),Ec,ESH_tot)    
    call allocate_blk_dns(ESH,nbl)
    call zcsr2blk_sod(ESH_tot,ESH,str%mat_PL_start)
    call destroy(ESH_tot)

    ! Add the retarded contact self-energies to the relevant blocks
    do i=1,ncont
      ESH(str%cblk(i),str%cblk(i))%val = ESH(str%cblk(i),str%cblk(i))%val-SelfEneR(i)%val
    enddo

    !--------------------------------------------------------------------------!
    ! Retarded Green function calculation
    !--------------------------------------------------------------------------!

    call allocate_gsm_dns(gsmr,nbl)
    call Make_gsmr_mem_dns(ESH,nbl,2) 
    call allocate_blk_dns(Gr,nbl)
    call Make_Gr_mem_dns(ESH,1)
    call Make_Gr_mem_dns(ESH,2,nbl)

    !Computation of transmission(s) between contacts ni(:) -> nf(:)
    do icpl=1,size_ni

      nit=ni(icpl)
      nft=nf(icpl)

      if (ncont.eq.2) then
        call trasmission_dns(nit,nft,ESH,SelfEneR,str%cblk,tun) 
      else
        call trasmission_old(nit,nft,ESH,SelfEneR,str%cblk,tun)
      endif

      TUN_MAT(icpl) = tun 

    enddo
    
    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)

    !Distruzione dei blocchi fuori-diagonale
    do i=2,nbl
      call destroy(Gr(i-1,i))
      call destroy(Gr(i,i-1))
    enddo
    call destroy_ESH(ESH)
    call deallocate_blk_dns(ESH)

    call create(Grm,str%mat_PL_start(nbl+1)-1,str%mat_PL_start(nbl+1)-1,0)

    Grm%rowpnt(:)=1

    do i=1,nbl
      call create(GrCSR,Gr(i,i)%nrow,Gr(i,i)%ncol,Gr(i,i)%nrow*Gr(i,i)%ncol)
      call dns2csr(Gr(i,i),GrCSR)
      !Concatena direttamente la parte immaginaria per il calcolo della DOS
      zc=(-1.d0,0.d0)/pi
      
      call concat(Grm,zc,GrCSR,Im,str%mat_PL_start(i),str%mat_PL_start(i))
      call destroy(Gr(i,i))
      call destroy(GrCSR)
    enddo

    call deallocate_blk_dns(Gr)

    !Compute LDOS on the specified intervals
    if (nLDOS.gt.0) then
      call log_allocate(diag, Grm%nrow)
      call getdiag(Grm,diag)
      do iLDOS=1,nLDOS
        do i = 1, size(LDOS(iLDOS)%indexes)
          i2 = LDOS(iLDOS)%indexes(i)
          if (i2 .le. str%central_dim) then
            LEDOS(iLDOS) = LEDOS(iLDOS) + diag(i2)  
          end if
        end do
      enddo
      call log_deallocate(diag)
    endif

    call destroy(Grm)

  end subroutine tun_and_dos  

  !------------------------------------------------------------------------------------------------!

  !> Calculation of the retarded Green function (extended diagonal) at a single energy point.
  !> Note that the function A is not full GF.
  !> (public)
  subroutine iterativeGreenRetarded(pnegf,E,SelfEneR,Tlc,Tcl,gsurfR,A,outer)

    !----------------------------------------------------------------------------------------------!
    !> Input
    !! pnegf:    negf data container
    !! E:        energy point
    !! SelfEneR: matrices array containing contacts self-energy
    !! Tlc:      matrices array containing contacts-device interaction blocks (ES-H)
    !! Tcl:      matrices array containing device-contacts interaction blocks (ES-H)
    !! gsurfR:   matrices array containing contacts surface green
    !! outer:    parameter (0,1,2)
    !!
    !> Output
    !! A: Retarded Green function (Device + Contacts overlap regions -> effective conductor)
    !!    outer = 0  no outer parts are computed
    !!    outer = 1  only D/C part is computed
    !!    outer = 2  D/C and C/D parts are computed
    !!               (needed for K-points calculations)
    !----------------------------------------------------------------------------------------------!

    implicit none

    !In/Out
    TYPE(Tnegf), intent(inout) :: pnegf
    COMPLEX(dp), intent(in) :: E
    TYPE(z_DNS), DIMENSION(:), intent(in) :: SelfEneR
    TYPE(z_DNS), DIMENSION(:), intent(in) :: Tlc, Tcl, gsurfR
    TYPE(z_CSR), intent(out) :: A
    INTEGER, intent(in) :: outer

    !Work
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_CSR) :: ESH_tot, Ain
    INTEGER :: i,ierr, nbl, ncont,ii,n                                
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk

    if (id0.and.pnegf%verbose.gt.130) print *, 'debug: iterativeGreenRetarded is started'

    nbl = pnegf%str%num_PLs
    ncont = pnegf%str%num_conts
    cblk => pnegf%str%cblk
    indblk => pnegf%str%mat_PL_start

    ! Create ESH : Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(pnegf%H,pnegf%S,(-1.d0, 0.d0),E,ESH_tot)
    call allocate_blk_dns(ESH,nbl)
    call zcsr2blk_sod(ESH_tot,ESH,indblk)
    call destroy(ESH_tot)

    ! Add the retarded contact self-energies to the relevant blocks
    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    end do

    ! Add the retarded dephasing and/or many-body self-energies 
    call add_sigma_r(pnegf,ESH) 

    !--------------------------------------------------------------------------!
    ! Retarded Green function calculation
    !--------------------------------------------------------------------------!
    
    call allocate_gsm_dns(gsmr,nbl)      ! allocates type z_DNS(:)
    call Make_gsmr_mem_dns(ESH,nbl,2)    ! g_small right (gsmr) calculation
    call allocate_blk_dns(Gr,nbl)        ! allocates type z_DNS(:,:)
    call Make_Gr_mem_dns(ESH,1)          ! Diagonal block Gr(1,1) of Green Retarded 
    call Make_Gr_mem_dns(ESH,2,nbl)      ! Gr(i,i), Gr(i-1,i), Gr(i,i-1) for i=2,nbl
    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)

    ! SAVE ON FILES/MEMORY (for elph).........................
    !if (pnegf%elph%numselmodes.gt.0 .and. pnegf%elph%model .eq. -1) then
      ! save diagonal blocks of Gn = -i G<
    !  DO i = 1, nbl
    !    !print*,'G_r ',minval(abs(Gr(i,i)%val)), maxval(abs(Gr(i,i)%val))
    !    call write_blkmat(Gr(i,i),pnegf%scratch_path,'G_r_',i,i,pnegf%iE)
    !  ENDDO
    !endif
    
    ! Update vibronic dephasing retarded self-energy if any
    if (allocated(pnegf%inter)) call pnegf%inter%set_Gr(Gr, pnegf%iE, pnegf%inter%Mixing)

    call blk2csr(Gr,pnegf%str,pnegf%S,A)

    SELECT CASE (outer)
    CASE(0)
    CASE(1)
      CALL Outer_Gr_mem_dns(Tlc,Tcl,gsurfR,pnegf%str,.FALSE.,A)   
    CASE(2)
      CALL Outer_Gr_mem_dns(Tlc,Tcl,gsurfR,pnegf%str,.TRUE.,A) 
    END SELECT
    
    ! Destroy
    call destroy_blk(Gr)
    deallocate(Gr)
    call destroy_ESH(ESH)
    call deallocate_blk_dns(ESH)

    if (id0.and.pnegf%verbose.gt.130) print *, 'debug: iterativeGreenRetarded is finished'

  end subroutine iterativeGreenRetarded
  
  !------------------------------------------------------------------------------------------------!


  !****************************************************************************
  !
  ! Driver for computing G_n contributions due to all contacts but reference:
  !
  !   Sum   [f_j(E)-f_r(E)] Gr Gam_j Ga
  !   j!=r
  !
  ! NOTE: The subroutine assumes that 
  !
  !****************************************************************************

  !> Calculation of the lesser (-iG<) Green function (extended diagonal) at a single energy point.
  !> ONLY for NONINTERACTING (Landauer, DFTB+NEGF)
  !> (public)
  subroutine iterativeGreenLesser_Landauer(pnegf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,Glout,out)

    !****************************************************************************
    !
    !Input
    !pnegf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !Tlc:      matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl:      matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR:   matrices array containing contacts surface green
    !out:      optional parameter (0,1,2).  
    !
    !Output:
    !Gn: NE GF (Device + Contacts overlap regions -> effective conductor)
    !   out = 0  no outer parts are computed
    !   out = 1  only D/C part is computed
    !   out = 2  D/C and C/D parts are computed
    !              (needed for K-points calculations)
    !
    !*****************************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(Tnegf), intent(inout) :: pnegf
    TYPE(z_CSR), intent(inout)  :: Glout
    TYPE(z_DNS), DIMENSION(:), intent(in)  :: SelfEneR, gsurfR, Tlc, Tcl
    REAL(dp), intent(in)  :: E
    REAL(dp), DIMENSION(:), intent(in)  :: frm
    INTEGER, intent(in)  :: out

    !Work
    INTEGER :: ref
    COMPLEX(dp) :: Ec
    INTEGER :: i,ierr,ncont,nbl, lbl, rbl
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Gn
    TYPE(z_CSR) :: ESH_tot, Gl
    LOGICAL :: mask(MAXNCONT)

    if (id0.and.pnegf%verbose.gt.130) print *, 'debug: iterativeGreenLesser_Landauer is started'

    if (pnegf%elph%model .ne. 0) then
      !TODO: output logfile
      write(*,*) "Warning:  iterativeGreenLesser_Landauer is not compatible with el-ph models"; stop
    endif

    nbl = pnegf%str%num_PLs
    ncont = pnegf%str%num_conts
    indblk => pnegf%str%mat_PL_start
    cblk => pnegf%str%cblk
    ref = pnegf%refcont

    Ec=cmplx(E,0.d0,dp)

    ! Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(pnegf%H,pnegf%S,(-1.d0, 0.d0),Ec,ESH_tot)
    call allocate_blk_dns(ESH,nbl)
    call zcsr2blk_sod(ESH_tot,ESH,indblk)
    call destroy(ESH_tot)

    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    end do

    call allocate_gsm_dns(gsmr,nbl)
    call allocate_gsm_dns(gsml,nbl)

    ! Compute blocks for gsmr and gsml
    mask = .true.
    mask(ref) = .false. 
    rbl = minval(cblk(1:ncont),mask(1:ncont)) + 1  
    lbl = maxval(cblk(1:ncont),mask(1:ncont)) - 1

    ! Fix to a bug when there are 2PLs
    ! later Make_Gr tries to compute Gr(1,1) but needs gsmr(2,2)
    ! 
    IF (nbl.eq.2) then
      CALL Make_gsmr_mem_dns(ESH,nbl,rbl-1)
      CALL Make_gsml_mem_dns(ESH,1,lbl+1)    
    ELSE
      CALL Make_gsmr_mem_dns(ESH,nbl,rbl)
      CALL Make_gsml_mem_dns(ESH,1,lbl)    
    ENDIF

    call allocate_blk_dns(Gr,nbl)

    ! -------------------------------------------------------------
    ! 1. rbl>lbl  => lbl+1=rbl-1 => compute first Gr(rbl-1,rbl-1)
    ! 2. rbl<lbl  => lbl=rbl-2 has been computed
    ! Make_Gr does not compute if sbl>nbl or sbl<1
    CALL Make_Gr_mem_dns(ESH,rbl-1)
    CALL Make_Gr_mem_dns(ESH,rbl,nbl)
    CALL Make_Gr_mem_dns(ESH,rbl-2,1)

    !Chiamata di Make_Grcol_mem per i contatti necessari 
    DO i=1,ncont
      IF (i.NE.ref) THEN
        CALL Make_Grcol_mem_dns(ESH,cblk(i),indblk)
      ENDIF
    ENDDO

    !Distruzione delle gsmall................................
    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)
    call destroy_gsm(gsml)
    call deallocate_gsm_dns(gsml)
    !.......................................................

    !Calcolo della G_n nel device
    call allocate_blk_dns(Gn,nbl)

    call init_blkmat(Gn,ESH)

    CALL Make_Gn_mem_dns(ESH,SelfEneR,frm,ref,pnegf%str,Gn)

    call blk2csr(Gn,pnegf%str,pnegf%S,Glout)

    !Calcolo degli outer blocks 
    SELECT CASE (out)
    CASE(0)
    CASE(1)
      CALL Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.false.,Glout)
    CASE(2)
      CALL Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.true.,Glout)
    END SELECT

    CALL destroy_blk(Gn)
    DEALLOCATE(Gn)

    CALL destroy_blk(Gr)
    DEALLOCATE(Gr)

    CALL destroy_ESH(ESH)
    DEALLOCATE(ESH)

    if (id0.and.pnegf%verbose.gt.130) print *, 'debug: iterativeGreenLesser_Landauer is finished'

  end subroutine iterativeGreenLesser_Landauer

  !****************************************************************************
  !
  ! Driver for computing G_n = -iG< contributions including el-ph interactions
  !
  !   Sum   f_j(E) Gr Gam_j Ga +   Gr Sigma_ph< Ga
  !    j
  !
  ! NOTE: The subroutine assumes that
  !
  !****************************************************************************

  !> Calculation of the lesser (-iG<) Green function (extended diagonal) at a single energy point.
  !> (public)
  subroutine iterativeGreenLesser(tranas,E,SelfEneR,Tlc,Tcl,gsurfR,frm,Glout,out)

    !****************************************************************************
    !
    !Input
    !H: sparse matrix contaning Device Hamiltonian
    !S: sparse matrix containing Device Overlap
    !SelfEneR: sparse matrices array containing contacts Self Energy
    !Tlc: sparse matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl: sparse matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR: sparse matrices array containing contacts surface green
    !frm: array containing Fermi distribution values for all contacts
    !ref: reference contact excluded from summation
    !
    !Output:
    !Aout: G_n contributions (Device + Contacts overlap regions -> effective conductor)
    !
    !*****************************************************************************

    implicit none

    !In/Out
    type(TTraNaS), target :: tranas
    type(Tnegf), pointer :: pnegf
    TYPE(z_CSR) :: Glout
    TYPE(z_DNS), DIMENSION(:) :: SelfEneR, gsurfR, Tlc, Tcl
    REAL(dp) :: E
    REAL(dp), DIMENSION(:) :: frm
    INTEGER :: out

    !Work
    COMPLEX(dp) :: Ec
    INTEGER :: i,ierr,ncont,nbl,lbl,ii,n,m                              
    INTEGER :: ref, iter
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH, Gn, Gp
    TYPE(z_CSR) :: ESH_tot, Gl
    LOGICAL :: mask(MAXNCONT)
    REAL(dp), DIMENSION(:), allocatable :: cfrm

    pnegf => tranas%negf
    
!    type(z_dns), dimension(:,:), allocatable :: Gr_previous !, Gn_previous !DAR
!    save Gr_previous !, Gn_previous !DAR

    if (id0.and.pnegf%verbose.gt.130) print *, 'debug: iterativeGreenLesser is started'

    nbl = pnegf%str%num_PLs
    ncont = pnegf%str%num_conts
    indblk => pnegf%str%mat_PL_start
    cblk => pnegf%str%cblk
    ref = pnegf%refcont
    iter=0
    if(allocated(pnegf%inter)) iter = pnegf%inter%scba_iter

    Ec=cmplx(E,0.d0,dp)

    ! Create ESH : Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(pnegf%H,pnegf%S,(-1.d0, 0.d0),Ec,ESH_tot)
    call allocate_blk_dns(ESH,nbl)
    call zcsr2blk_sod(ESH_tot,ESH,indblk)
    call destroy(ESH_tot)

    ! Add the retarded contact self-energies
    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    end do

    ! Add the retarded dephasing and/or many-body self-energies
    call add_sigma_r(pnegf,ESH) 

    !--------------------------------------------------------------------------!
    ! Retarded Green function calculation
    !--------------------------------------------------------------------------!
    
    call allocate_gsm_dns(gsmr,nbl)
    call allocate_gsm_dns(gsml,nbl)
    call Make_gsmr_mem_dns(ESH,nbl,2)
    !Chiamata di Make_gsml_mem solo per i blocchi 1..lbl dove  
    ! lbl = maxval(cblk,mask) - 2 
    lbl = nbl - 2    ! ALL BLOCKS are needed!!
    if(ncont.gt.1) call Make_gsml_mem_dns(ESH,1,lbl)    
    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'
    call allocate_blk_dns(Gr,nbl)
    call Make_Gr_mem_dns(ESH,1)
    call Make_Gr_mem_dns(ESH,2,nbl)

    ! Update vibronic dephasing retarded self energy if any
    if (allocated(pnegf%inter)) call pnegf%inter%set_Gr(Gr, pnegf%iE, pnegf%inter%Mixing)

    ! With el-ph we need all columns (?)
    do i=1,nbl
      call Make_Grcol_mem_dns(ESH,i,indblk)
    enddo

    ! Destroy gsmr, gsml
    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)
    call destroy_gsm(gsml)
    call deallocate_gsm_dns(gsml)
    
    !--------------------------------------------------------------------------!
    ! Lesser Green function calculation
    !--------------------------------------------------------------------------!
    
    ! Outer blocks
    select case (out)
      case(0)
      case(1)
        call Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.false.,Glout)
      case(2)
        call Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.true.,Glout)
    end select
    allocate(Gn(nbl,nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate Gn(nbl,nbl)'
    call init_blkmat(Gn,ESH)
    call Make_Gn_mem_dns(ESH,SelfEneR,frm,ref,pnegf%str,Gn)
    if(allocated(pnegf%inter)) call Make_Gn_ph(pnegf,ESH,iter,Gn)
    !if(pnegf%tMBNGF) call Make_Gn_mbngf(pnegf,ESH,iter,Gn)
    
    ! Update vibronic dephasing lesser self energy if any
    if (allocated(pnegf%inter)) call pnegf%inter%set_Gn(Gn, pnegf%iE, pnegf%inter%Mixing)
        
    ! The following destroys Gn
    call blk2csr(Gn,pnegf%str,pnegf%S,Gl)
    deallocate(Gn)

    ! Destroy
    call destroy_blk(Gr)
    deallocate(Gr)
    call destroy_ESH(ESH)
    deallocate(ESH)

    ! Concatenazione di Gl in Glout
    SELECT CASE(out)
    CASE(0)
      call clone(Gl,Glout)
    CASE(1:2)
      call concat(Glout,Gl,1,1)
    END SELECT

    call destroy(Gl)

    if (id0.and.pnegf%verbose.gt.130) print *, 'debug: iterativeGreenLesser is finished'

  end SUBROUTINE iterativeGreenLesser

  !---------------------------------------------------------------------
  !>
  !  Iterative algorithm implementing Meir Wingreen formula for a given
  !  electrode
  !  Note: self consistent born approximation is not accounted for here
  !  It is assumed that the el-ph container already includes the 
  !  desired values. SCBA loop should be run outside
  !  Many operations from calls_neq_ph are repeated here, as it is
  !  assumed that A and Gn are not available at the time of the call
  !
  !  It implements the form without the reference electrode
  !  Iop = \Sigma_{i}^{n}A - \Gamma_{i}G^{n} =
  !      = \Gamma_{i}[f_{i}A - f_{ref}A - G^{n,l\neq ref} +
  !         - G^{n,l\neq ref}_{\phi}]
  !  where G^{n,l\neq ref} is the component including no el-ph
  !  If i=ref it reduces to
  !  Iop = \Gamma_{i}[-G^{n,l\neq ref}  - G^{n,l\neq ref}_{\phi}]
  !---------------------------------------------------------------------
  subroutine iterativeMeirWingreen(pnegf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,ni,tun_mat)
    IMPLICIT NONE

    !In/Out
    TYPE(Tnegf), intent(inout) :: pnegf
    integer, dimension(:) :: ni
    TYPE(z_DNS), DIMENSION(:) :: SelfEneR, gsurfR, Tlc, Tcl
    Real(dp), Dimension(:), allocatable :: tun_mat
    REAL(dp) :: E
    REAL(dp), DIMENSION(:) :: frm
    INTEGER :: out, size_ni
    complex(dp) :: tmp

    !Work
    COMPLEX(dp) :: Ec
    INTEGER :: i,ierr,ncont,nbl,lbl,ii,n                              !DAR +ii,n
    INTEGER :: ref, iter, lead, lead_blk, ref_blk
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH, Gn, Gp
    type(z_DNS) :: work1, work2, work3, Gamma, A
    TYPE(z_CSR) :: ESH_tot, Gl
    LOGICAL :: mask(MAXNCONT)
    REAL(dp), DIMENSION(:), allocatable :: cfrm

    if (id0.and.pnegf%verbose.gt.90) print *, 'iterativeMeirWingreen is started'
    
    nbl = pnegf%str%num_PLs
    ncont = pnegf%str%num_conts
    indblk => pnegf%str%mat_PL_start
    cblk => pnegf%str%cblk
    ref = pnegf%refcont
    ref_blk = pnegf%str%cblk(ref)
    iter = 0
    if(allocated(pnegf%inter)) iter = pnegf%inter%scba_iter

    !print *, 'debug: nbl=',nbl,' ncont=',ncont
    !print *, 'debug: indblk=',indblk
    !print *, 'debug: cblk=',cblk
    !print *, 'debug: ref=',ref,' ref_blk=',ref_blk,' iter=',iter

    Ec=cmplx(E,0.d0,dp)

    ! Create ESH : Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(pnegf%H,pnegf%S,(-1.d0, 0.d0),Ec,ESH_tot)
    call allocate_blk_dns(ESH,nbl)
    call zcsr2blk_sod(ESH_tot,ESH,indblk)
    call destroy(ESH_tot)

    ! Add the retarded contact self-energies
    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    end do

    !--------------------------------------------------------------------------!
    ! Retarded Green function calculation
    !--------------------------------------------------------------------------!
    
    ! Add the retarded dephasing and /or many-body self-energies
    call add_sigma_r(pnegf,ESH) !DAR
    
    ! Allocation of gsmr and gsml
    call allocate_gsm_dns(gsmr,nbl)
    call allocate_gsm_dns(gsml,nbl)

    ! Make gsmr
    call Make_gsmr_mem_dns(ESH,nbl,2)

    !Chiamata di Make_gsml_mem solo per i blocchi 1..lbl dove  
    ! lbl = maxval(cblk,mask) - 2 
    lbl = nbl - 2    ! ALL BLOCKS are needed!!

    if( ncont.gt.1 ) then
      CALL Make_gsml_mem_dns(ESH,1,lbl)    
    endif
    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    call allocate_blk_dns(Gr,nbl)

    CALL Make_Gr_mem_dns(ESH,1)
    CALL Make_Gr_mem_dns(ESH,2,nbl)

    !! Give Gr to interaction model if any
    if (allocated(pnegf%inter)) call pnegf%inter%set_Gr(Gr, pnegf%iE, pnegf%inter%Mixing)
    !---------------------------------------------------
    !With el-ph we need all columns
    DO i=1,nbl
      CALL Make_Grcol_mem_dns(ESH,i,indblk)
    ENDDO

    !Distruzione delle gsmall
    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)
    call destroy_gsm(gsml)
    call deallocate_gsm_dns(gsml)

    !--------------------------------------------------------------------------!
    ! Lesser Green function calculation
    !--------------------------------------------------------------------------!

    !! Never calculate outer blocks

    !Calcolo della G_n nel device
    !if (debug) write(*,*) 'Compute G_n' 
    !Allocazione degli array di sparse
    ALLOCATE(Gn(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gn(nbl,nbl)'

    call init_blkmat(Gn,ESH)

    !! TEMPORARY AND INEFFICIENT: CALCULATE THE FULL GN WHEN A BLOCK 
    !! (THE CONTACT ONE) IS ENOUGH
    !! WE HAVE Gn WITHOUT REFERENCE CONTACT?? (usual neq contributions)
    CALL Make_Gn_mem_dns(ESH,SelfEneR,frm,ref,pnegf%str,Gn)

    if(allocated(pnegf%inter)) call Make_Gn_ph(pnegf,ESH,iter,Gn) !DAR
    !if(pnegf%tMBNGF) call Make_Gn_mbngf(pnegf,ESH,iter,Gn) !DAR

    ! I probably need the next set_Gn in interaction for current conservation here

    !--------------------------------------------------------------------------!
    ! Current spectral density calculation
    !--------------------------------------------------------------------------!

    do i=1,size(ni)       
      if (ni(i) .eq. 0) then
        cycle
      end if    
      lead = ni(i)
      lead_blk = pnegf%str%cblk(lead)
      !print *, 'debug: lead=',lead,' lead_blk=',lead_blk
      call zspectral(SelfEneR(lead),SelfEneR(lead), 0, Gamma)
      if (lead.eq.ref) then
        call prealloc_mult(Gamma, Gn(lead_blk, lead_blk), (-1.d0, 0.d0), work1)
        tun_mat(i) = - real(trace(work1))                      !DAR "-" is added
        call destroy(work1)
      else       
        !call zspectral(Gr(ref, ref), Gr(ref, ref), 0, A)                   !DAR
        call zspectral(Gr(lead_blk, lead_blk), Gr(lead_blk, lead_blk), 0, A)!DAR
        tmp = frm(lead)-frm(ref)
        call prealloc_sum(A, Gn(lead_blk, lead_blk), tmp, (-1.d0, 0.d0), work1)
        call destroy(A)
        call prealloc_mult(Gamma, work1, work2)
        call destroy(work1)
        tun_mat(i) = - real(trace(work2))                      !DAR "-" is added
        call destroy(work2)
      endif
      call destroy(Gamma)
    enddo

    DEALLOCATE(Gn)

    !--------------------------------------------------------------------------!
    !--------------------------------------------------------------------------!

    if(pnegf%tZeroCurrent) call transmission_BP_corrected(pnegf,SelfEneR)
    
    !Distruzione dell'array Gr
    CALL destroy_blk(Gr)
    DEALLOCATE(Gr)
 
    CALL destroy_ESH(ESH)
    DEALLOCATE(ESH)

    if (id0.and.pnegf%verbose.gt.90) print *, 'iterativeMeirWingreen is finished'

  end subroutine iterativeMeirWingreen

  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  
  subroutine add_sigma_r(pnegf,ESH)

    TYPE(Tnegf), intent(inout) :: pnegf
    type(z_DNS), dimension(:,:), allocatable, intent(inout) :: ESH
    integer :: n,ii,nbl

    if(allocated(pnegf%inter)) call pnegf%inter%add_sigma_r(ESH)

    if(pnegf%tDephasingBP) then
       nbl = pnegf%str%num_PLs
       do n=1,nbl
          associate(pl_start=>pnegf%str%mat_PL_start(n),pl_end=>pnegf%str%mat_PL_end(n))
          forall(ii = 1:pl_end - pl_start + 1) 
             ESH(n,n)%val(ii,ii) = ESH(n,n)%val(ii,ii)+ &
                  (0.d0,1.d0)*0.5_dp*pnegf%deph%bp%coupling(pl_start + ii - 1)  
          end forall
          end associate
       end do
    end if
  
    if(pnegf%tMBNGF) call pnegf%mbngf%add_SelfEnergyR(ESH)
       
  end subroutine add_sigma_r

  !----------------------------------------------------------------------------!
  
  subroutine transmission_BP_corrected(negf,SelfEner)

    TYPE(Tnegf), intent(inout) :: negf
    TYPE(z_DNS), DIMENSION(:) :: SelfEneR  
    integer :: nn,mm,k,l,m,n,ii,jj
    integer :: NumLevels
    real(dp) :: Transmission
    real(dp), allocatable  :: Trans(:,:),Gam1(:,:),Gam2(:,:),W(:,:),R(:)
    complex(dp), allocatable :: GreenR(:,:),SigmaR(:,:,:)

    integer :: INFO                                       ! LAPACK
    integer, allocatable :: IPIV(:)                       ! LAPACK
    complex, allocatable :: GG(:,:),WORK(:)               ! LAPACK

    if (id0.and.negf%verbose.gt.50) print *, 'transmission_BP_corrected is started'

    NumLevels=negf%str%central_dim

    allocate(Trans(0:NumLevels+1,0:NumLevels+1))
    allocate(Gam1(NumLevels,NumLevels),Gam2(NumLevels,NumLevels))
    allocate(W(NumLevels,NumLevels),R(0:NumLevels+1))
    Trans=0.
    allocate(GG(NumLevels,NumLevels),WORK(NumLevels),IPIV(NumLevels))
    allocate(GreenR(NumLevels,NumLevels),SigmaR(2,NumLevels,NumLevels))
    SigmaR=(0.0_dp,0.0_dp)

    do n=1,negf%str%num_PLs
       do m=1,negf%str%num_PLs   
          associate(pl_start1=>negf%str%mat_PL_start(n),pl_end1=>negf%str%mat_PL_end(n), &
               pl_start2=>negf%str%mat_PL_start(m),pl_end2=>negf%str%mat_PL_end(m))
          do ii = 1,pl_end1 - pl_start1 + 1
             do jj = 1,pl_end2 - pl_start2 + 1
                !print *,n,m,ii,jj,pl_start1,pl_end1,pl_start2,pl_end2,Gr(n,m)%val(ii,jj)
                GreenR(pl_start1 + ii - 1,pl_start2 + jj - 1) = Gr(n,m)%val(ii,jj)
                !print *,n,m,ii,jj,pl_start1,pl_end1,pl_start2,pl_end2,GreenR(pl_start1 + ii - 1,pl_start2 + jj - 1)
             end do
          end do
          end associate
       end do
    end do

    !print *, GreenR

    do n=1,2
       associate(pl_start=>negf%str%mat_PL_start(negf%str%cblk(n)),pl_end=>negf%str%mat_PL_end(negf%str%cblk(n)))
       do ii = 1,pl_end - pl_start + 1
          do jj = 1,pl_end - pl_start + 1
             !print *, n, negf%str%cblk(n), pl_start, pl_end
             !print *,n,ii,jj,pl_start,pl_end,SelfEner(n)%val(ii,jj)
             SigmaR(n,pl_start + ii - 1,pl_start + jj - 1) = SelfEner(n)%val(ii,jj)
             !print *,n,ii,jj,pl_start,pl_end,SigmaR(n,pl_start + ii - 1,pl_start + jj - 1)
          end do
       end do
       end associate
    end do

    !print *, SigmaR(1,:,:)
    !print *, SigmaR(2,:,:)
    
    do nn=0,NumLevels+1
       do mm=0,NumLevels+1     

          Gam1=0.
          if(nn.eq.0) then
             Gam1(:,:)=-2*aimag(SigmaR(1,:,:))
          else if(nn.eq.NumLevels+1) then
             Gam1(:,:)=-2*aimag(SigmaR(2,:,:))      
          else
             Gam1(nn,nn)=negf%deph%bp%coupling(nn)
          end if
          Gam2=0.
          if(mm.eq.0) then
             Gam2(:,:)=-2*aimag(SigmaR(1,:,:))
          else if(mm.eq.NumLevels+1) then
             Gam2(:,:)=-2*aimag(SigmaR(2,:,:))       
          else
             Gam2(mm,mm)=negf%deph%bp%coupling(mm)
          end if

          Trans(nn,mm)=0.
          do n=1,NumLevels
             do m=1,NumLevels
                do k=1,NumLevels
                   do l=1,NumLevels
                      Trans(nn,mm)=Trans(nn,mm)+Gam1(n,m)*GreenR(m,k)*Gam2(k,l)*conjg(GreenR(n,l))
                   end do
                end do
             end do
          end do

          !   Trans(nn,mm)=0.
          !   do n=1,NumLevels
          !      GG=matmul(Gam1,matmul(GreenR,matmul(Gam2,transpose(conjg(GreenR)))))
          !      Trans(nn,mm)=Trans(nn,mm)+real(GG(n,n))
          !   end do
          
          !write(*,"('Trans(',I0,',',I0,',',I0,')=',F12.6)")nn,mm,i,Trans(nn,mm,i)
          if (id0.and.negf%verbose.gt.80) write(*,"('    T(',I0,',',I0,') is calculated')")nn,mm
       end do
    end do

    do n=0,NumLevels+1
       R(n)=1.
       do m=0,NumLevels+1
          if(m.ne.n) R(n)=R(n)-Trans(m,n)
       end do
       !write(*,"('R(',I0,',',I0,')=',F12.6)")n,i,R(n,i)
       if (id0.and.negf%verbose.gt.80) write(*,"('    R(',I0,') is calculated')")n
    end do

    do n=1,NumLevels
       do m=1,NumLevels   
          W(n,m)=-Trans(n,m)   
       end do
       W(n,n)=W(n,n)+1-R(n)+Trans(n,n)
    end do
    if (id0.and.negf%verbose.gt.80) write(*,"('    W is calculated')")

    GG=W

    call CGETRF(NumLevels,NumLevels,GG,NumLevels,IPIV,INFO)
    call CGETRI(NumLevels,GG,NumLevels,IPIV,WORK,NumLevels,INFO)

    W=GG

    if (id0.and.negf%verbose.gt.80) write(*,"('    W is inverted')")

    Transmission=Trans(0,NumLevels+1)
    do n=1,NumLevels
       do m=1,NumLevels
          Transmission=Transmission+Trans(0,n)*W(n,m)*Trans(m,NumLevels+1)
       end do
    end do

    if (id0.and.negf%verbose.gt.80) print *, 'Transmission(',negf%iE,')=',Transmission

    negf%tunn_mat_bp(negf%iE,1)=Transmission

    deallocate(Trans,Gam1,Gam2,W,R,GG,WORK,IPIV,GreenR,SigmaR)

    if (id0.and.negf%verbose.gt.50) print *, 'transmission_BP_corrected is finished'

  end subroutine transmission_BP_corrected
 
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  
  !**********************************************************************
  SUBROUTINE init_blkmat(Matrix,S)
    type(z_DNS), dimension(:,:) :: Matrix,S

    integer :: nbl, j

    nbl = SIZE(Matrix,1)

    call create(Matrix(1,1),S(1,1)%nrow,S(1,1)%ncol)      
    Matrix(1,1)%val=(0.d0,0.d0)
    DO j=2,nbl-1
      call create(Matrix(j-1,j),S(j-1,j)%nrow,S(j-1,j)%ncol)      
      Matrix(j-1,j)%val=(0.d0,0.d0)
      call create(Matrix(j,j),S(j,j)%nrow,S(j,j)%ncol)      
      Matrix(j,j)%val=(0.d0,0.d0)
      call create(Matrix(j,j-1),S(j,j-1)%nrow,S(j,j-1)%ncol)      
      Matrix(j,j-1)%val=(0.d0,0.d0)
    ENDDO
    IF (nbl.gt.1) then
      call create(Matrix(nbl,nbl),S(nbl,nbl)%nrow,S(nbl,nbl)%ncol)      
      Matrix(nbl,nbl)%val=(0.d0,0.d0)
    ENDIF

  END SUBROUTINE init_blkmat


  !**********************************************************************
  SUBROUTINE destroy_gsm(gsm)
    type(z_DNS), DIMENSION(:) :: gsm
    integer :: i, i1, nbl

    nbl=size(gsm,1)

    do i=1,nbl
      if (allocated(gsm(i)%val)) call destroy(gsm(i))
    end do

  END SUBROUTINE destroy_gsm

  !**********************************************************************
  SUBROUTINE destroy_blk(M)
    type(z_DNS), DIMENSION(:,:) :: M
    integer :: i, i1, nbl

    nbl=size(M,1)

    DO i=1,nbl
      DO i1=1,nbl
        IF (ALLOCATED(M(i1,i)%val)) THEN
          !print*,'kill Gr',i1,i
          CALL destroy(M(i1,i))
        END IF
      END DO
    END DO

  END SUBROUTINE destroy_blk

  !**********************************************************************
  SUBROUTINE destroy_ESH(ESH)

    integer :: i, nbl
    type(z_DNS), dimension(:,:) :: ESH

    nbl=size(ESH,1)

    DO i=1,nbl
      CALL destroy(ESH(i,i))
    END DO
    DO i=2,nbl
      CALL destroy(ESH(i-1,i))
      CALL destroy(ESH(i,i-1))
    END DO

  END SUBROUTINE destroy_ESH

  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !----------------------------------------- Makes for Gr -----------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!

  !> g_small right (gsmr) calculation
  !> Iterative calculation of the retarded GF from start block (sbl) to end block (ebl<sbl).
  subroutine Make_gsmr_mem_dns(ESH,sbl,ebl)

    implicit none

    ! In
    type(z_DNS), dimension(:,:), intent(in) :: ESH                    ! dense matrices array ESH(nbl,nbl)
    integer, intent(in) :: sbl,ebl                                    ! block indexes           

    ! Out gsmr(:) is global module variable !!!
    
    ! Work
    type(z_DNS) :: work1, work2
    integer :: nrow
    integer :: i, nbl

    if (sbl.lt.ebl) return

    nbl = size(ESH,1)

    if (nbl.gt.1) then
      nrow=ESH(sbl,sbl)%nrow
      call create(gsmr(sbl),nrow,nrow)                                ! module mat_def
      call compGreen(gsmr(sbl),ESH(sbl,sbl),nrow)                     ! module inversions
    endif

    do i=sbl-1,ebl,-1
      call prealloc_mult(ESH(i,i+1),gsmr(i+1),(-1.d0, 0.d0),work1)    ! module sparsekit_drv  
      call prealloc_mult(work1,ESH(i+1,i),work2)
      call destroy(work1)
      call prealloc_sum(ESH(i,i),work2,work1)
      call destroy(work2)
      call create(gsmr(i),work1%nrow,work1%nrow)
      call compGreen(gsmr(i),work1,work1%nrow)
      call destroy(work1)
    end do

  end subroutine Make_gsmr_mem_dns

  !------------------------------------------------------------------------------------------------!

  !> g_small left (gsml) calculation
  !> Iterative calculation of the retarded GF from start block (sbl) to end block (ebl>sbl).
  SUBROUTINE Make_gsml_mem_dns(ESH,sbl,ebl)

    IMPLICIT NONE

    type(z_DNS), dimension(:,:), intent(in) :: ESH                    ! dense matrices array ESH(nbl,nbl)
    integer, intent(in) :: sbl,ebl                                    ! block indexes           

    ! Out gsmr(:) is global module variable !!!

    !Work
    TYPE(z_DNS) :: work1, work2
    INTEGER :: nrow
    INTEGER :: i, nbl
    !TYPE(z_DNS) :: INV(sbl,sbl)

    if (sbl.gt.ebl) return

    !***
    !gsml(sbl)
    !***
    nbl = size(ESH,1)
    nrow=ESH(sbl,sbl)%nrow  

    CALL create(gsml(sbl),nrow,nrow)

    call compGreen(gsml(sbl),ESH(sbl,sbl),nrow)


    DO i=sbl+1,ebl

      nrow=ESH(i,i)%nrow   !indblk(i+1)-indblk(i)

      CALL prealloc_mult(ESH(i,i-1),gsml(i-1),(-1.d0, 0.d0),work1)

      CALL prealloc_mult(work1,ESH(i-1,i),work2)

      CALL destroy(work1)

      CALL prealloc_sum(ESH(i,i),work2,work1)

      CALL destroy(work2)

      call create(gsml(i),work1%nrow,work1%nrow)

      call compGreen(gsml(i),work1,work1%nrow)

      CALL destroy(work1)

    ENDDO

    if (debug) then
      WRITE(*,*) '********************'
      WRITE(*,*) 'Make_gsml_mem done'
      WRITE(*,*) '********************'
    endif

  END SUBROUTINE Make_gsml_mem_dns

  !------------------------------------------------------------------------------------------------!

  !> Diagonal, Subdiagonal, Superdiagonal blocks of Green Retarded Gr(nbl,nbl).  
  subroutine Make_Gr_mem_dns(ESH,sbl,ebl)

    ! Input:
    ! ESH: dense matrices array ESH(nbl,nbl)
    ! sbl, ebl : block indexes
    ! If only sbl is specified, it calculates Gr(sbl, sbl)
    ! If sbl > ebl, it calculates Gr(i,i), Gr(i+1,i), Gr(i,i+1) for i=sbl,ebl (need gsml)          
    ! If sbl < ebl, it calculates Gr(i,i), Gr(i-1,i), Gr(i,i-1) for i=sbl,ebl (need gsmr)          
    !
    ! Output:
    ! dense matrices array global module variable Gr(i,j), i,i=1,nbl
    ! single blocks are allocated internally, array Gr(nbl,nbl) must be allocated externally

    implicit none 

    ! In
    type(z_DNS), dimension(:,:) :: ESH
    integer :: sbl
    integer, optional :: ebl

    ! Out Gr(:,:) is global module variable !!!
    
    ! Work
    integer :: i,nrow,nbl
    type(z_DNS) :: work1, work2, work3

    nbl = size(ESH,1)

    if (sbl.gt.nbl) return 
    if (sbl.lt.1) return

    !--------------------------------------------------------------------------!
    ! One diagonal block only
    !--------------------------------------------------------------------------!

    if (.not.present(ebl)) then
      if (nbl.eq.1) then
        nrow = ESH(sbl,sbl)%nrow     
        call create(Gr(sbl,sbl),nrow,nrow)
        call compGreen(Gr(sbl,sbl),ESH(sbl,sbl),nrow)
      else
        nrow = ESH(sbl,sbl)%nrow     
        call create(work1,nrow,nrow)
        work1%val = ESH(sbl,sbl)%val
        if (sbl+1.le.nbl) then
          call prealloc_mult(ESH(sbl,sbl+1),gsmr(sbl+1),work2)
          call prealloc_mult(work2,ESH(sbl+1,sbl),work3)
          call destroy(work2)
          call prealloc_sum(work1,work3,(-1.d0, 0.d0),work2)
          call destroy(work3)
          work1%val = work2%val
          call destroy(work2)
        end if
        if (sbl-1.ge.1) then
          call prealloc_mult(ESH(sbl,sbl-1),gsml(sbl-1),work2)
          call prealloc_mult(work2,ESH(sbl-1,sbl),work3)
          call destroy(work2)
          call prealloc_sum(work1,work3,(-1.d0, 0.d0),work2)
          call destroy(work3)
          work1%val = work2%val
          call destroy(work2)
        end if
        call create(Gr(sbl,sbl),nrow,nrow)
        call compGreen(Gr(sbl,sbl),work1,nrow)
        call destroy(work1)
      end if
      return
    end if

    !--------------------------------------------------------------------------!
    ! Diagonal, Subdiagonal and Superdiagonal blocks
    !--------------------------------------------------------------------------!
    
    if ((ebl.ge.sbl).and.(ebl.gt.1).and.(sbl.gt.1)) then
      do i=sbl,ebl,1
        call prealloc_mult(gsmr(i),ESH(i,i-1),work1)
        call prealloc_mult(work1,Gr(i-1,i-1),(-1.d0,0.d0),Gr(i,i-1))
        call destroy(work1)
        call prealloc_mult(ESH(i-1,i),gsmr(i),work2)
        call prealloc_mult(Gr(i-1,i-1),work2,(-1.d0, 0.d0),Gr(i-1,i))
        call prealloc_mult(Gr(i,i-1),work2,(-1.d0,0.d0),work1)
        call destroy(work2)
        call prealloc_sum(gsmr(i),work1,Gr(i,i))
        call destroy(work1) 
      end do
    else
      do i=sbl,ebl,-1
        call prealloc_mult(gsml(i),ESH(i,i+1),work1)
        call prealloc_mult(work1,Gr(i+1,i+1),(-1.d0,0.d0),Gr(i,i+1))
        call destroy(work1)
        call prealloc_mult(ESH(i+1,i),gsml(i),work2)
        call prealloc_mult(Gr(i+1,i+1),work2,(-1.d0, 0.d0),Gr(i+1,i))
        call prealloc_mult(Gr(i,i+1),work2,(-1.d0,0.d0),work1)
        call destroy(work2)
        call prealloc_sum(gsml(i),work1,Gr(i,i))
        call destroy(work1) 
      end do
    end if

  end subroutine Make_Gr_mem_dns

  !------------------------------------------------------------------------------------------------!
  
  !> Calculates Green Retarded block column "n"
  subroutine Make_Grcol_mem_dns(ESH,n,indblk)

    ! Input:
    ! ESH: dense matrices array ESH(nbl,nbl)
    ! n: n umber of column to be calculated
    !
    ! global variables needed: Gr(i,j) diagonal, subadiagonal and superdiagonal,
    ! gsmr(:) for downgoing and gsml(:) for upgoing
    ! 
    ! Output:
    ! dense matrices array global module variable Gr(:,n)
    ! single blocks are allocated internally, array Gr(nbl,nbl) must be allocated externally

    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    INTEGER, intent(in) :: n
    INTEGER, DIMENSION(:), intent(in) :: indblk 

    !Work
    INTEGER :: i,nrow,ncol,nbl
    TYPE(z_DNS) :: work1
    REAL(dp) :: max

    nbl = size(ESH,1)

    IF (n.GT.nbl) THEN
      STOP 'Error in Make_Grcol_mem : n is greater than nbl'
    ENDIF

    !***************************************
    !  Downgoing (j>=n+2 && n<nbl-1)
    !
    !   G_j,n = -gR_jj T_j,j-1 G_j-1,n
    !
    !***************************************
    ncol=indblk(n+1)-indblk(n)

    IF (n.LT.(nbl-1)) THEN

      DO i=n+2,nbl

        nrow=indblk(i+1)-indblk(i)

        max=MAXVAL(ABS(Gr(i-1,n)%val))

        IF (max.GT.EPS) THEN

          CALL prealloc_mult(gsmr(i),ESH(i,i-1),(-1.d0, 0.d0),work1)
          CALL prealloc_mult(work1,Gr(i-1,n),Gr(i,n))
          CALL destroy(work1)
        ENDIF

        ! NOW WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
        !    CALL create(Gr(i,n),nrow,ncol)
        !    Gr(i,n)%val(:,:)=(0.d0,0.d0)
        ! ENDIF

      ENDDO

    ENDIF
    !*************************************
    !   Up-going (j<=n-2 && n>2)
    !
    !   G_j,n = -gL_jj T_j,j+1 G_j+1,n
    !
    !*************************************

    IF (n.GT.2) THEN

      DO i=n-2,1,(-1)

        nrow=indblk(i+1)-indblk(i)

        max=MAXVAL(ABS(Gr(i+1,n)%val))

        IF (max.GT.EPS) THEN

          CALL prealloc_mult(gsml(i),ESH(i,i+1),(-1.d0, 0.d0),work1)
          CALL prealloc_mult(work1,Gr(i+1,n),Gr(i,n))
          CALL destroy(work1)

        ENDIF

        ! NOW WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
        !   CALL create(Gr(i,n),nrow,ncol)
        !   Gr(i,n)%val(:,:)=(0.d0,0.d0)
        !ENDIF

      ENDDO

    ENDIF

    if (debug) then
      WRITE(*,*) '********************'
      WRITE(*,*) 'Make_Grcol_mem done column',n
      WRITE(*,*) '********************'
    endif

  END SUBROUTINE Make_Grcol_mem_dns

  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !----------------------------------------- Makes for Gn -----------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  
  !****************************************************************************
  !
  ! Calculate G_n contributions for all contacts (except reference)
  ! Writing on memory
  !
  !****************************************************************************
  
  SUBROUTINE Make_Gn_mem_dns(ESH,SelfEneR,frm,ref,struct,Gn)

    !******************************************************************************
    !Input:  
    !ESH(nbl,nbl): sparse matrices array ES-H
    !SelfEneR(ncont): sparse matrices array with contacts Self Energy 
    !frm(ncont): Fermi diistribution value for each contact 
    !ref:  reference contact 
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont), Gr(:,:) 
    !diagonal, subdiagonal, overdiagonal and in colums Gr(:,cb) where cb are
    !layers interacting with all contacts but collector
    !
    !Output:
    !Gl: sparse matrix containing G  contacts contribution
    !*******************************************************************************  

    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: Gn
    TYPE(z_DNS), DIMENSION(:), intent(in) :: SelfEneR
    TYPE(Tstruct_info), intent(in) :: struct
    REAL(dp), DIMENSION(:), intent(in) :: frm
    INTEGER, intent(in) :: ref 

    !Work
    Type(z_DNS) :: Gam
    TYPE(z_DNS) :: work1,Ga
    INTEGER :: i,j,cb
    INTEGER :: ncont, nbl
    INTEGER, DIMENSION(:), POINTER :: cblk
    COMPLEX(dp) :: frmdiff

    ncont = struct%num_conts    
    nbl = struct%num_PLs
    cblk => struct%cblk

    !*******************************************
    ! Contact Iteration
    !*******************************************
    DO j=1,ncont

      ! NOTE: this soubroutine uses prealloc_mult that performs 
      ! C = C + A*B 
      IF (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN

        cb=cblk(j) ! block corresponding to contact j

        CALL zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

        frmdiff = cmplx(frm(j)-frm(ref),0.d0,dp)
        ! Computation of Gl(1,1) = Gr(1,cb) Gam(cb) Ga(cb,1)

        if (allocated(Gr(1,cb)%val)) then
          CALL prealloc_mult(Gr(1,cb),Gam,work1)    
          CALL zdagger(Gr(1,cb),Ga)
          CALL prealloc_mult(work1,Ga,frmdiff,Gn(1,1))
          CALL destroy(work1, Ga)
        else
           Gn(1,1)%val=(0.d0,0.d0)
        endif

        ! Computation of all tridiagonal blocks
        DO i=2,nbl

          ! Computation of Gl(i,j) = Gr(i,cb) Gam(cb) Ga(cb,j)
          ! Both Gr(i,cb) and Gr(j,cb) must be non-zero
          if (Gr(i-1,cb)%nrow.gt.0 .and. Gr(i,cb)%nrow.gt.0) then
            CALL prealloc_mult(Gr(i-1,cb),Gam,work1)
            CALL zdagger(Gr(i,cb),Ga)
            CALL prealloc_mult(work1,Ga,frmdiff,Gn(i-1,i))
            CALL destroy(work1)

            CALL prealloc_mult(Gr(i,cb),Gam,work1)
            CALL prealloc_mult(work1,Ga,frmdiff,Gn(i,i))

            CALL destroy(work1, Ga)

            CALL prealloc_mult(Gr(i,cb),Gam,work1)
            CALL zdagger(Gr(i-1,cb),Ga)
            CALL prealloc_mult(work1,Ga,frmdiff,Gn(i,i-1))

            CALL destroy(work1, Ga)
          else
            Gn(i-1,i)%val=(0.d0,0.d0)
            Gn(i,i-1)%val=(0.d0,0.d0)
          endif

        ENDDO

        call destroy(Gam)

      ENDIF

    ENDDO

  END SUBROUTINE Make_Gn_mem_dns

  !****************************************************************************
  !
  ! Calculate G_n contributions for all contacts (except reference)
  ! This version computes Grcol on the fly
  !
  !****************************************************************************
  
  SUBROUTINE Make_Gn_mem_dns2(ESH,SelfEneR,frm,ref,struct,Gn)

    !******************************************************************************
    !Input:  
    !ESH(nbl,nbl): sparse matrices array ES-H
    !SelfEneR(ncont): sparse matrices array with contacts Self Energy 
    !frm(ncont): Fermi diistribution value for each contact 
    !ref:  reference contact 
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont), Gr(:,:) 
    !diagonal, subdiagonal, overdiagonal and in colums Gr(:,cb) where cb are
    !layers interacting with all contacts but collector
    !
    !Output:
    !Gl: sparse matrix containing G  contacts contribution
    !*******************************************************************************  

    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: Gn
    TYPE(z_DNS), DIMENSION(:), intent(in) :: SelfEneR
    TYPE(Tstruct_info), intent(in) :: struct
    REAL(dp), DIMENSION(:), intent(in) :: frm
    INTEGER, intent(in) :: ref 

    !Work
    Type(z_DNS) :: Gam
    TYPE(z_DNS) :: work1,Ga
    INTEGER :: i,j,cb
    INTEGER :: ncont, nbl
    INTEGER, DIMENSION(:), POINTER :: cblk
    COMPLEX(dp) :: frmdiff

    ncont = struct%num_conts    
    nbl = struct%num_PLs
    cblk => struct%cblk

    !*******************************************
    ! Contact Iteration
    !*******************************************
    DO j=1,ncont

      ! NOTE: this soubroutine uses prealloc_mult that performs 
      ! C = C + A*B 
      IF (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN

        cb=cblk(j) ! block corresponding to contact j

        CALL zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

        frmdiff = cmplx(frm(j)-frm(ref),0.d0,dp)
        ! Computation of Gl(1,1) = Gr(1,cb) Gam(cb) Ga(cb,1)
        if (Gr(1,cb)%nrow.gt.0) then
          CALL prealloc_mult(Gr(1,cb),Gam,work1)    
          CALL zdagger(Gr(1,cb),Ga)
          CALL prealloc_mult(work1,Ga,frmdiff,Gn(1,1))
          CALL destroy(work1, Ga)
        else
          Gn(1,1)%val=(0.d0,0.d0)
        endif

        ! Computation of all tridiagonal blocks
        DO i=2,nbl

          ! Computation of Gr(i,cb) assuming Gr(i-1,cb) exists
          ! Assume downgoing: i > cb
          IF (Gr(i-1,cb)%nrow.GT.0) THEN
            CALL prealloc_mult(gsmr(i),ESH(i,i-1),(-1.d0, 0.d0),work1)
            CALL destroy(gsmr(i))
            CALL prealloc_mult(work1,Gr(i-1,cb),Gr(i,cb))
            CALL destroy(work1)
            IF (MAXVAL(ABS(Gr(i,cb)%val)).lt.EPS) call destroy(Gr(i,cb))
          ENDIF

          ! Computation of Gl(i,j) = Gr(i,cb) Gam(cb) Ga(cb,j)
          ! Both Gr(i,cb) and Gr(j,cb) must be non-zero
          IF (Gr(i-1,cb)%nrow.gt.0 .and. Gr(i,cb)%nrow.gt.0) THEN 
            CALL prealloc_mult(Gr(i-1,cb),Gam,work1)
            CALL zdagger(Gr(i,cb),Ga)
            CALL prealloc_mult(work1,Ga,frmdiff,Gn(i-1,i))
            CALL destroy(work1)

            CALL prealloc_mult(Gr(i,cb),Gam,work1)
            CALL prealloc_mult(work1,Ga,frmdiff,Gn(i,i))

            CALL destroy(work1, Ga)

            CALL prealloc_mult(Gr(i,cb),Gam,work1)
            CALL zdagger(Gr(i-1,cb),Ga)
            CALL prealloc_mult(work1,Ga,frmdiff,Gn(i,i-1))

            CALL destroy(work1, Ga)
          ELSE
            Gn(i-1,i)%val=(0.d0,0.d0)
            Gn(i,i-1)%val=(0.d0,0.d0)
          ENDIF

          IF (Gr(i-1,cb)%nrow.gt.0) call destroy(Gr(i-1,cb))

        ENDDO

        call destroy(Gam)

      ENDIF

    ENDDO

  END SUBROUTINE Make_Gn_mem_dns2
  
  !****************************************************************************
  !>
  ! Calculate G_n contributions due to el-ph
  ! Writing on memory
  !
  !****************************************************************************
  SUBROUTINE Make_Gn_ph(pnegf,ESH,iter,Gn)

    TYPE(Tnegf), intent(in) :: pnegf
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: Gn
    INTEGER, intent(in) :: iter

    INTEGER, DIMENSION(:), POINTER :: indblk
    Type(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_n, sigma_blk
    Type(z_DNS) :: Ga, work1, work2, sigma_tmp
    INTEGER :: n, k, nbl, nrow, ierr, ii, jj, norbs, nblk, indstart, indend

    !print *, 'debug: subroutine make_Gn_ph is started'
    
    !! If this is the first scba cycle, there's nothing to do
    if (pnegf%inter%scba_iter .eq. 0) then
      return
    endif
    nbl = pnegf%str%num_PLs
    ALLOCATE(Sigma_ph_n(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'

    ! The block sigma n is made available from el-ph model
    ! Note: the elph models could not keep a copy and calculate it 
    ! on the fly. You have to rely on the local copy
    if (allocated(pnegf%inter)) then
    call pnegf%inter%get_sigma_n(Sigma_ph_n, pnegf%ie)
    end if

    !! Make the sigma_ph_n available.
    !! The exact way depends on the el-ph model
    !! I could do it on the fly to save some memory
    select case(pnegf%elph%model)
    case (0)
      continue
    case (1)
      stop
    case (2)
      stop "Deprecated"
    case (3)
      stop "Deprecated"

    case default
      write(*,*) 'Not yet implemented'
      stop 0

      !! old alex implementation, not active
      !if (iter .gt. 0) then
      !  call read_blkmat(Sigma_ph_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
      !else 
      !  call write_blkmat(Sigma_ph_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
      !endif
      !! old alex implementation, not active
    end select

    !! Calculate the diagonal and off diagonal (if needed) blocks of Gn
    !! in the assumption of diagonal self energy
    !! G(k,k) = Gr(k,i)Sigma_n(i,i)Ga(i,k)
    !! G(k,k+1) = Gr(k,i)Sigma_n(i,i)Ga(i,k+1)
    !! G(k,k-1) = Gr(k,i)Sigma_n(i,i)Ga(i,k-1)
    !! All the rows of Gr need to be available
    !print *, 'debug: Calculate the diagonal and off diagonal (if needed) blocks of Gn'
    DO n = 1, nbl-1
      DO k = 1, nbl
        if (Gr(n,k)%nrow.gt.0) then
          CALL zdagger(Gr(n,k),Ga)
          !print *, 'debug: Gr(n,k)%ncol=',Gr(n,k)%ncol,' Sigma_ph_n(k,k)%ncol=', &
          !     Sigma_ph_n(k,k)%ncol,' work1%ncol=',work1%ncol
          CALL prealloc_mult(Gr(n,k), Sigma_ph_n(k,k), work1)
          CALL prealloc_mult(work1, Ga, work2)
          ! Computing diagonal blocks of Gn(n,n)
          Gn(n,n)%val = Gn(n,n)%val + work2%val
          call destroy(work2,Ga)                        
        endif
        ! Computing blocks of Gn(n,n+1)
        ! Only if S is not identity: Gn is initialized on ESH therefore
        ! we need to check the number of rows (or column)
        if (Gr(n+1,k)%nrow.gt.0 .and. Gn(n,n+1)%nrow .gt. 0) then
          CALL zdagger(Gr(n+1,k),Ga)
          !print *, 'debug: work1%ncol=',work1%ncol,' Ga%ncol=',Ga%ncol,' work2%ncol=',work2%ncol
          CALL prealloc_mult(work1, Ga, work2)
          Gn(n,n+1)%val = Gn(n,n+1)%val + work2%val
          !call destroy(work1,work2,Ga)                                     !DAR 
          Gn(n+1,n)%val = conjg(transpose(Gn(n,n+1)%val))
        endif
        call destroy(work1,work2,Ga)                                        !DAR
      END DO
    END DO
    !print *, 'debug: check2'
    DO k = 1, nbl
      if (Gr(nbl,k)%nrow.gt.0) then
        CALL zdagger(Gr(nbl,k),Ga)
        CALL prealloc_mult(Gr(nbl,k), Sigma_ph_n(k,k), work1)
        CALL prealloc_mult(work1, Ga, work2)
        Gn(nbl,nbl)%val = Gn(nbl,nbl)%val + work2%val
        call destroy(work1,work2,Ga)

      endif
    END DO
    !print *, 'debug: End Gn calculation'
    !! End Gn calculation


    DO n = 1, nbl
      CALL destroy(Sigma_ph_n(n,n))
    END DO

    DEALLOCATE(Sigma_ph_n)

    !print *, 'debug: subroutine make_Gn_ph is finished'

  END SUBROUTINE Make_Gn_ph

  !> Calculates G_n=-iG< contributions due to many-body interactions.
  !  Writing on memory.
  !----------------------------------------------------------------------------!
  SUBROUTINE Make_Gn_mbngf(pnegf,ESH,iter,Gn)

    TYPE(Tnegf), intent(inout) :: pnegf
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: Gn
    INTEGER, intent(in) :: iter

    INTEGER, DIMENSION(:), POINTER :: indblk
    Type(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_n, sigma_blk
    Type(z_DNS) :: Ga, work1, work2, sigma_tmp
    INTEGER :: n, k, nbl, nrow, ierr, ii, jj, norbs, nblk, indstart, indend

    !print *, 'debug: subroutine make_Gn_ph is started'
    
    !! If this is the first scba cycle, there's nothing to do
    if (pnegf%inter%scba_iter .eq. 0) then
      return
    endif
    nbl = pnegf%str%num_PLs
    ALLOCATE(Sigma_ph_n(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'

    ! The block sigma n is made available from el-ph model
    ! Note: the elph models could not keep a copy and calculate it 
    ! on the fly. You have to rely on the local copy
    if (allocated(pnegf%inter)) then
    call pnegf%inter%get_sigma_n(Sigma_ph_n, pnegf%ie)
    end if
 
    call add_sigma_n(pnegf)                                                 !DAR
    
    !! Make the sigma_ph_n available.
    !! The exact way depends on the el-ph model
    !! I could do it on the fly to save some memory
    select case(pnegf%elph%model)
    case (0)
      continue
    case (1)
      stop
    case (2)
      stop "Deprecated"
    case (3)
      stop "Deprecated"

    case default
      write(*,*) 'Not yet implemented'
      stop 0

      !! old alex implementation, not active
      !if (iter .gt. 0) then
      !  call read_blkmat(Sigma_ph_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
      !else 
      !  call write_blkmat(Sigma_ph_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
      !endif
      !! old alex implementation, not active
    end select

    !! Calculate the diagonal and off diagonal (if needed) blocks of Gn
    !! in the assumption of diagonal self energy
    !! G(k,k) = Gr(k,i)Sigma_n(i,i)Ga(i,k)
    !! G(k,k+1) = Gr(k,i)Sigma_n(i,i)Ga(i,k+1)
    !! G(k,k-1) = Gr(k,i)Sigma_n(i,i)Ga(i,k-1)
    !! All the rows of Gr need to be available
    !print *, 'debug: Calculate the diagonal and off diagonal (if needed) blocks of Gn'
    DO n = 1, nbl-1
      DO k = 1, nbl
        if (Gr(n,k)%nrow.gt.0) then
          CALL zdagger(Gr(n,k),Ga)
          !print *, 'debug: Gr(n,k)%ncol=',Gr(n,k)%ncol,' Sigma_ph_n(k,k)%ncol=', &
          !     Sigma_ph_n(k,k)%ncol,' work1%ncol=',work1%ncol
          CALL prealloc_mult(Gr(n,k), Sigma_ph_n(k,k), work1)
          CALL prealloc_mult(work1, Ga, work2)
          ! Computing diagonal blocks of Gn(n,n)
          Gn(n,n)%val = Gn(n,n)%val + work2%val
          call destroy(work2,Ga)                        
        endif
        ! Computing blocks of Gn(n,n+1)
        ! Only if S is not identity: Gn is initialized on ESH therefore
        ! we need to check the number of rows (or column)
        if (Gr(n+1,k)%nrow.gt.0 .and. Gn(n,n+1)%nrow .gt. 0) then
          CALL zdagger(Gr(n+1,k),Ga)
          !print *, 'debug: work1%ncol=',work1%ncol,' Ga%ncol=',Ga%ncol,' work2%ncol=',work2%ncol
          CALL prealloc_mult(work1, Ga, work2)
          Gn(n,n+1)%val = Gn(n,n+1)%val + work2%val
          !call destroy(work1,work2,Ga)                                     !DAR 
          Gn(n+1,n)%val = conjg(transpose(Gn(n,n+1)%val))
        endif
        call destroy(work1,work2,Ga)                                        !DAR
      END DO
    END DO
    !print *, 'debug: check2'
    DO k = 1, nbl
      if (Gr(nbl,k)%nrow.gt.0) then
        CALL zdagger(Gr(nbl,k),Ga)
        CALL prealloc_mult(Gr(nbl,k), Sigma_ph_n(k,k), work1)
        CALL prealloc_mult(work1, Ga, work2)
        Gn(nbl,nbl)%val = Gn(nbl,nbl)%val + work2%val
        call destroy(work1,work2,Ga)

      endif
    END DO
    !print *, 'debug: End Gn calculation'
    !! End Gn calculation


    DO n = 1, nbl
      CALL destroy(Sigma_ph_n(n,n))
    END DO

    DEALLOCATE(Sigma_ph_n)

    !print *, 'debug: subroutine make_Gn_ph is finished'

  END SUBROUTINE Make_Gn_mbngf

  !----------------------------------------------------------------------------!

  subroutine add_sigma_n(negf)

    TYPE(Tnegf), intent(inout) :: negf

  end subroutine add_sigma_n
  
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------- blk2csr --------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!

  SUBROUTINE blk2csr(G,struct,P,Gcsr)

    TYPE(z_DNS), DIMENSION(:,:) :: G
    TYPE(Tstruct_info), intent(in) :: struct
    TYPE(z_CSR) :: Gcsr
    TYPE(z_CSR) :: P, G_sp

    INTEGER, DIMENSION(:), POINTER :: indblk, cblk
    INTEGER :: nbl, oldx, row, col, iy, ix, x, y, ii, jj, nrows

    nbl = struct%num_PLs
    cblk => struct%cblk
    indblk => struct%mat_PL_start
    nrows = struct%mat_PL_end(nbl)

    !create Gcsr with same pattern of P
    CALL create(Gcsr,P%nrow,P%ncol,P%nnz)
    Gcsr%rowpnt = P%rowpnt
    Gcsr%colind = P%colind
    Gcsr%nzval = (0.d0, 0.d0)

    !Cycle upon all rows
    x = 1
    DO ii = 1, nrows
      !Search block x containing row ii
      oldx = x
      IF (oldx.EQ.nbl) THEN 
        x = oldx
      ELSE
        DO ix = oldx, oldx+1
          IF ( (ii.GE.indblk(ix)).AND.(ii.LT.indblk(ix+1)) ) x = ix
        ENDDO
      ENDIF

      !Offset: row is the index for separate blocks
      row = ii - indblk(x) + 1

      !Cycle upon columns of Gcsr (which has been ALREADY MASKED by S)
      DO jj = Gcsr%rowpnt(ii), Gcsr%rowpnt(ii+1) -1
        IF (Gcsr%colind(jj).gt.nrows) CYCLE
        !Choose which block column we're dealing with
        y = 0
        IF (x.eq.1) then
          IF ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then 
            y = 1
          ELSEIF ( (Gcsr%colind(jj).GE.indblk(x + 1)).AND.(Gcsr%colind(jj).LT.indblk(x + 2)) ) then 
            y = 2
          ENDIF
        elseif (x.eq.nbl) then
          IF ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then 
            y = nbl
          ELSEIF ( (Gcsr%colind(jj).GE.indblk(x - 1)).AND.(Gcsr%colind(jj).LT.indblk(x)) ) then 
            y = nbl - 1
          ENDIF
        else
          DO iy = x-1, x+1
            if ( (Gcsr%colind(jj).GE.indblk(iy)).AND.(Gcsr%colind(jj).LT.indblk(iy + 1)) ) y = iy
          ENDDO
        ENDIF

        IF (y.EQ.0) THEN
          write(*,*)     
          write(*,*) 'ERROR in blk2csr: probably wrong PL size', x
          write(*,*) 'row',ii,nrows,Gcsr%colind(jj)
          write(*,*) 'block indeces:',indblk(1:nbl)
          STOP
        ENDIF

        col = Gcsr%colind(jj) - indblk(y) + 1

        IF (allocated(G(x,y)%val)) THEN
          Gcsr%nzval(jj) = Gcsr%nzval(jj) + G(x,y)%val(row,col)
        ENDIF

      ENDDO

    ENDDO


  END SUBROUTINE blk2csr

  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !--------------------------------------------- Outer --------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!

  !****************************************************************************
  !
  !  Calculate Green Retarded in the 
  !  contacts regions, where overlap with device orbitals is non-zero
  !  (only the upper part, needed in both K=0 and K points calculations,
  !   or both upper and lower parts)
  !  writing on memory
  !
  !****************************************************************************

  SUBROUTINE Outer_Gr_mem_dns(Tlc,Tcl,gsurfR,struct,lower,Aout)

    !****************************************************************************
    !Input:
    !Tlc: sparse arraymatrices array containing contact-device interacting blocks 
    !     ESH-H
    !Tlc: sparse arraymatrices array containing device-contact interacting blocks 
    !     ESH-H
    !gsurfR: sparse matrices array containing contacts surface green
    !lower: if .true., also lower parts are calculated and concatenated
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont), 
    !cindblk(ncont), Gr in the diagonal block related to the interaction layer
    ! 
    !Output:
    !Aout: sparse matrix containing density matrix in the region 
    !      corresponding to non-zero overlap
    !
    !****************************************************************************

    IMPLICIT NONE 

    !In/Out
    TYPE(z_DNS), DIMENSION(:) :: Tlc,Tcl,gsurfR
    LOGICAL :: lower
    TYPE(Tstruct_info), intent(in) :: struct
    TYPE(z_CSR) :: Aout


    !Work
    TYPE(z_DNS) :: work1, Grcl, Grlc 
    TYPE(z_CSR) :: GrCSR, TCSR
    INTEGER :: i,cb,nrow_tot,i1,j1
    INTEGER :: ncont, nbl
    INTEGER, DIMENSION(:), POINTER :: indblk, cblk

    ncont = struct%num_conts
    nbl = struct%num_PLs
    nrow_tot = struct%total_dim
    indblk => struct%mat_PL_start
    cblk => struct%cblk

    IF (.not.allocated(Aout%nzval)) THEN
      CALL create(Aout,nrow_tot,nrow_tot,0)
      Aout%rowpnt(:)=1
    ENDIF

    DO i=1,ncont

      !Numero di blocco del contatto
      cb=cblk(i)
      CALL prealloc_mult(Gr(cb,cb),Tlc(i),(-1.d0, 0.d0),work1)
      CALL prealloc_mult(work1,gsurfR(i),Grlc)

      CALL destroy(work1)

      j1=nzdrop(Grlc,EPS) 
      CALL create(GrCSR,Grlc%nrow,Grlc%ncol,j1)
      CALL dns2csr(Grlc,GrCSR)
      CALL destroy(Grlc)
      j1=nzdrop(Tlc(i),EPS)
      CALL create(TCSR,Tlc(i)%nrow,Tlc(i)%ncol,j1)
      CALL dns2csr(Tlc(i),TCSR)
      CALL zmask_realloc(GrCSR,TCSR)
      CALL destroy(TCSR)

      !Concatenazione di Asub nella posizione corrispondente
      i1=indblk(cb)
      j1=struct%mat_B_start(i)

      CALL concat(Aout,GrCSR,i1,j1)   

      CALL destroy(GrCSR)

      IF (lower) THEN

        CALL prealloc_mult(gsurfR(i),Tcl(i),(-1.d0, 0.d0), work1)
        CALL prealloc_mult(work1, Gr(cb,cb), Grcl)

        CALL destroy(work1)

        j1=nzdrop(Grcl,EPS) 
        CALL create(GrCSR,Grcl%nrow,Grcl%ncol,j1)
        CALL dns2csr(Grcl,GrCSR)
        CALL destroy(Grcl)
        j1=nzdrop(Tcl(i),EPS)
        CALL create(TCSR,Tcl(i)%nrow,Tcl(i)%ncol,j1)
        CALL dns2csr(Tcl(i),TCSR)
        CALL zmask_realloc(GrCSR,TCSR)
        CALL destroy(TCSR)

        i1 = struct%mat_B_start(i)-struct%central_dim+indblk(nbl+1)-1
        j1 = indblk(cb)

        CALL concat(Aout,GrCSR,i1,j1)

        CALL destroy(GrCSR)

      ENDIF

    ENDDO

    if (debug) then
      WRITE(*,*) '********************'
      WRITE(*,*) 'Outer_GreenR_mem done'
      WRITE(*,*) '********************'
    endif

  END SUBROUTINE Outer_Gr_mem_dns



  !****************************************************************************
  !
  !  Calculate Gn contributions for all contacts except collector, in the
  !  outer region where contacts-device overlap is non-zero - writing on 
  !  memory
  !
  !****************************************************************************

  SUBROUTINE Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,struct,frm,ref,lower,Glout)

    !****************************************************************************
    !Input:
    !Tlc: sparse matrices array containing contacts-device interaction blocks (ES-H)
    !gsurfR: sparse matrices array containing contacts surface green
    !SelfEneR: sparse matrices array containing contacts Self Energy
    !frm: array containing Fermi distribution values for all contacts
    !ref: reference contact
    !
    !global variables needed: nbl, indblk(nbl+1), cindblk(ncont), ncont, 
    !cblk(ncont), Gr(:,:), diagonal, subdiagonal, overdiagonal and 
    !in colums Gr(:,cb) where cb are layers interacting with all contacts 
    !but collector 
    !
    !Output:
    !Glout: sparse matrix containing G_n  contributions in the region 
    !       corresponding to non-zero overlap
    !
    !****************************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:) :: Tlc, gsurfR, SelfEneR
    REAL(dp), DIMENSION(:) :: frm
    TYPE(Tstruct_info), intent(in) :: struct
    INTEGER :: ref 
    LOGICAL :: lower
    TYPE(z_CSR) :: Glout

    !Work
    TYPE(z_DNS) :: Gam, gsurfA, Ga, work1, work2, work3, Glsub
    TYPE(z_CSR) :: GlCSR, TCSR
    INTEGER :: j,k,cb,cbj,i1,j1,nrow_tot
    INTEGER :: ncont, nbl
    INTEGER, DIMENSION(:), POINTER :: indblk, cblk
    COMPLEX(dp) :: frmdiff


    ncont = struct%num_conts
    nbl = struct%num_PLs
    nrow_tot = struct%total_dim
    indblk => struct%mat_PL_start
    cblk => struct%cblk

    !Allocazione della Glout
    !Righe totali del conduttore effettivo nrow_tot 
    !nrow_tot=indblk(nbl+1)-1
    !DO i=1,ncont
    !   nrow_tot=nrow_tot+ncdim(i) !gsurfR(i)%nrow
    !ENDDO
    IF (.not.allocated(Glout%nzval)) THEN
      CALL create(Glout,nrow_tot,nrow_tot,0)
      Glout%rowpnt(:)=1
    ENDIF
    !***
    !Iterazione su tutti i contatti "k" 
    !***
    DO k=1,ncont

      !Esegue le operazioni relative al contatto solo se e` valida la condizione
      !sulle distribuzioni di Fermi e se non si tratta del contatto iniettante (ref)
      IF ((ABS(frm(k)-frm(ref)).GT.EPS).AND.(k.NE.ref)) THEN

        !Calcolo della Gamma corrispondente
        cb=cblk(k)
        !nota:cb indica l'indice di blocco corrispondente al contatto j-esimo

        CALL zspectral(SelfEneR(k),SelfEneR(k),0,Gam)

        !***
        !Calcolo del contributo sulla proria regione
        !***          
        frmdiff=cmplx(frm(ref)-frm(k))

        CALL zspectral(gsurfR(k),gsurfR(k),0,work1)

        !print*, 'work2=Tlc*work1=Tlc*j(gsurfR-gsurfA)'
        CALL prealloc_mult(Tlc(k),work1,work2)
        CALL destroy(work1)

        !print *, 'work1=-Gr(cb,cb)*work2=-Gr(cb,cb)*Tlc*j(gsurfR-gsurfA)'
        CALL prealloc_mult(Gr(cb,cb),work2,frmdiff,work1)
        CALL destroy(work2)

        !print*,'work2=Tlc*gsurfA'
        CALL zdagger(gsurfR(k),gsurfA)
        CALL prealloc_mult(Tlc(k),gsurfA,work2)
        CALL destroy(gsurfA)

        !print*,'work3=Ga*work2=Ga*Tlc*gsurfA'           
        CALL zdagger(Gr(cb,cb),Ga)
        CALL prealloc_mult(Ga,work2,work3)

        CALL destroy(Ga)
        CALL destroy(work2)

        !print*,'work2=Gam*work3=Gam*Ga*Tlc*gsurfA'          
        CALL prealloc_mult(Gam,work3,work2)
        CALL destroy(work3)

        !print *,'work3=-Gr*work2=-Gr*Gam*Ga*Tlc*gsurfA'        
        CALL prealloc_mult(Gr(cb,cb),work2,frmdiff,work3)
        CALL destroy(work2)

        !Contributo totale sulla propria regione
        CALL prealloc_sum(work3,work1,Glsub)
        CALL destroy(work1)
        CALL destroy(work3)

        call mask(Glsub,Tlc(k)) 
        i1=nzdrop(Glsub,EPS) 

        IF (i1.gt.0) THEN
          CALL create(GlCSR,Glsub%nrow,Glsub%ncol,i1)
          CALL dns2csr(Glsub,GlCSR)

          !Concatenazione di Glsub nella matrice globale Glout
          i1=indblk(cb)
          j1=struct%mat_B_start(k)-struct%central_dim+indblk(nbl+1)-1
          CALL concat(Glout,GlCSR,i1,j1)

          ! compute lower outer part using (iG<)+ = iG<    
          IF (lower) THEN
            call zdagger(GlCSR,TCSR)
            call concat(Glout,TCSR,j1,i1)
            call destroy(TCSR)
          ENDIF

          CALL destroy(GlCSR)
        END IF

        CALL destroy(Glsub)
        !***
        !Ciclo per il calcolo dei contributi sulle regioni degli altri contatti
        !***

        DO j=1,ncont

          cbj=cblk(j)
          !Esegue le operazioni del ciclo solo se il j.ne.k o se
          !il blocco colonna di Gr e` non nullo (altrimenti il contributo e` nullo) 

          IF ((j.NE.k).AND.(Gr(cbj,cb)%nrow.NE.0 .AND. (Gr(cbj,cb)%ncol.NE.0))) THEN

            !print*,'work1=Tlc*gsurfA'  
            CALL zdagger(gsurfR(j),gsurfA)
            CALL prealloc_mult(Tlc(j),gsurfA,work1)
            CALL destroy(gsurfA)

            !print*,'work2=Ga*work1=Ga*Tlc*gsurfA'  
            CALL zdagger(Gr(cbj,cb),Ga)
            CALL prealloc_mult(Ga,work1,work2)

            CALL destroy(Ga)
            CALL destroy(work1)

            !print*,'work1=Gam*work2=Gam*Ga*Tlc*gsurfA'  
            CALL prealloc_mult(Gam,work2,work1)
            CALL destroy(work2)

            !print*,'Glsub=-Gr*work1=-Gr*Gam*Ga*Tlc*gsurfA'  
            CALL prealloc_mult(Gr(cbj,cb),work1,frmdiff,Glsub)
            CALL destroy(work1)

            call mask(Glsub,Tlc(j)) 
            i1=nzdrop(Glsub,EPS) 

            IF (i1.gt.0) THEN
              CALL create(GlCSR,Glsub%nrow,Glsub%ncol,i1)
              CALL dns2csr(Glsub,GlCSR)

              !Concatenazione di Glsub nella posizione corrispondente al contatto "j"
              i1=indblk(cbj)
              j1=struct%mat_B_start(j)-struct%central_dim+indblk(nbl+1)-1

              CALL concat(Glout,GlCSR,i1,j1)

              ! compute lower outer part using (iG<)+ = iG<    
              IF (lower) THEN
                call zdagger(GlCSR,TCSR)
                call concat(Glout,TCSR,j1,i1)
                call destroy(TCSR)
              ENDIF

              CALL destroy(GlCSR)
            ENDIF

            CALL destroy(Glsub)

          ENDIF
        ENDDO

        CALL destroy(Gam)

      ENDIF

    ENDDO

    if (debug) then
      WRITE(*,*) '********************'
      WRITE(*,*) 'Outer_Gl_mem done'
      WRITE(*,*) '********************'
    endif

  END SUBROUTINE Outer_Gn_mem_dns

  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------- allocate -------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!

  subroutine allocate_gsm_dns(gsm,nbl)
    
    type(z_DNS), dimension(:), allocatable :: gsm 
    integer :: nbl, ierr

    allocate(gsm(nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate gsm'

  end subroutine allocate_gsm_dns

  !---------------------------------------------------

  subroutine allocate_blk_dns(blkM,nbl)
    
    type(z_DNS), dimension(:,:), allocatable :: blkM 
    integer :: nbl, ierr

    allocate(blkM(nbl,nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate block-Matrix'

  end subroutine allocate_blk_dns

  !---------------------------------------------------

  subroutine deallocate_gsm_dns(gsm)
    
    type(z_DNS), dimension(:), allocatable :: gsm 
    integer :: ierr

    deallocate(gsm,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate gsmr'

  end subroutine deallocate_gsm_dns

  !---------------------------------------------------

  subroutine deallocate_blk_dns(blkM)
    
    type(z_DNS), dimension(:,:), allocatable :: blkM 
    integer :: ierr

    deallocate(blkM,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate block-Matrix'

  end subroutine deallocate_blk_dns

  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------!

END MODULE tranas_ngf_iterative


