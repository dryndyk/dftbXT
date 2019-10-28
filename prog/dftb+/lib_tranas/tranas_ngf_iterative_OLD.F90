
  !public :: sigma_ph_n
  !public :: sigma_ph_p
  !public :: sigma_ph_r
  !public :: sigma_ph_r_z
  !public :: check_sigma_ph_r
  !public :: check_Gl_Gr
  !public :: complete_sigma_ph_r

  !public :: create_scratch
  !public :: destroy_scratch

  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: GGn
  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: GGp
  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: GGr
  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: Sgn
  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: Sgp
  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: Sgr
  LOGICAL, PARAMETER :: memory = .true.

  ! READ Matrices
  SUBROUTINE read_blkmat(Matrix, path, name, i, j, iE)
    TYPE(z_DNS), intent(inout) :: Matrix
    CHARACTER(*), intent(in) :: path
    CHARACTER(*), intent(in) :: name
    INTEGER, intent(in) :: i, j, iE

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    logical :: lex
    integer :: i1,i2
    complex(dp) :: mat_el

    if (memory) then
      select case(trim(name))
      case('G_n_')     
        Matrix%val = GGn(i,j,iE)%val 
      case('G_r_')     
        Matrix%val = GGr(i,j,iE)%val
      case('G_p_')     
        Matrix%val = GGp(i,j,iE)%val
      case('Sigma_ph_r_')     
        Matrix%val = Sgr(i,j,iE)%val 
      case('Sigma_ph_n_')     
        Matrix%val = Sgn(i,j,iE)%val
      case('Sigma_ph_p_')     
        Matrix%val = Sgp(i,j,iE)%val
      case default 
        stop 'internal error: read_blkmat does not correspond'
      end select
      return
    endif

    Matrix%val = (0.d0,0.d0)

    if (i.le.9999) write(ofblki,'(i4.4)') i
    if (i.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (j.le.9999) write(ofblkj,'(i4.4)') j 
    if (j.gt.9999) stop 'ERROR: too many blks (> 9999)'

    if (iE.le.99999) write(ofpnt,'(i5.5)') iE

    filename = trim(name)//trim(ofblki)//'_'//trim(ofblkj)//'_'//trim(ofpnt)//'.dat'

    inquire(file=trim(path)//trim(filename),EXIST=lex)
    if (.not.lex) then
      RETURN
      !WRITE(*,*) 'ERROR: FILE '//trim(filename)//' DOES NOT EXIST'
      !STOP
    endif

    open(9091,file=trim(path)//trim(filename), access='STREAM')

    call inmat_c(9091,.false.,Matrix%val,Matrix%nrow,Matrix%ncol)

    !open(9001,file=trim(path)//trim(filename), access='DIRECT', recl=4)
    !call direct_in_c(9001,Matrix%val,Matrix%nrow)

    close(9091)

  END SUBROUTINE read_blkmat

  ! WRITE Matrices
  SUBROUTINE write_blkmat(Matrix, path, name, i, j, iE)
    TYPE(z_DNS), intent(in) :: Matrix
    CHARACTER(*), intent(in) :: path
    CHARACTER(*), intent(in) :: name
    INTEGER, intent(in) :: i, j, iE

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    integer :: m,n

    if (memory) then

      select case(trim(name))
      case('G_n_')     
        if (.not.allocated(GGn(i,j,iE)%val)) then
          call create(GGn(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGn(i,j,iE)%val = Matrix%val
      case('G_r_')     
        if (.not.allocated(GGr(i,j,iE)%val)) then
          call create(GGr(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGr(i,j,iE)%val = Matrix%val
      case('G_p_')     
        if (.not.allocated(GGp(i,j,iE)%val)) then
          call create(GGp(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGp(i,j,iE)%val = Matrix%val 
      case('Sigma_ph_r_')     
        if (.not.allocated(Sgr(i,j,iE)%val)) then
          call create(Sgr(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgr(i,j,iE)%val = Matrix%val
      case('Sigma_ph_n_')     
        if (.not.allocated(Sgn(i,j,iE)%val)) then
          call create(Sgn(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgn(i,j,iE)%val = Matrix%val
      case('Sigma_ph_p_')     
        if (.not.allocated(Sgp(i,j,iE)%val)) then
          call create(Sgp(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgp(i,j,iE)%val = Matrix%val
      case default 
        stop 'internal error: write_blkmat does not correspond'
      end select
      return
    endif

    if (i.le.9999) write(ofblki,'(i4.4)') i
    if (i.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (j.le.9999) write(ofblkj,'(i4.4)') j
    if (j.gt.9999) stop 'ERROR: too many blks (> 9999)'

    if (iE.le.99999) write(ofpnt,'(i5.5)') iE

    filename = trim(name)//trim(ofblki)//'_'//trim(ofblkj)//'_'//trim(ofpnt)//'.dat'

    open(9001,file=trim(path)//trim(filename), access='STREAM', status='REPLACE')

    call outmat_c(9001,.false.,Matrix%val,Matrix%nrow,Matrix%ncol) !,1.0d-36)

    !open(9001,file=trim(path)//trim(filename), status='REPLACE', access='DIRECT', recl=4)
    !call direct_out_c(9001,Matrix%val,Matrix%nrow)

    close(9001)

  END SUBROUTINE write_blkmat

  !****************************************************************************
  !
  !  Calculate Gn contributions for all contacts but collector, in the 
  !  contacts regions, where overlap with device orbitals is non-zero
  !  writing on memory
  !
  !****************************************************************************

  !SUBROUTINE Outer_A_mem_dns(Tlc,Tcl,gsurfR,struct,Aout)

  !****************************************************************************
  !Input:
  !Tlc: sparse arraymatrices array containing contact-device interacting blocks 
  !     ESH-H
  !Tlc: sparse arraymatrices array containing device-contact interacting blocks 
  !     ESH-H
  !gsurfR: sparse matrices array containing contacts surface green
  !
  !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont), 
  !cindblk(ncont), Gr in the diagonal block related to the interaction layer
  ! 
  !Output:
  !Aout: sparse matrix containing density matrix in the region 
  !      corresponding to non-zero overlap
  !
  !****************************************************************************

  !    IMPLICIT NONE 

  !In/Out
  !    TYPE(z_CSR), DIMENSION(:) :: Tlc,Tcl,gsurfR
  !    TYPE(Tstruct_Info), intent(in) :: struct
  !    TYPE(z_CSR), intent(out) :: Aout

  !Work
  !    TYPE(z_CSR) :: work1,GrCSR
  !    TYPE(z_CSR) :: Asub, Grlc, Grcl
  !    INTEGER :: i,cb,i1,j1
  !    INTEGER :: ncont,nrow_tot,nbl
  !    INTEGER, DIMENSION(:), POINTER :: indblk, cblk

  !    ncont = struct%num_conts
  !    nbl = struct%num_PLs
  !    nrow_tot = struct%total_dim
  !    indblk => struct%mat_PL_start
  !    cblk => struct%cblk

  !Righe totali del conduttore effettivo nrow_tot 
  !nrow_tot=indblk(nbl+1)-1
  !DO i=1,ncont
  !   nrow_tot=nrow_tot+ncdim(i)   !gsurfR(i)%nrow
  !ENDDO
  !    CALL create(Aout,nrow_tot,nrow_tot,0)
  !    Aout%rowpnt(:)=1



  !    DO i=1,ncont

  !Numero di blocco del contatto
  !       cb=cblk(i)

  !converto Gr(cb,cb) da denso a sparso
  !       call create(GrCSR,Gr(cb,cb)%nrow,Gr(cb,cb)%ncol,Gr(cb,cb)%nrow*Gr(cb,cb)%ncol)
  !       call dns2csr(Gr(cb,cb),GrCSR)

  !Calcolo di Grlc
  !       CALL prealloc_mult(GrCSR,Tlc(i),(-1.d0, 0.d0),work1)
  !Nota: numero colonne di Tlc = numero di righe di gsurf(i)

  !       CALL prealloc_mult(work1,gsurfR(i),Grlc)
  !       CALL destroy(work1)

  !Calcolo di Grcl
  !       CALL prealloc_mult(gsurfR(i),Tcl(i),(-1.d0, 0.d0),work1)

  !       CALL prealloc_mult(work1,GrCSR,Grcl)
  !       CALL destroy(work1)

  !Calcolo della spectral density
  !       CALL zspectral(Grlc,Grcl,0,Asub)

  !Dellocazione delle Green Retarded corrispondenti
  !       CALL destroy(Grlc)
  !       CALL destroy(Grcl)

  !Concatenazione di Asub nella posizione corrispondente
  !       i1=indblk(cb)
  !       j1=struct%mat_B_start(i)-struct%central_dim+indblk(nbl+1)-1
  !       CALL concat(Aout,Asub,i1,j1)
  !       CALL destroy(Asub)  

  !       call destroy(GrCSR)

  !    ENDDO



  !    if (debug) then
  !       WRITE(*,*) '********************'
  !       WRITE(*,*) 'Outer_A_mem done'
  !       WRITE(*,*) '********************'
  !    endif

  !  END SUBROUTINE Outer_A_mem_dns

!****************************************************************************
!
!  Calculate Spectral Density - writing on memory
!
!****************************************************************************

!  SUBROUTINE Make_Spectral_mem(struct,A)

!****************************************************************************
!Input:
!nrow_tot: total device matrix rows
!
!global variable needed: nbl, indblk(nbl+1), Gr(:,:) 
!
!Output:
!A: sparse matrix containing spectral density of device (allocated internally)
!****************************************************************************

!    IMPLICIT NONE

!In/Out
!    TYPE(z_CSR) :: A
!    TYPE(Tstruct_info) :: struct

!Work
!    INTEGER :: i,nrow,nrow_prev,i1,j1,ierr,nrow_tot,nbl
!    TYPE(z_CSR), DIMENSION(:,:), ALLOCATABLE :: Asub
!    INTEGER, DIMENSION(:), POINTER :: indblk

!    nbl = struct%num_PLs
!    indblk => struct%mat_PL_start
!    nrow_tot = struct%total_dim


!Allocazione dell'array di sparse Asub e di A
!    ALLOCATE(Asub(nbl,nbl),stat=ierr)
!    IF (ierr.NE.0) THEN
!       STOP 'ALLOCATION ERROR: could not allocate Asub(nbl,nbl)'
!    ENDIF

!    CALL create(A,nrow_tot,nrow_tot,0)
!    A%rowpnt(:)=1

!***
!A(1,1) 
!***
!    nrow=indblk(2)-indblk(1)

!    CALL zspectral(Gr(1,1),Gr(1,1),0,Asub(1,1))

!    CALL concat(A,Asub(1,1),1,1)
!    CALL destroy(Asub(1,1))

!***
!Diagonal, Subdiagonal and Super diagonal blocks
!***

!    DO i=2,nbl

!       nrow=indblk(i+1)-indblk(i)

!       i1=indblk(i)
!       j1=indblk(i)

!       CALL zspectral(Gr(i,i),Gr(i,i),0,Asub(i,i))

!       CALL concat(A,Asub(i,i),i1,j1)
!       CALL destroy(Asub(i,i))

!       i1=indblk(i-1)
!       j1=indblk(i)

!       CALL zspectral(Gr(i-1,i),Gr(i,i-1),0,Asub(i-1,i))
!       CALL concat(A,Asub(i-1,i),i1,j1)
!       CALL destroy(Asub(i-1,i))  

!       i1=indblk(i)
!       j1=indblk(i-1)

!       CALL zspectral(Gr(i,i-1),Gr(i-1,i),0,Asub(i,i-1))

!Blocks concatenation
!       CALL concat(A,Asub(i,i-1),i1,j1)
!       CALL destroy(Asub(i,i-1))

!    ENDDO

!    DEALLOCATE(Asub)

!    if (debug) then    
!       WRITE(*,*) '**********************'
!       WRITE(*,*) 'Make_Spectral_mem done'
!       WRITE(*,*) '**********************'
!    endif

!  END SUBROUTINE Make_Spectral_mem





!****************************************************************************
!
!  Calculate Green Retarded - writing on memory
!
!****************************************************************************

!!$  SUBROUTINE Make_GreenR_mem(nrow_tot,A)
!!$
!!$    !****************************************************************************
!!$    !Input:
!!$    !nrow_tot: total device matrix rows
!!$    !
!!$    !global variable needed: nbl, indblk(nbl+1), Gr(:,:) 
!!$    !
!!$    !Output:
!!$    !A: sparse matrix containing Green Retarded of device (allocated internally)
!!$    !****************************************************************************
!!$
!!$    IMPLICIT NONE
!!$
!!$    !In/Out
!!$    INTEGER :: nrow_tot
!!$    TYPE(z_CSR) :: A
!!$
!!$    !TYPE(z_CSR), DIMENSION(nbl,nbl) :: ESH
!!$
!!$    !Work
!!$    INTEGER :: i,nrow,i1,j1,ierr
!!$
!!$    CALL create(A,nrow_tot,nrow_tot,0)
!!$    A%rowpnt(:)=1
!!$
!!$    !write(*,*) 'A created'
!!$
!!$    !***
!!$    !A(1,1)
!!$    !***
!!$    nrow=indblk(2)-indblk(1)
!!$
!!$    CALL concat(A,Gr(1,1),1,1)
!!$
!!$    !***
!!$    !Diagonal, Subdiagonal and Superdiagonal blocks
!!$    !***
!!$    DO i=2,nbl
!!$
!!$       nrow=indblk(i+1)-indblk(i)
!!$
!!$       !       write(*,*) 'nrow=',nrow
!!$
!!$       i1=indblk(i)
!!$       j1=indblk(i)
!!$
!!$       CALL concat(A,Gr(i,i),i1,j1)
!!$
!!$       !       write(*,*) 'Gr',i,i,'concat'
!!$
!!$       i1=indblk(i-1)
!!$       j1=indblk(i)
!!$
!!$       CALL concat(A,Gr(i-1,i),i1,j1)
!!$
!!$       !       write(*,*) 'Gr',i-1,i,'concat'
!!$
!!$       i1=indblk(i)
!!$       j1=indblk(i-1)
!!$
!!$       CALL concat(A,Gr(i,i-1),i1,j1)
!!$
!!$       !       write(*,*) 'Gr',i,i-1,'concat'
!!$
!!$
!!$       !       write(*,*) 'Gr dealloc'
!!$
!!$    ENDDO
!!$
!!$    !if (debug) call writePeakInfo(6)    
!!$    if (debug) then
!!$       WRITE(*,*) '**********************'
!!$       WRITE(*,*) 'Make_GreenR_mem done'
!!$       WRITE(*,*) '**********************'
!!$    endif
!!$
!!$  END SUBROUTINE Make_GreenR_mem

  ! ******************************************************************************
  ! Computes Sigma_ph_n and save it file
  ! ******************************************************************************
  SUBROUTINE Sigma_ph_n(pnegf,Epnt)

    TYPE(Tnegf) :: pnegf
    REAL(dp), DIMENSION(:),allocatable :: Epnt

    !Local variables
    TYPE(z_DNS) :: work1,work2
    TYPE(z_DNS) :: Mq_mat
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_n, G_n_interP, G_n_interN

    REAL(dp), DIMENSION(:), POINTER :: Wq
    REAL(dp), DIMENSION(:), POINTER :: Nq
    REAL(dp), DIMENSION(:), POINTER :: Mq
    LOGICAL, DIMENSION(:), POINTER :: selmodes

    INTEGER :: i, m, iE, i1,i2
    INTEGER :: nummodes, numselmodes, nbl
    INTEGER, DIMENSION(:), pointer :: indblk
    REAL(dp) :: E1, E2, En

    nbl = pnegf%str%num_PLs
    indblk => pnegf%str%mat_PL_start

    selmodes => pnegf%elph%selmodes
    Wq => pnegf%elph%Wq
    Mq => pnegf%elph%Mq
    Nq => pnegf%elph%Nq
    numselmodes = pnegf%elph%numselmodes
    nummodes = pnegf%elph%nummodes

    iE = pnegf%iE

    allocate(G_n_interP(nbl,nbl))
    allocate(G_n_interN(nbl,nbl))
    allocate(Sigma_n(nbl,nbl))

    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_n(i,i), m, m)
      Sigma_n(i,i)%val = (0.d0, 0.d0)
      call create(G_n_interP(i,i), m, m)
      call create(G_n_interN(i,i), m, m)
    enddo

    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle
      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,Wq(m),Epnt,i1,i2,E1,E2)

      En = real(pnegf%Epnt)+Wq(m)

      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
          'G_n_', G_n_interP)

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,-Wq(m),Epnt,i1,i2,E1,E2)

      En = real(pnegf%Epnt)-Wq(m)

      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
          'G_n_', G_n_interN)


      do i = 1, nbl

        i1 =  G_n_interN(i,i)%nrow

        call create(work1, i1, i1)       

        work1%val = (Nq(m)+1.0_dp)*G_n_interP(i,i)%val + Nq(m)* G_n_interN(i,i)%val

        call create_id(Mq_mat,i1,Mq(m))

        call prealloc_mult( Mq_mat, work1, work2)

        call destroy(work1)

        call prealloc_mult( work2, Mq_mat, work1)

        call destroy(work2)

        Sigma_n(i,i)%val = Sigma_n(i,i)%val + work1%val 

        call destroy(work1, Mq_mat)

      end do

    end do

    !Throw away all non diagonal parts
    !do i = 1, nbl
    !   i1 = Sigma_n(i,i)%nrow
    !   do m = 1, i1 
    !      Sigma_n(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp) 
    !      Sigma_n(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp) 
    !   enddo
    !enddo
    !

    do i = 1, nbl
      call write_blkmat(Sigma_n(i,i),pnegf%scratch_path,'Sigma_ph_n_',i,i,iE)
      call destroy(Sigma_n(i,i))
      call destroy(G_n_interN(i,i))
      call destroy(G_n_interP(i,i))
    enddo

    deallocate(Sigma_n,G_n_interN, G_n_interP)

  END SUBROUTINE Sigma_ph_n

  
  ! ******************************************************************************
  ! Computes Sigma_ph> and save it file
  ! ******************************************************************************
  SUBROUTINE Sigma_ph_p(pnegf,Epnt)

    TYPE(Tnegf) :: pnegf
    REAL(dp), DIMENSION(:),allocatable :: Epnt

    !Local variables
    TYPE(z_DNS) :: work1,work2
    TYPE(z_DNS) :: Mq_mat
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_p, G_p_interP, G_p_interN

    REAL(dp), DIMENSION(:), POINTER :: Wq
    REAL(dp), DIMENSION(:), POINTER :: Nq
    REAL(dp), DIMENSION(:), POINTER :: Mq
    LOGICAL, DIMENSION(:), POINTER :: selmodes

    INTEGER :: i, m, iE, i1,i2
    INTEGER :: nummodes, numselmodes, nbl
    INTEGER, DIMENSION(:), pointer :: indblk
    REAL(dp) :: E1, E2, En

    nbl = pnegf%str%num_PLs
    indblk => pnegf%str%mat_PL_start

    selmodes => pnegf%elph%selmodes
    Wq => pnegf%elph%Wq
    Mq => pnegf%elph%Mq
    Nq => pnegf%elph%Nq
    numselmodes = pnegf%elph%numselmodes
    nummodes = pnegf%elph%nummodes

    iE = pnegf%iE

    allocate(G_p_interP(nbl,nbl))
    allocate(G_p_interN(nbl,nbl))
    allocate(Sigma_p(nbl,nbl))

    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_p(i,i), m, m)
      Sigma_p(i,i)%val = (0.d0, 0.d0)
      call create(G_p_interP(i,i), m, m)
      call create(G_p_interN(i,i), m, m)
    enddo

    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle
      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,Wq(m),Epnt,i1,i2,E1,E2)

      En = real(pnegf%Epnt)+Wq(m)

      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
          'G_p_', G_p_interP)

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,-Wq(m),Epnt,i1,i2,E1,E2)
      En = real(pnegf%Epnt)-Wq(m)
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
          'G_p_', G_p_interN)

      do i = 1, nbl

        !print*
        !print*,'(sigma_r) G_p+',maxval(abs(G_p_interP(i,i)%val))
        !print*,'(sigma_r) G_p-',maxval(abs(G_p_interN(i,i)%val))

        i1 =  G_p_interN(i,i)%nrow

        call create(work1, i1, i1)       

        work1%val = (Nq(m)+1.0_dp)*G_p_interN(i,i)%val + Nq(m)* G_p_interP(i,i)%val

        call create_id(Mq_mat,i1,Mq(m))

        call prealloc_mult( Mq_mat, work1, work2)

        call destroy(work1)

        call prealloc_mult ( work2, Mq_mat, work1)

        call destroy(work2)

        Sigma_p(i,i)%val = Sigma_p(i,i)%val + work1%val 

        call destroy(work1, Mq_mat)

      end do

    end do

    ! throw away all non diagonal parts
    !do i = 1, nbl
    !   i1 = Sigma_p(i,i)%nrow
    !   do m = 1, i1 
    !      Sigma_p(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp) 
    !      Sigma_p(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp) 
    !   enddo
    !enddo
    !

    do i = 1, nbl
      call write_blkmat(Sigma_p(i,i),pnegf%scratch_path,'Sigma_ph_p_',i,i,iE)
      call destroy(Sigma_p(i,i))
      call destroy(G_p_interN(i,i))
      call destroy(G_p_interP(i,i))
    enddo

    deallocate(Sigma_p,G_p_interN, G_p_interP)

  END SUBROUTINE Sigma_ph_p

  !--------------------------------------------------------------------------------
  SUBROUTINE Sigma_ph_r(pnegf,Epnt)

    TYPE(Tnegf) :: pnegf
    REAL(dp), DIMENSION(:),allocatable :: Epnt

    !Local variables
    TYPE(z_DNS) :: work1,work2
    TYPE(z_DNS) :: Mq_mat
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r, G_r_interP, G_r_interN
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: G_n_interP, G_n_interN

    REAL(dp), DIMENSION(:), POINTER :: Wq
    REAL(dp), DIMENSION(:), POINTER :: Nq
    REAL(dp), DIMENSION(:), POINTER :: Mq
    LOGICAL, DIMENSION(:), POINTER :: selmodes
    INTEGER :: i, m, iE, i1, i2
    INTEGER :: nummodes, numselmodes, nbl
    INTEGER, DIMENSION(:), pointer :: indblk
    REAL(dp) :: E1, E2, En

    nbl = pnegf%str%num_PLs
    indblk => pnegf%str%mat_PL_start


    selmodes => pnegf%elph%selmodes
    Mq => pnegf%elph%Mq
    Wq => pnegf%elph%Wq
    Nq => pnegf%elph%Nq
    nummodes = pnegf%elph%nummodes
    numselmodes = pnegf%elph%numselmodes

    iE=pnegf%iE

    allocate(G_r_interP(nbl,nbl))
    allocate(G_r_interN(nbl,nbl))
    allocate(G_n_interP(nbl,nbl))
    allocate(G_n_interN(nbl,nbl))
    allocate(Sigma_r(nbl,nbl))
    !allocate(G_r(nbl,nbl))



    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_r(i,i), m, m)
      Sigma_r(i,i)%val = (0.d0, 0.d0)
      call create(G_r_interP(i,i), m, m)
      call create(G_r_interN(i,i), m, m)
      call create(G_n_interP(i,i), m, m)
      call create(G_n_interN(i,i), m, m)
      !call create(G_r(i,i), m, m)
      !call read_blkmat(G_r(i,i),pnegf%scratch_path,'G_r_',i,i,iE)
    enddo

    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,Wq(m),Epnt,i1,i2,E1,E2)
      En = real(pnegf%Epnt)+Wq(m)
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
          'G_n_', G_n_interP)
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
          'G_r_', G_r_interP)

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,-Wq(m),Epnt,i1,i2,E1,E2)
      En = real(pnegf%Epnt)-Wq(m)
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
          'G_n_', G_n_interN)
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
          'G_r_', G_r_interN)

      do i = 1, nbl
        !print*
        !print*,'(sigma_r) G_r+',maxval(abs(G_r_interP(i,i)%val))
        !print*,'(sigma_r) G_r-',maxval(abs(G_r_interN(i,i)%val))

        !print*,'(sigma_r) G_n+',maxval(abs(G_n_interP(i,i)%val))
        !print*,'(sigma_r) G_n-',maxval(abs(G_n_interN(i,i)%val))

        i1 = Sigma_r(i,i)%nrow

        call create(work1,i1,i1)

        work1%val = (0.d0,0.d0)

        if (pnegf%elph%selfene_gr) then
          !print*,'SelfEneR_Gr'    
          ! Via I: should be exact for E >> Ef_max
          work1%val = (Nq(m)+1.0_dp)*G_r_interN(i,i)%val + Nq(m)*G_r_interP(i,i)%val
          ! Via II: should be exact for E << Ef_min
          !work1%val = (Nq(m)+1.0_dp)*G_r_interP(i,i)%val + Nq(m)*G_r_interN(i,i)%val
          ! Via III: should work as a compromise
          !work1%val = Nq(m)*(G_r_interP(i,i)%val + G_r_interN(i,i)%val) + G_r(i,i)%val
        endif

        if (pnegf%elph%selfene_gless) then
          !print*,'SelfEneR_G<'    
          ! should be 1/2 [G<- - G<+] == i/2 [Gn- - Gn+]    
          work1%val = work1%val + &
              (0.0_dp,0.5_dp) * (G_n_interN(i,i)%val - G_n_interP(i,i)%val)
        endif

        if (pnegf%elph%selfene_gless .or. pnegf%elph%selfene_gr) then

          call create_id(Mq_mat,i1,Mq(m))

          call prealloc_mult(Mq_mat, work1, work2)

          call destroy(work1)

          call prealloc_mult(work2, Mq_mat, work1)

          call destroy(work2)

          Sigma_r(i,i)%val = Sigma_r(i,i)%val + work1%val 

          call destroy(work1, Mq_mat)

        else

          Sigma_r(i,i)%val = (0.d0,0.d0) 
          call destroy(work1)

        endif

      enddo

    enddo

    ! throw away all non diagonal parts
    !do i = 1, nbl
    !    i1 = Sigma_r(i,i)%nrow
    !    do m = 1, i1 
    !       Sigma_r(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp) 
    !       Sigma_r(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp) 
    !    enddo
    ! enddo
    !

    do i = 1, nbl
      !print*          
      !print*,'(sigma_r) Sigma_ph_r',maxval(abs(Sigma_r(i,i)%val)), iE
      call write_blkmat(Sigma_r(i,i),pnegf%scratch_path,'Sigma_ph_r_',i,i,iE)
      call destroy(Sigma_r(i,i))
      call destroy(G_r_interP(i,i))
      call destroy(G_r_interN(i,i))
      call destroy(G_n_interP(i,i))
      call destroy(G_n_interN(i,i))
      !call destroy(G_r(i,i))
    enddo

    deallocate(Sigma_r,G_r_interP,G_r_interN,G_n_interP,G_n_interN)
    !deallocate(G_r)

  END SUBROUTINE Sigma_ph_r

  SUBROUTINE Sigma_ph_r_z(pnegf,z)

    TYPE(Tnegf) :: pnegf
    COMPLEX(dp), intent(in) :: z

    !Local variables
    TYPE(z_DNS) :: work1,work2
    TYPE(z_DNS) :: Mq_mat
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r, G_r
    REAL(dp), DIMENSION(:), POINTER :: Wq
    REAL(dp), DIMENSION(:), POINTER :: Nq
    REAL(dp), DIMENSION(:), POINTER :: Mq
    LOGICAL, DIMENSION(:), POINTER :: selmodes
    INTEGER :: i, m, iE,i1
    INTEGER :: nummodes, numselmodes, nbl
    INTEGER, DIMENSION(:), pointer :: indblk

    nbl = pnegf%str%num_PLs
    indblk => pnegf%str%mat_PL_start


    selmodes => pnegf%elph%selmodes
    Mq => pnegf%elph%Mq
    Wq => pnegf%elph%Wq
    Nq => pnegf%elph%Nq
    nummodes = pnegf%elph%nummodes
    numselmodes = pnegf%elph%numselmodes

    iE=pnegf%iE

    allocate(Sigma_r(nbl,nbl))
    allocate(G_r(nbl,nbl))

    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_r(i,i), m, m)
      Sigma_r(i,i)%val = (0.d0, 0.d0)
      call create(G_r(i,i), m, m)
    enddo

    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle

      do i = 1, nbl

        call read_blkmat(G_r(i,i),pnegf%scratch_path,'G_r_',i,i,iE)

        i1 = Sigma_r(i,i)%nrow

        if (pnegf%elph%selfene_gr) then

          call create(work1,i1,i1)

          work1%val = (0.d0,0.d0)
          !print*,'SelfEneR_Gr'    
          work1%val = (2.0_dp*Nq(m)+1.0_dp)*G_r(i,i)%val

          call create_id(Mq_mat,i1,Mq(m))

          call prealloc_mult(Mq_mat, work1, work2)

          call destroy(work1)

          call prealloc_mult(work2, Mq_mat, work1)

          call destroy(work2)

          Sigma_r(i,i)%val = Sigma_r(i,i)%val +  work1%val 

          call destroy(work1, Mq_mat)

        else
          Sigma_r(i,i)%val = (0.d0,0.d0) 
        endif

      enddo

    enddo

    do i = 1, nbl
      !print*,'(sigma_r) Sigma_ph_r', pnegf%iE, maxval(abs(Sigma_r(i,i)%val))
      call write_blkmat(Sigma_r(i,i),pnegf%scratch_path,'Sigma_ph_r_',i,i,iE)
      call destroy(Sigma_r(i,i))
      call destroy(G_r(i,i))
    enddo

  END SUBROUTINE Sigma_ph_r_z
  ! ----------------------------------------------------------
  SUBROUTINE check_Gl_Gr(pnegf)
    TYPE(Tnegf) :: pnegf
    REAL(dp), DIMENSION(:),allocatable :: Epnt
    integer :: ioffset

    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: G_r, G_p, G_n
    TYPE(z_DNS) :: A, T 
    integer :: nbl, n, sizebl, i_start, i_stop, psize, maxpos(2)
    real(dp) :: Wmax, maxdev, tmp, maxG
    INTEGER, DIMENSION(:), pointer :: indblk

    nbl = pnegf%str%num_PLs
    indblk => pnegf%str%mat_PL_start

    allocate(G_n(nbl,nbl))
    allocate(G_p(nbl,nbl))
    allocate(G_r(nbl,nbl))

    maxdev = 0.0_dp
    maxG = 0.0_dp
    psize = 0

    do n = 1, nbl
      sizebl = indblk(n+1)-indblk(n)
      call create(G_r(n,n),sizebl,sizebl)
      call create(G_n(n,n),sizebl,sizebl)
      call create(G_p(n,n),sizebl,sizebl)
      call read_blkmat(G_r(n,n), pnegf%scratch_path, 'G_r_',n,n, pnegf%iE)
      call read_blkmat(G_n(n,n),pnegf%scratch_path,'G_n_',n,n,pnegf%iE)
      call read_blkmat(G_p(n,n),pnegf%scratch_path,'G_p_',n,n,pnegf%iE)

      call zspectral(G_r(n,n),G_r(n,n),0,A)

      call create(T,sizebl,sizebl)

      !print*,'CHECK Sigma_ph_n+Sigma_ph_p',maxval(abs(Sigma_n(n,n)%val+Sigma_p(n,n)%val))

      T%val = G_n(n,n)%val + G_p(n,n)%val - A%val 

      tmp = maxval(abs(A%val))
      if (tmp .gt. maxG) maxG=tmp

      tmp = maxval(abs(T%val))/maxval(abs(A%val)) 
      if (tmp .gt. maxdev) then
        maxdev = tmp
        maxpos = maxloc(abs(T%val)) + psize
      endif

      !print*

      psize = psize + sizebl

      call destroy(G_r(n,n))
      call destroy(G_n(n,n))
      call destroy(G_p(n,n))
      call destroy(A)
      call destroy(T)

    enddo

    print*,'CHECK Gn+Gp=Gr-Ga',pnegf%iE, maxG, maxdev

    deallocate(G_n, G_p, G_r)

  END SUBROUTINE check_Gl_Gr

  ! ----------------------------------------------------------
  SUBROUTINE check_sigma_ph_r(pnegf)
    TYPE(Tnegf) :: pnegf

    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r, Sigma_p, Sigma_n
    TYPE(z_DNS) :: Gam, T 
    integer :: nbl, n, sizebl, psize, maxpos(2)
    real(dp) :: maxdev, tmp, maxG
    INTEGER, DIMENSION(:), pointer :: indblk

    nbl = pnegf%str%num_PLs
    indblk => pnegf%str%mat_PL_start

    allocate(Sigma_n(nbl,nbl))
    allocate(Sigma_p(nbl,nbl))
    allocate(Sigma_r(nbl,nbl))

    maxdev = 0.0_dp
    maxG = 0.0_dp
    psize = 0

    do n = 1, nbl
      sizebl = indblk(n+1)-indblk(n)
      call create(Sigma_r(n,n),sizebl,sizebl)
      call create(Sigma_n(n,n),sizebl,sizebl)
      call create(Sigma_p(n,n),sizebl,sizebl)
      call read_blkmat(Sigma_r(n,n), pnegf%scratch_path, 'Sigma_ph_r_',n,n, pnegf%iE)
      call read_blkmat(Sigma_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
      call read_blkmat(Sigma_p(n,n),pnegf%scratch_path,'Sigma_ph_p_',n,n,pnegf%iE)

      call zspectral(Sigma_r(n,n),Sigma_r(n,n),0,Gam)

      call create(T,sizebl,sizebl)

      !print*,'CHECK Sigma_ph_n+Sigma_ph_p',maxval(abs(Sigma_n(n,n)%val+Sigma_p(n,n)%val))

      T%val = Sigma_n(n,n)%val + Sigma_p(n,n)%val - Gam%val 

      tmp = maxval(abs(Gam%val))
      if (tmp .gt. maxG) maxG=tmp

      tmp = maxval(abs(T%val))/maxval(abs(Gam%val)) 
      if (tmp .gt. maxdev) then
        maxdev = tmp
        maxpos = maxloc(abs(T%val)) + psize
      endif

      !print*

      psize = psize + sizebl

      call destroy(Sigma_r(n,n))
      call destroy(Sigma_n(n,n))
      call destroy(Sigma_p(n,n))
      call destroy(Gam)
      call destroy(T)

    enddo

    print*,'CHECK Sigma_ph_r',pnegf%iE, maxG, maxdev

    deallocate(Sigma_n, Sigma_p, Sigma_r)


  END SUBROUTINE check_sigma_ph_r


  ! ----------------------------------------------------------
  ! Search points for interpolations. 
  ! Wq has a sign (+/- Wq)
  !
  SUBROUTINE search_points(pnegf, Wq, Epnt, i1, i2, E1, E2)
    use energy_mesh, only : elem 
    type(TNegf) :: pnegf
    real(dp) :: Wq                
    real(dp), dimension(:), allocatable :: Epnt
    integer, intent(out) :: i1,i2 
    real(dp), intent(out) :: E1, E2

    integer :: iE, iel, istart, iend, ip
    Type(elem), pointer :: pel
    real(dp) :: En

    if (Wq.eq.0) then
      i1 = pnegf%iE
      i2 = pnegf%iE
      E1 = pnegf%Epnt
      E2 = pnegf%Epnt
      return
    endif

    if (allocated(Epnt)) then
      ! Remove offset such that search can work on Epnt(1..N)      
      iE = pnegf%iE - pnegf%Np_n(1) - pnegf%Np_n(2) - pnegf%n_poles
      En = real(pnegf%Epnt) + Wq !Wq carry the right sign      
      !print*
      !print*,'iE', iE, real(pnegf%Epnt) + Wq 

      if (sign(1.0_dp,Wq) .gt. 0) then
        i2 = iE + 1
        if (i2.gt.size(Epnt)) then
          i2 = size(Epnt)
        else    
          do while (Epnt(i2) .lt. En)
            i2 = i2 + 1
          end do
        endif
        i1 = i2 - 1
      else
        i1 = iE - 1
        if (i1.lt.1) then
          i1 = 1 
        else
          do while (Epnt(i1) .gt. En)
            i1 = i1 - 1
          end do
        endif
        i2 = i1 + 1
      endif
      E1 = Epnt(i1)
      E2 = Epnt(i2)
      ! add back the offset to the point
      i1 = i1 + pnegf%Np_n(1) + pnegf%Np_n(2) + pnegf%n_poles
      i2 = i2 + pnegf%Np_n(1) + pnegf%Np_n(2) + pnegf%n_poles

    else  

      !if (.not.allocated(pnegf%emesh%pactive)) STOP 'emesh not initialized'

      En = real(pnegf%Epnt) + Wq       

      if (sign(1.0_dp,Wq) .gt. 0) then
        istart = pnegf%emesh%iactive
        iend = pnegf%emesh%maxind

        elloop1: do iel = istart, iend
          pel => pnegf%emesh%pactive(iel)%pelem
          do ip = 1, 3
            if (pel%pnt(ip) .gt. En) then
              exit elloop1 
            endif
          end do
        end do elloop1
        i1 = pel%map(ip-1)
        i2 = pel%map(ip)
        E1 = pel%pnt(ip-1)
        E2 = pel%pnt(ip)
      else
        istart = pnegf%emesh%iactive
        iend = 1

        elloop2: do iel = istart, iend, -1
          pel => pnegf%emesh%pactive(iel)%pelem
          do ip = 3, 1, -1
            if (pel%pnt(ip) .lt. En) then
              exit elloop2 
            endif
          end do
        end do elloop2
        i1 = pel%map(ip)
        i2 = pel%map(ip+1)
        E1 = pel%pnt(ip)
        E2 = pel%pnt(ip+1)
      end if
    end if
    !print*
    !print*,E1,En,E2
    !print*,'interpolate between:',i1,i2

  END SUBROUTINE search_points


  SUBROUTINE interpolation(i1,i2, E1, E2, E, path, name, G_interp)
    INTEGER, intent(in) :: i1, i2
    REAL(dp) :: E1, E2, E
    TYPE(z_DNS), DIMENSION(:,:) :: G_interp
    CHARACTER(*) :: path
    CHARACTER(*) :: name

    !local variables
    TYPE(z_DNS) :: work1,work2
    INTEGER :: i

    do i = 1, size(G_interp,1)

      call create(work1,G_interp(i,i)%nrow,G_interp(i,i)%ncol)
      work1%val = (0.d0,0.d0)
      call read_blkmat(work1, path, name, i, i, i1)

      if (E1.ne.E2 .and. i1.ne.i2) then

        call create(work2,G_interp(i,i)%nrow,G_interp(i,i)%ncol)
        work2%val = (0.d0,0.d0)
        call read_blkmat(work2, path, name, i, i, i2)

        G_interp(i,i)%val = ((E-E1)*work2%val + (E2-E)*work1%val)/(E2-E1)

        call destroy(work2)

      else

        G_interp(i,i)%val = work1%val

      endif

      call destroy(work1)

    end do


  END SUBROUTINE interpolation


  !*********************************************************************
  !ADD Hilbert-transform part to Sigma_ph_r  
  !Need to set back FFT transforms
  !********************************************************************
  SUBROUTINE complete_sigma_ph_r(pnegf, Epnt, ioffset)

    TYPE(Tnegf) :: pnegf
    REAL(dp), DIMENSION(:) :: Epnt
    INTEGER :: ioffset

    ! Locals
    INTEGER :: i,j,n,m, iE, bl, sizebl,i_start,i_stop, ierr, nummodes
    REAL(dp), DIMENSION(:), POINTER :: Wq
    REAL(dp), DIMENSION(:), POINTER :: Mq
    REAL(dp) :: Wmax, dE, tmp

    character(4) :: ofblki
    complex(dp), DIMENSION(:), allocatable :: temp1 
    type(z_DNS), DIMENSION(:), allocatable :: Gn_E, Sigma_r_E
    !type(z_DNS) :: work1, work2, work3, Mq_mat 
    LOGICAL, DIMENSION(:), POINTER :: selmodes

    Mq => pnegf%elph%Mq
    Wq => pnegf%elph%Wq
    selmodes => pnegf%elph%selmodes

    n = size(Epnt)
    nummodes = size(Wq)

    Wmax = maxval(Wq)
    dE = Epnt(2)-Epnt(1)
    ! ENERGY-INTERVAL FOR SELF-ENERGIES:
    i_start = aint(Wmax/dE) + 1
    i_stop =  (n - i_start) 
    i_start = i_start + 1

    ! PROBLEM: 
    ! The Hilb-Transf should be resticted on the sub interval
    ! Emin+m*Wmax ... Emax-m*Wmax
    ! However the number of points should be always 2^p
    ! This condition is difficult to fulfill.
    ! Possible solution: interpolation between grids.

    ! CREATE THE ARRAY OF MATRICES (one for each E-point)
    allocate(Gn_E(n),stat=ierr)
    allocate(Sigma_r_E(n),stat=ierr)
    if(ierr.ne.0) STOP 'ERROR in allocation of Gn_E'
    call log_allocate(temp1,n)

    print*
    print*,'HILBERT TRANSFORM memory:',sizebl*sizebl*n*16
    print*,'HILBERT TRANSFORM interval:',i_start,i_stop

    ! LOOP ON blocks
    do bl = 1, pnegf%str%num_PLs

      sizebl =  pnegf%str%mat_PL_start(bl+1) - pnegf%str%mat_PL_start(bl)
      if (bl.le.9999) write(ofblki,'(i4.4)') bl 
      if (bl.gt.9999) stop 'ERROR: too many blks (> 9999)'

      ! LOAD ALL G< FROM FILES and store in Gn_E(iE)
      do iE = 1, n
        call create(Gn_E(iE),sizebl,sizebl)
        call read_blkmat(Gn_E(iE), pnegf%scratch_path, 'G_n_', bl, bl,iE+ioffset)
      enddo

      ! LOAD ALL Sigma_r in the right interval
      do iE = i_start, i_stop
        call create(Sigma_r_E(iE),sizebl,sizebl)
        call read_blkmat(Sigma_r_E(iE), pnegf%scratch_path, 'Sigma_ph_r_', bl, bl, iE+ioffset)
      enddo

      do m = 1, nummodes 
        if (.not.selmodes(m)) cycle

        ! COMPUTE   Mq G< Mq   (assume diagonal now) 
        do iE = 1, n
          !call create_id(Mq_mat,sizebl,Mq(m))
          !call prealloc_mult(Mq_mat, Gn_E(iE), work1)
          !call prealloc_mult(work1, Mq_mat, work2)
          !call destroy(work1, Mq_mat)
          !Gn_E(iE)%val = work2%val
          !call destroy(work2)
          tmp = Mq(m)*Mq(m)
          Gn_E(iE)%val = tmp*Gn_E(iE)%val
        enddo

        ! PERFORM Hilbert in Energy for each (i,j)     
        do j = 1, sizebl  
          do i = 1, sizebl

            !          k = (j-1)*sizebl+i
            !          do iE = 1, n
            !             if (iE.le.99999) write(ofpnt,'(i5.5)') iE+ioffset
            !             filename = 'G_n_'//trim(ofblki)//'_'//trim(ofblki)//'_'//trim(ofpnt)//'.dat'
            !             open(100,file=trim(pnegf%scratch_path)//trim(filename), access='DIRECT', recl=4)
            !             READ(100, rec = k)  temp1(iE)
            !             close(100)
            !          enddo

            ! SETUP a vector out of all G<_ij(E)
            ! Here we could perform an efficient interpolation on a regular grid 2^p
            do iE = 1, n
              temp1(iE) = Gn_E(iE)%val(i,j)
            enddo

            Wq(m) = Wq(m)*2.0_dp*pi/((n-1)*dE)

            !call Hilbert_shift(temp1, Wq(m))

            Wq(m) = Wq(m)*((n-1)*dE)/(2.0_dp*pi)

            ! UPDATE the self-energies with the Hilbert part. 
            do iE = i_start, i_stop
              ! Should be  -i/2 H[ G<+ - G<-  ] = 1/2 H[ Gn+ - Gn- ]
              Sigma_r_E(iE)%val(i,j) =  Sigma_r_E(iE)%val(i,j) + (0.5_dp, 0.0)* temp1(iE) 
            enddo


            !          do iE = i_start+1, i_stop
            !             if (iE.le.99999) write(ofpnt,'(i5.5)') iE+ioffset
            !             filename = 'Sigma_ph_r_'//trim(ofblki)//'_'//trim(ofblki)//'_'//trim(ofpnt)//'.dat'
            !             open(100,file=trim(pnegf%scratch_path)//trim(filename), access='DIRECT', recl=4)
            !             READ (100, rec = k) temp2 
            !             temp2 = temp2 - (0.0_dp, 0.5_dp)* temp1(iE)   
            !             WRITE (100, rec = k) temp2
            !             close(100)
            !          enddo

          enddo ! Loop on block size   
        enddo ! Loop on block size 

      enddo !Loop on modes 

      do iE = i_start, i_stop
        call write_blkmat(Sigma_r_E(iE), pnegf%scratch_path, 'Sigma_ph_r_', bl, bl, iE+ioffset)
      enddo

      do iE = 1, n
        call destroy(Gn_E(iE))
      enddo
      do iE = i_start, i_stop
        call destroy(Sigma_r_E(iE))
      enddo

    enddo !Loop on blocks 

    call log_deallocate(temp1)
    deallocate(Gn_E)
    deallocate(Sigma_r_E)

  END SUBROUTINE complete_sigma_ph_r

  SUBROUTINE create_scratch(nbl, npoints)
    integer :: nbl, npoints

    integer :: i,j,k,err

    call destroy_scratch(nbl, npoints)

    ALLOCATE(GGn(nbl,nbl,npoints),stat=err)
    ALLOCATE(GGr(nbl,nbl,npoints),stat=err)
    ALLOCATE(GGp(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgn(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgr(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgp(nbl,nbl,npoints),stat=err)

    ! Initialize everything to 0 
    do k=1,npoints
      do j=1,nbl
        do i=1,nbl  
          GGn(i,j,k)%nrow=0
          GGn(i,j,k)%ncol=0
          GGr(i,j,k)%nrow=0
          GGr(i,j,k)%ncol=0
          GGp(i,j,k)%nrow=0
          GGp(i,j,k)%ncol=0
          Sgn(i,j,k)%nrow=0
          Sgn(i,j,k)%ncol=0
          Sgr(i,j,k)%nrow=0
          Sgr(i,j,k)%ncol=0
          Sgp(i,j,k)%nrow=0
          Sgp(i,j,k)%ncol=0
        enddo
      enddo
    enddo

    if(err.ne.0) then
      STOP 'ERROR: Cannot allocate GG'
    endif
    print*,'Created memory scratch',nbl,'x',nbl,'x',npoints

  END SUBROUTINE create_scratch

  SUBROUTINE destroy_scratch(nbl, npoints)
    integer :: nbl, npoints

    integer :: i, j, iE,  err

    err = 0


    do i = 1, nbl
      do j = 1, nbl
        do iE = 1, npoints

          if (allocated(GGn)) then
            if (allocated(GGn(i,j,iE)%val)) call destroy(GGn(i,j,iE))
          endif
          if (allocated(GGr)) then
            if (allocated(GGr(i,j,iE)%val)) call destroy(GGr(i,j,iE))
          endif
          if (allocated(GGp)) then
            if (allocated(GGp(i,j,iE)%val)) call destroy(GGp(i,j,iE))
          endif
          if (allocated(Sgn)) then
            if (allocated(Sgn(i,j,iE)%val)) call destroy(Sgn(i,j,iE))     
          endif
          if (allocated(Sgr)) then
            if (allocated(Sgr(i,j,iE)%val)) call destroy(Sgr(i,j,iE))
          endif
          if (allocated(Sgp)) then
            if (allocated(Sgp(i,j,iE)%val)) call destroy(Sgp(i,j,iE))     
          endif

        end do
      end do
    end do

    if (allocated(GGn)) DEALLOCATE(GGn,stat=err)
    if (allocated(GGr)) DEALLOCATE(GGr,stat=err)
    if (allocated(GGp)) DEALLOCATE(GGp,stat=err)
    if (allocated(Sgn)) DEALLOCATE(Sgn,stat=err)
    if (allocated(Sgr)) DEALLOCATE(Sgr,stat=err)
    if (allocated(Sgp)) DEALLOCATE(Sgp,stat=err)

    if(err.ne.0) then
      STOP 'ERROR: Cannot deallocate GG'
    endif

  END SUBROUTINE destroy_scratch
