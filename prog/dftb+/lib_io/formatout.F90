!--------------------------------------------------------------------------------------------------!
! DFTB+XT open software package for quantum nanoscale modeling (TraNaS OpenSuite)                  !
! Copyright (C) 2018-2020 Dmitry A. Ryndyk                                                         !
! DFTB+: general package for performing fast atomistic simulations                                 !
! Copyright (C) 2006-2020 DFTB+ developers group                                                   !
!--------------------------------------------------------------------------------------------------!
! GNU Lesser General Public License version 3 or (at your option) any later version.               !
! See the LICENSE file for terms of usage and distribution.                                        !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines for formatted output of data
module dftbp_formatout
  use dftbp_globalenv
  use dftbp_environment
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_message
  use dftbp_constants
  use dftbp_lapackroutines, only: matinv
  use dftbp_sparse2dense
  use dftbp_fileid !!DAR!!
  implicit none
  private

  public :: clearFile, writeGenFormat, writeXYZFormat
  public :: printDFTBHeader
  public :: writeSparseAsSquare, writeSparseAsSquare_old, writeSparse !!DAR!! + writeSparseAsSquare_old


  !> Clears contents of a file
  interface clearFile
    module procedure clearFile_fname
  end interface clearFile


  !> Writes geometry information in gen format to a file
  interface writeGenFormat
    module procedure writeGenFormat_fname
    module procedure writeGenFormat_fid
  end interface writeGenFormat


  !> Writes geometry information in xyz format to a file
  interface writeXYZFormat
    module procedure writeXYZFormat_fname
    module procedure writeXYZFormat_fid
  end interface writeXYZFormat


  !> Writes DFTB+ type sparse matrix in square form to disc
  interface writeSparseAsSquare
    module procedure writeSparseAsSquare_real
    module procedure writeSparseAsSquare_cplx
  end interface writeSparseAsSquare
  
  interface writeSparseAsSquare_old !!DAR!!
    module procedure writeSparseAsSquare_real_old
    module procedure writeSparseAsSquare_cplx_old
  end interface writeSparseAsSquare_old
  
contains


  !> Clears contents of file
  subroutine clearFile_fname(fileName)

    !> name of the file which should be cleared
    character(len=*), intent(in) :: fileName

    integer :: fd

    open(newunit=fd, file=fileName, status="replace", position="rewind")
    close(fd)

  end subroutine clearFile_fname


  !> A wrapper around writeGenFormat_fid to open a file first.
  subroutine writeGenFormat_fname(fileName, coord, species, speciesName, latVec, tFracCoord, append)

    !> File name of the file which should be created
    character(len=*), intent(in) :: fileName

    !> Coordinates in atomic units
    real(dp),         intent(in) :: coord(:,:)

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Name of the different species
    character(mc),    intent(in) :: speciesName(:)

    !> Lattice vectors
    real(dp), intent(in), optional :: latVec(3,3)

    !> Print out fractional coordinates?
    logical, intent(in), optional :: tFracCoord

    !> Whether geometry should be appended (default: it is overwritten)
    logical, intent(in), optional :: append

    integer :: fd

    logical :: append0

    if (present(append)) then
      append0 = append
    else
      append0 = .false.
    end if

    @:ASSERT((.not.(present(tFracCoord).neqv.present(latVec))) .or.(present(latVec)))

    if (append0) then
      open(newunit=fd, file=fileName, form="formatted", action="write", status="old",&
          & position="append")
    else
      open(newunit=fd, file=fileName, form="formatted", action="write", status="replace")
    end if
    call writeGenFormat(fd, coord, species, speciesName, latVec, tFracCoord)
    close(fd)

  end subroutine writeGenFormat_fname


  !> Writes coordinates in the famous GEN format to a file
  subroutine writeGenFormat_fid(fd, coord, species, speciesName, latVec, tFracCoord)

    !> File id of an open file where output should be written
    integer, intent(in) :: fd

    !> Coordinates in atomic units
    real(dp),          intent(in) :: coord(:,:)

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Name of the different species
    character(mc),     intent(in) :: speciesName(:)

    !> Lattice vectors
    real(dp), intent(in), optional :: latVec(:,:)

    !> Print out fractional coordinates?
    logical, intent(in), optional :: tFracCoord

    integer :: nAtom, nSpecies
    character(6) :: formatSpecies
    integer :: ii, jj
    logical :: tFractional
    real(dp) :: invLatVec(3,3)

100 format(I5," ",A2)
101 format("(",I2.2,"A3)")
102 format(I5,I2,3E20.10)
103 format(3E20.10)

    nAtom = size(coord, dim=2)
    nSpecies = maxval(species)

    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(size(species) == nAtom)
    @:ASSERT(size(speciesName) == nSpecies)
#:call ASSERT_CODE
    if (present(latVec)) then
      @:ASSERT(all(shape(latVec) == (/3, 3 /)))
    end if
#:endcall ASSERT_CODE
    @:ASSERT((.not.(present(tFracCoord).neqv.present(latVec))) .or.(present(latVec)))

    tFractional = .false.
    if (present(latVec)) then
      if (present(tFracCoord) ) then
        tFractional = tFracCoord
      end if
      if (tFractional) then
        write(fd, 100) nAtom, "F"
      else
        write(fd, 100) nAtom, "S"
      end if
    else
      write(fd, 100) nAtom, "C"
    end if
    write(formatSpecies, 101) nSpecies
    write(fd, formatSpecies) (trim(speciesName(ii)), ii = 1, nSpecies)

    if (tFractional) then
      invLatVec(:,:) = latVec(:,:)
      call matinv(invLatVec)
      do ii = 1, nAtom
        write(fd, 102) ii, species(ii), matmul(invLatVec,coord(:, ii))
      end do
    else
      do ii = 1, nAtom
        write(fd, 102) ii, species(ii), (coord(jj, ii) * Bohr__AA, jj = 1, 3)
      end do
    end if
    if (present(latVec)) then
      write(fd, 103) 0.0_dp, 0.0_dp, 0.0_dp
      do ii = 1, 3
        write(fd, 103) (latVec(jj, ii) * Bohr__AA, jj = 1, 3)
      end do
    end if
  end subroutine writeGenFormat_fid


  !> Writes coordinates in the XYZ format
  subroutine writeXYZFormat_fname(fileName, coord, species, speciesName, charges, velocities,&
      & comment, append)

    !> File name of a file to be created
    character(len=*), intent(in) :: fileName

    !> Coordinates in atomic units
    real(dp), intent(in) :: coord(:,:)

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Name of the different species
    character(mc), intent(in) :: speciesName(:)

    !> Optional vector with charges for each atom.
    real(dp), intent(in), optional :: charges(:)

    !> Optional array of velocity vectors for each atom.
    real(dp), intent(in), optional :: velocities(:,:)

    !> Optional comment for line 2 of the file
    character(len=*), intent(in), optional :: comment

    !> Whether geometry should be appended (default: it is overwritten)
    logical, intent(in), optional :: append

    integer :: fd
    logical :: append0

    if (present(append)) then
      append0 = append
    else
      append0 = .false.
    end if

    if (append0) then
      open(newunit=fd, file=fileName, action="write", form="formatted", status="old",&
          & position="append")
    else
      open(newunit=fd, file=fileName, action="write", form="formatted", status="replace")
    end if
    call writeXYZFormat(fd, coord, species, speciesName, charges, velocities, comment)
    close(fd)

  end subroutine writeXYZFormat_fname


  !> Writes coordinates in the XYZ format with additional charges and vectors
  subroutine writeXYZFormat_fid(fd, coords, species, speciesNames, charges, velocities, comment)

    !> File id of an open file where output should be written
    integer, intent(in) :: fd

    !> Coordinates in atomic units
    real(dp), intent(in) :: coords(:,:)

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Name of the different species
    character(mc), intent(in) :: speciesNames(:)

    !> Optional vector with charges for each atom.
    real(dp), intent(in), optional :: charges(:)

    !> Optional array of velocity vectors for each atom.
    real(dp), intent(in), optional :: velocities(:,:)

    !> Optional comment for line 2 of the file
    character(len=*), intent(in), optional :: comment

    integer :: nAtom, nSpecies
    integer :: ii, jj

200 format(I5)
201 format(A5,3F16.8)
202 format(A5,6F16.8)
203 format(A5,4F16.8)
204 format(A5,7F16.8)

    nAtom = size(coords, dim=2)
    nSpecies = maxval(species)
    @:ASSERT(size(coords, dim=1) == 3)
    @:ASSERT(size(species) == nAtom)
    @:ASSERT(size(speciesNames) == nSpecies)
#:call ASSERT_CODE
    if (present(charges)) then
      @:ASSERT(size(charges) == nAtom)
    end if
    if (present(velocities)) then
      @:ASSERT(all(shape(velocities) == (/ 3, nAtom /)))
    end if
#:endcall ASSERT_CODE

    write(fd, 200) nAtom
    if (present(comment)) then
      write(fd, "(A)") trim(comment)
    elseif (present(velocities)) then
      write(fd, *) "Velocity in AA/ps"
    else
      write(fd, *) ""
    end if

    if (present(charges) .and. present(velocities)) then
      write(fd, 204) (trim(speciesNames(species(ii))), coords(:, ii) * Bohr__AA, charges(ii),&
          & velocities(:,ii) * Bohr__AA / au__fs * 1000.0_dp, ii = 1, nAtom)
    elseif (present(charges) .and. .not. present(velocities)) then
      write(fd, 203) (trim(speciesNames(species(ii))), coords(:, ii) * Bohr__AA,&
          & charges(ii), ii = 1, nAtom)
    elseif (.not. present(charges) .and. present(velocities)) then
      write(fd, 202) (trim(speciesNames(species(ii))), coords(:, ii) * Bohr__AA,&
          & velocities(:,ii) * Bohr__AA / au__fs * 1000.0_dp, ii = 1, nAtom)
    else
      write(fd, 201) (trim(speciesNames(species(ii))),&
          & (coords(jj, ii) * Bohr__AA, jj = 1, 3), ii = 1, nAtom)
    end if

  end subroutine writeXYZFormat_fid


  !> Writes the greeting message of dftb+ on stdout
  subroutine printDFTBHeader(release, year, outunit)

    !> release version of the code
    character(len=*), intent(in) :: release

    !> release year
    integer, intent(in) :: year, outunit

    character, parameter :: vBar = '|'
    character, parameter :: hBar = '='
    integer, parameter :: headerWidth = 80

    write(outunit,*)
    write(outunit, '(A)') repeat(hBar, headerWidth)
    write(outunit, '(5A)') hBar,hBar, repeat(' ', headerWidth - 4),hBar,hBar
    write(outunit, '(5A)') hbar,hbar,'                              TraNaS OpenSuite                              ',hbar,hbar
    write(outunit, '(5A)') hbar,hbar,'          (integrated open software suite for nanoscale modeling)           ',hbar,hbar 
    write(outunit, '(5A)') hBar,hBar, repeat(' ', headerWidth - 4),hBar,hBar
    write(outunit, '(A)') repeat(hBar, headerWidth)
    write(outunit, '(5A)') hBar,hBar, repeat(' ', headerWidth - 4),hBar,hBar
    write(outunit, '(5A)') hBar,hBar,'                                DFTB+XT 1.03 (10.05.2020)                   ',hBar,hBar 
    write(outunit, '(5A)') hBar,hBar, repeat(' ', headerWidth - 4),hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' DFTB+ eXTended version for quantum nanoscale modeling                      ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' Copyright (C) 2006-2020 DFTB+ developers group                             ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' Copyright (C) 2018-2020 Dmitry A. Ryndyk                                   ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar, repeat(' ', headerWidth - 4),hBar,hBar
    write(outunit, '(A)') repeat(hBar, headerWidth)
    write(outunit, '(5A)') hBar,hBar, repeat(' ', headerWidth - 4),hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' Please cite as:                                                            ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' DFTB+XT code as a part of TraNaS OpenSuite [1], partially based on         ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' the DFTB+ software package [2,3] and the CP2K software package [4].        ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar, repeat(' ', headerWidth - 4),hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' [1] TraNaS OpenSuite, tranas.org/opensuite                                 ',hBar,hBar 
    write(outunit, '(5A)') hBar,hBar,' [2] B. Hourahine et al., J. Chem. Phys. 152, 124101 (2020)                 ',hBar,hBar 
    write(outunit, '(5A)') hBar,hBar,' [3] A. Pecchia et al., New Journal of Physics 10, 065022 (2008)            ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' [4] CP2K (Development Version), cp2k.org                                   ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar, repeat(' ', headerWidth - 4),hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' You should also cite additional publications crediting the parametrization ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' data you use. Please consult the documentation of the SK-files for the     ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' references. Other references to the particular methods are given in the    ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar,' manual and in the end of the output file.                                  ',hBar,hBar
    write(outunit, '(5A)') hBar,hBar, repeat(' ', headerWidth - 4),hBar,hBar
    write(outunit, '(A)') repeat(hBar, headerWidth)

  end subroutine printDFTBHeader


  !> Converts a sparse matrix to its square form and write it to a file.
  subroutine writeSparseAsSquare_real(env, fname, sparse, iNeighbour, nNeighbourSK, iAtomStart,&
      & iPair, img2CentCell)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Name of the file to write the matrix to.
    character(len=*), intent(in) :: fname

    !> Sparse matrix.
    real(dp), intent(in) :: sparse(:)

    !> Neighbour list index.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours.
    integer, intent(in) :: nNeighbourSK(:)

    !> Offset array in the square matrix.
    integer, intent(in) :: iAtomStart(:)

    !> Offset array in the sparse matrix
    integer, intent(in) :: iPair(0:,:)

    !> Pair indexing array.
    integer, intent(in) :: img2CentCell(:)

    !> Mapping of the atoms to the central cell.
    real(dp), allocatable :: square(:,:)
    character(mc) :: strForm
    integer :: fd, nOrb

    if (withMpi) then
      call error("Writing of HS not working with MPI yet")
    end if

    nOrb = iAtomStart(size(nNeighbourSK) + 1) - 1

    allocate(square(nOrb, nOrb))
    open(newunit=fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10,A10,A10,A10)") "#", "REAL", "NALLORB", "NKPOINT"
    write(fd, "(1X,L10,I10,I10,I10)") .true., nOrb, 1

    write (strForm, "(A,I0,A)") "(", nOrb, "ES24.15)"
    call unpackHS(square, sparse, iNeighbour, nNeighbourSK, iAtomStart, iPair, img2CentCell)
    call blockSymmetrizeHS(square, iAtomStart)
    write(fd, "(A1,A10,A10)") "#", "IKPOINT"
    write(fd, "(1X,I10,I10)") 1
    write(fd, "(A1,A)") "#", " MATRIX"
    write(fd, strForm) square
    close(fd)

  end subroutine writeSparseAsSquare_real


  !> Converts a sparse matrix to its square form and write it to a file.
  subroutine writeSparseAsSquare_cplx(env, fname, sparse, kPoints, iNeighbour, nNeighbourSK,&
      & iAtomStart, iPair, img2CentCell, iCellVec, cellVec)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Name of the file to write the matrix into.
    character(len=*), intent(in) :: fname

    !> Sparse matrix.
    real(dp), intent(in) :: sparse(:)

    !> List of k-points.
    real(dp), intent(in) :: kPoints(:,:)

    !> Neighbour list index.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours.
    integer, intent(in) :: nNeighbourSK(:)

    !> Offset array in the square matrix.
    integer, intent(in) :: iAtomStart(:)

    !> Pair indexing array.
    integer, intent(in) :: iPair(0:,:)

    !> Mapping of the atoms to the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Index of the cell translation vectors for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Cell translation vectors.
    real(dp), intent(in) :: cellVec(:,:)

    complex(dp), allocatable :: square(:,:)
    character(mc) :: strForm
    integer :: fd, nOrb, nKPoint
    integer :: iK

    if (withMpi) then
      call error("Writing of HS not working with MPI yet")
    end if

    nOrb = iAtomStart(size(nNeighbourSK) + 1) - 1
    nKPoint = size(kPoints, dim =2)

    allocate(square(nOrb, nOrb))
    open(newunit=fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10,A10,A10,A10)") "#", "REAL", "NALLORB", "NKPOINT"
    write(fd, "(1X,L10,I10,I10)") .false., nOrb, nKPoint

    write (strForm, "(A,I0,A)") "(", 2 * nOrb, "ES24.15)"
    do iK = 1, nKPoint
      call unpackHS(square, sparse, kPoints(:,iK), iNeighbour, nNeighbourSK, iCellVec, cellVec,&
          & iAtomStart, iPair, img2CentCell)
      call blockHermitianHS(square, iAtomStart)
      write(fd, "(A1,A10,A10)") "#", "IKPOINT"
      write(fd, "(1X,I10,I10)") iK
      write(fd, "(A1,A)") "#", " MATRIX"
      write(fd, strForm) square
    end do
    close(fd)

  end subroutine writeSparseAsSquare_cplx
  
  !> Converts a sparse matrix to its square form and write it to a file.
  subroutine writeSparseAsSquare_real_old(fname, sparse, iNeighbor, nNeighbor, iAtomStart, iPair, & !!DAR!!
      & img2CentCell) 

    !> Name of the file to write the matrix to.
    character(len=*), intent(in) :: fname

    !> Sparse matrix.
    real(dp), intent(in) :: sparse(:)

    !> Neighbor list index.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Number of neighbors.
    integer, intent(in) :: nNeighbor(:)

    !> Offset array in the square matrix.
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)

    !> Pair indexing array.
    integer, intent(in) :: img2CentCell(:)

    !> Mapping of the atoms to the central cell.
    real(dp), allocatable :: square(:,:)
    character(mc) :: strForm
    integer :: fd, nOrb

    nOrb = iAtomStart(size(nNeighbor) + 1) - 1

    allocate(square(nOrb, nOrb))
    fd = getFileId()
    open(fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10,A10,A10,A10)") "#", "REAL", "NALLORB", "NKPOINT"
    write(fd, "(1X,L10,I10,I10,I10)") .true., nOrb, 1

    write (strForm, "(A,I0,A)") "(", nOrb, "ES24.15)"
    call unpackHS(square, sparse, iNeighbor, nNeighbor, iAtomStart, iPair, &
        &img2CentCell)
    call blockSymmetrizeHS(square, iAtomStart)
    write(fd, "(A1,A10,A10)") "#", "IKPOINT"
    write(fd, "(1X,I10,I10)") 1
    write(fd, "(A1,A)") "#", " MATRIX"
    write(fd, strForm) square
    close(fd)

  end subroutine writeSparseAsSquare_real_old  


  !> Converts a sparse matrix to its square form and write it to a file.
  subroutine writeSparseAsSquare_cplx_old(fname, sparse, kPoints, iNeighbor, nNeighbor, iAtomStart, & !!DAR!!
      & iPair, img2CentCell, iCellVec, cellVec)

    !> Name of the file to write the matrix into.
    character(len=*), intent(in) :: fname

    !> Sparse matrix.
    real(dp), intent(in) :: sparse(:)

    !> List of k-points.
    real(dp), intent(in) :: kPoints(:,:)

    !> Neighbor list index.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Number of neighbors.
    integer, intent(in) :: nNeighbor(:)

    !> Offset array in the square matrix.
    integer, intent(in) :: iAtomStart(:)

    !> Pair indexing array.
    integer, intent(in) :: iPair(0:,:)

    !> Mapping of the atoms to the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Index of the cell translation vectors for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Cell translation vectors.
    real(dp), intent(in) :: cellVec(:,:)

    complex(dp), allocatable :: square(:,:)
    character(mc) :: strForm
    integer :: fd, nOrb, nKPoint
    integer :: iK

    nOrb = iAtomStart(size(nNeighbor) + 1) - 1
    nKPoint = size(kPoints, dim =2)

    allocate(square(nOrb, nOrb))
    fd = getFileId()
    open(fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10,A10,A10,A10)") "#", "REAL", "NALLORB", "NKPOINT"
    write(fd, "(1X,L10,I10,I10)") .false., nOrb, nKPoint

    write (strForm, "(A,I0,A)") "(", 2 * nOrb, "ES24.15)"
    do iK = 1, nKPoint
      call unpackHS(square, sparse, kPoints(:,iK), iNeighbor, nNeighbor, &
          &iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
      call blockHermitianHS(square, iAtomStart)
      write(fd, "(A1,A10,A10)") "#", "IKPOINT"
      write(fd, "(1X,I10,I10)") iK
      write(fd, "(A1,A)") "#", " MATRIX"
      write(fd, strForm) square
    end do
    close(fd)

  end subroutine writeSparseAsSquare_cplx_old


  !> Writes a sparse matrix to a file.
  subroutine writeSparse(fname, sparse, iNeighbour, nNeighbourSK, iAtomStart, iPair, img2CentCell,&
      & iCellVec, cellVec)

    !> Name of the file to write the matrix to.
    character(len=*), intent(in) :: fname

    !> Sparse matrix.
    real(dp), intent(in) :: sparse(:)

    !> Neighbour list index.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours.
    integer, intent(in) :: nNeighbourSK(:)

    !> Offset array in the square matrix.
    integer, intent(in) :: iAtomStart(:)

    !> Pair indexing array.
    integer, intent(in) :: iPair(0:,:)

    !> Mapping of the atoms to the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Index of the cell translation vectors for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Cell translation vectors.
    real(dp), intent(in) :: cellVec(:,:)

    integer :: fd, nAtom
    integer :: iAt1, iAt2, iAt2f, iNeigh, iOrig, nOrb1, nOrb2
    character(mc) :: strForm

    if (.not. tIoProc) then
      return
    end if

    nAtom = size(nNeighbourSK)

    open(newunit=fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10)") "#", "NATOM"
    write(fd, "(1X,I10)") nAtom
    write(fd, "(A1,A10,A10,A10)") "#", "IATOM", "NNEIGH", "NORB"
    do iAt1 = 1, nAtom
      write(fd, "(1X,I10,I10,I10)") iAt1, nNeighbourSK(iAt1) + 1, iAtomStart(iAt1+1)&
          & - iAtomStart(iAt1)
    end do

    do iAt1 = 1, nAtom
      nOrb1 = iAtomStart(iAt1+1) - iAtomStart(iAt1)
      do iNeigh = 0, nNeighbourSK(iAt1)
        iOrig = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = iAtomStart(iAt2f+1) - iAtomStart(iAt2f)
        write(strForm, "(A,I0,A)") "(", nOrb2, "ES24.15)"
        write(fd, "(A1,A10,A10,A10,3A10)") "#", "IATOM1", "INEIGH", "IATOM2F", "ICELL(1)",&
            & "ICELL(2)", "ICELL(3)"
        write(fd, "(1X,I10,I10,I10,3I10)") iAt1, iNeigh, iAt2f, int(cellVec(:,iCellVec(iAt2)))
        write(fd, "(A1,A)") "#", " MATRIX"
        write(fd, strForm) sparse(iOrig:iOrig+nOrb1*nOrb2-1)
      end do
    end do
    close(fd)

  end subroutine writeSparse

end module dftbp_formatout
