program current

  integer, parameter :: dp=8 
  integer :: natoms, ngroups, natom_device    
  integer, dimension(:,:), allocatable :: neig     
  integer, dimension(:), allocatable :: nn     
  real(8), dimension(:,:), allocatable :: Inm    
  real(8), dimension(:,:), allocatable :: coord
  integer :: m, tmp, i, j, n, dir 
  real(8) :: z, Iz, lstep, dr(3)
  character(64) :: arg, filename


  if (iargc().lt.6) then 
     write(*,*) "current 'lcurrents.dat' natoms natom_device dir ngroups lstep"
     stop   
  endif

  CALL getarg(1, filename)
  
  CALL getarg(2, arg)
  read(arg,*) natoms
  
  CALL getarg(3, arg)
  read(arg,*) natom_device
  
  CALL getarg(4, arg)
  read(arg,*) dir
  
  CALL getarg(5, arg)
  read(arg,*) ngroups

  CALL getarg(6, arg)
  read(arg,*) lstep 

  !print*,natoms,dir,ngroups,lstep

  n = 40

  allocate(neig(natoms,n))
  allocate(nn(natoms))
  allocate(Inm(natoms,natoms))
  allocate(coord(3,natoms))

  open(105,file=trim(filename))
  do m=1, natoms    
     read(105,*) tmp, coord(1:3,m), nn(m), (neig(m,i), Inm(m,neig(m,i)), i=1,nn(m)) 
     !print*,nn(m)
  enddo
  close(105)

  !COMPUTE total current layer by layer
  z = minval(coord(dir,1:natom_device))+0.1_dp

  do j = 1, ngroups

    !print*,z
    Iz = 0.0_dp
    do m = 1, natom_device
      if (coord(dir,m).lt.z) cycle
      do i = 1, nn(m) 
        if (coord(dir,neig(m,i)).gt.z) cycle
        dr(:)=coord(:,neig(m,i))-coord(:,m)
        !print*,'<z >z',m,neig(m,i)
        Iz=Iz+Inm(m,neig(m,i)) !*abs(dr(dir))/sqrt(dot_product(dr,dr))
      enddo
    enddo

    print*,z*0.529177_dp,Iz
    z = z + lstep/0.529177_dp
  enddo

  deallocate(neig,nn,Inm,coord)

end program current      
