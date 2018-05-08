program makecube


integer :: i,j,k,m, nx,ny,nz, p1,p2, narg, ln, err
real(8), dimension(:,:,:), allocatable :: phi3d
real(8), dimension(:), allocatable :: x,y,z
real(8), dimension(3) :: or, dl
real(8), parameter :: au=0.529177
character(256) :: filebox,filex,filey,filez,filename
character(80) :: buffer
character(1) :: dir

narg=command_argument_count()

if (narg.lt.1 .or. (narg.gt.4 .and. narg.lt.8) ) then
 write(*,*) 'usage:'
 write(*,*) 'makecube pot_file z nx ny [boxfile xfile yfile zfile]'
 stop 
endif
 
call get_command_argument(1,filename,ln,err)
call get_command_argument(2,dir,ln,err)
call get_command_argument(3,buffer,ln,err)
read(buffer,*) p1
call get_command_argument(4,buffer,ln,err)
read(buffer,*) p2

if (narg.lt.5) then
  filebox="box3d.dat"
  filex="Xvector.dat"
  filey="Yvector.dat"
  filez="Zvector.dat"  
else
  call get_command_argument(5,filebox,ln,err)
  call get_command_argument(6,filex,ln,err)
  call get_command_argument(7,filey,ln,err)
  call get_command_argument(8,filez,ln,err)
endif 


open(110,file=filebox)
read(110,*) nx,ny,nz
close(110)

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))
allocate(phi3d(nx,ny,nz))

open(110,file=filex)
read(110,*) x
close(110)

open(110,file=filey)
read(110,*) y
close(110)

open(110,file=filez)
read(110,*) z
close(110)


open(110,file=filename)
read(110,'(a80)') buffer
read(buffer,'(3i)', iostat=err) i,j,k

if(err.ne.0) rewind(110)

do i=1,nx
 do j=1,ny 
  do k=1,nz
    read(110,*) phi3d(i,j,k)
  enddo
 enddo
enddo
close(110)

open(1000,file="pot_line.dat")
select case(dir)
 case("x")
   do k = 1, nx
     write(1000,*) x(k), phi3d(k,p1,p2) 
   enddo
 case("y")
   do k = 1, ny
     write(1000,*) y(k), phi3d(p1,k,p2) 
   enddo
 case("z")
   do k = 1, nz
     write(1000,*) z(k), phi3d(p1,p2,k) 
   enddo
end select



close(1000)

deallocate(x,y,z)
deallocate(phi3d)

end program makecube
