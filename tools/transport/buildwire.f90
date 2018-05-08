program buildwire

      integer :: i,d(3),pl_atm,tot_atm,atm_cont,num_pl,err,k
      character(30) :: gen_file
      character(2) :: period
      character(70) :: atm_spec
      character(100) :: cell_centre
      real(8) :: cell(3,3)
      integer, ALLOCATABLE, DIMENSION (:) :: n_atm,typ_atm
      real(8), ALLOCATABLE, DIMENSION (:) :: X,Y,Z      
      
      write(*,*) 'Insert PL .gen file name: '
      read(*,*) gen_file
      
      write(*,*) 'Insert direction (e.g. 0 0 1)'
      read(*,*) d(1:3)
    
      open(30,file=gen_file)
      
      read(30,*) pl_atm, period
      read(30,'(A70)') atm_spec

      ALLOCATE(X(pl_atm),stat=err)
      IF (err /= 0) STOP 'no space for allocation (X)'

      ALLOCATE(Y(pl_atm),stat=err)
      IF (err /= 0) STOP 'no space for allocation (Y)'

      ALLOCATE(Z(pl_atm),stat=err)
      IF (err /= 0) STOP 'no space for allocation (Z)'
      
      ALLOCATE(n_atm(pl_atm),stat=err)
      IF (err /= 0) STOP 'no space for allocation (n_atm)'
            
      ALLOCATE(typ_atm(pl_atm),stat=err)
      IF (err /= 0) STOP 'no space for allocation (typ_atm)'
                  
      
      do i = 1,pl_atm
         read(30,*) n_atm(i),typ_atm(i),X(i),Y(i),Z(i)
      end do
      
      read(30,'(A)') cell_centre
      do i = 1,3
         read(30,*) cell(i,1),cell(i,2),cell(i,3)
      end do

      close(30)

      write(*,*) 'Insert number of PLs in channel: '
      read(*,*) num_pl
                 
      tot_atm=pl_atm * (num_pl + 4)     


      open(30,file='Ordered_'//gen_file)

      write(30,'(I5,A4)') tot_atm,period
      write(30,'(A70)') atm_spec
      do k = 1,num_pl
         do i = 1,pl_atm
            write(30,'(I5,I5,F20.12,F20.12,F20.12)') n_atm(i),typ_atm(i), &
                                   X(i)+(k-1)*cell(1,1)*d(1),& 
                                   Y(i)+(k-1)*cell(2,2)*d(2), &
                                   Z(i)+(k-1)*cell(3,3)*d(3) 
         end do
      enddo

      ! build I contact
      do k = num_pl+1,num_pl+2
         do i = 1,pl_atm
            write(30,'(I5,I5,F20.12,F20.12,F20.12)') n_atm(i),typ_atm(i), &
                                   X(i)+(k-1)*cell(1,1)*d(1),& 
                                   Y(i)+(k-1)*cell(2,2)*d(2), &
                                   Z(i)+(k-1)*cell(3,3)*d(3) 
         end do
      enddo

      ! build II contact
      do k = 0,-1,-1
         do i = 1,pl_atm
            write(30,'(I5,I5,F20.12,F20.12,F20.12)') n_atm(i),typ_atm(i), &
                                   X(i)+(k-1)*cell(1,1)*d(1),& 
                                   Y(i)+(k-1)*cell(2,2)*d(2), &
                                   Z(i)+(k-1)*cell(3,3)*d(3) 
         end do
      enddo


      write(30,'(A)') cell_centre
      do i = 1,3
          if(d(1).eq.1) write(30,*) cell(i,1)*(num_pl+4)*d(1),cell(i,2),cell(i,3)
          if(d(2).eq.1) write(30,*) cell(i,1),cell(i,2)*(num_pl+4)*d(2),cell(i,3)
          if(d(3).eq.1) write(30,*) cell(i,1),cell(i,2),cell(i,3)*(num_pl+4)*d(3)
      end do

      close(30)

      write(*,*) 'structure built'
      write(*,*) '*iatc='
      write(*,*) pl_atm*num_pl+1, pl_atm*(num_pl+2), '0'
      write(*,*) pl_atm*(num_pl+2)+1, pl_atm*(num_pl+4), '0'

      write(*,*) '*PLs='
      do i=1,num_pl-1
         write(*,FMT='(i6)',ADVANCE='NO') (i-1)*pl_atm+1
      enddo
     
        write(*,FMT='(i6,a)') (num_pl-1)*pl_atm+1,';' 



end program
      
      
        
