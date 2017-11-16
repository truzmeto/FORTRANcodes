!###################################################################
!                    Talant Ruzmetov  11/21/2015 
! This code reads pdb trajectory and calculates mean square displacement 
! of a protein vs time. Later, output of the code can be used to evaluate
! diffusion constant of the macromolecule.
!###################################################################
program msdtime
  implicit none
  integer,parameter:: nmax = 60 !| shoud be more than res_num 
  real:: x,y,z,x_ave,y_ave,z_ave,dr2, mean2, error
  integer:: i,j,l,n,frame_n, idelta, keep_frames,step_save
  integer:: res_num
  character(LEN=72):: char72,char55
  character(LEN=4)::nameofatom
  character(LEN=3):: resname
  real, dimension(nmax,3)::rca
  real:: r_cm,mean
  character (LEN=3),dimension(nmax)::Lsequence
  character(LEN=80):: char500
  character(LEN=80):: char200
  character(LEN=7)::nameid
  integer,dimension(1000000):: counter
  character(LEN=32)::arg1
  real,dimension(1000000):: x_cm,y_cm,z_cm, msd, msd2
  
!########################--parameters which have to be given manually--############################
  res_num = 56      !## 1SRL
!##################################################################################################

step_save = 100  
!####################--argument passing from term--################################################
  do i = 1, iargc()
     call getarg(i, arg1)
     read(arg1,*) frame_n
  enddo
  !##################################################################################################
  
  !!  write(*,*) "i", "mean", "counter(i)", "mean2", "error"
  do l=1,frame_n
     
     i=0
     Lines: do
        read(5,'(a72)',end=10) char72
        
        if(char72(1:3)=='END')then 
           exit Lines
        endif
        
        if(char72(1:4) /= 'ATOM') cycle Lines
        
        read(char72,1000) nameofatom,resname,x,y,z
        if(nameofatom .eq. ' CA ') then
           i=i+1
           rca(i,1) = x
           rca(i,2) = y
           rca(i,3) = z
           Lsequence(i) = resname
        endif
1000    format(12x,a4,1x,a3,10x,3f8.3)
        
     end do Lines
     n=i
     
     x_ave = 0.0
     y_ave = 0.0
     z_ave = 0.0
     
     do i=1,n
        x = rca(i,1) 
        y = rca(i,2)
        z = rca(i,3) 
        
        x_ave = x_ave + x
        y_ave = y_ave + y
        z_ave = z_ave + z
     enddo
     
     x_cm(l) = x_ave / real(n)
     y_cm(l) = y_ave / real(n)
     z_cm(l) = z_ave / real(n)
     
     ! write(6,'(i8,2f10.2)') l, r_cm, msd/real(l)
     msd(l) = 0.0
     counter(l) = 0
  end do
10 continue 
  
  do i = 1,frame_n-1
     do j = i+1, frame_n
        idelta = j - i  
        dr2= (x_cm(j)-x_cm(i))**2 + (y_cm(j)-y_cm(i))**2 + (z_cm(j)-z_cm(i))**2 
        msd(idelta) = msd(idelta) + dr2
        msd2(idelta) = msd2(idelta) + dr2**2
        counter(idelta) = counter(idelta) + 1 
     enddo
  enddo

  keep_frames = frame_n / 5
  do i = 1, keep_frames
     mean = msd(i)/counter(i)
     mean2 = msd2(i)/counter(i)
     error = sqrt( ( mean2 - mean**2)/counter(i) )
     write(*,*) i*step_save, mean !, counter(i), mean2, error 
  enddo
  
end program msdtime
