!#########----- Program written by ------###############
!#########----- Talant Ruzmetov----------###############
!#########----- July, 2017---------------###############

!# This code reads movie(velo_pdb) file where, x,y,z is replaced with vx, vy, vz and
!# extracts relative velocity of each residue with respect to CM. Then, it outputs
!# V_CM alongside with all residual velocities

program VelocityCalculation
  implicit none
  integer,parameter::nmax = 400
  real:: vx,vy,vz,v_cm_x,v_cm_y,v_cm_z,vel_cm 
  real:: vr_x,vr_y,vr_z, vr
  integer:: i,j,iframe,m,frame_n
  integer:: res_num1, res_num2, res_num
  character(LEN=72):: char72
  character(LEN=4):: nameofatom
  character(LEN=3):: resname
  real,dimension(nmax,3):: rca
  real,dimension(nmax):: v_loc

!  integer,dimension(nmax,2)::u
  character (LEN=3),dimension(nmax)::Lsequence
  character(LEN=32)::arg1
   
  
!########## manually given params ###################################
  res_num1 = 81
  res_num2 = 28

  
!!###################--- Argument passing ---####################################
  do i = 1, iargc()
     call getarg(i, arg1)
     read(arg1,*) frame_n
  enddo
  
  
  !## parsing traj file and extracting velocities per residue
  
  do iframe = 1, frame_n
     i = 0
     Lines: do
        read(5,'(a72)',end=999) char72
        if(char72(1:3)=='END')then 
           exit Lines
        endif
        if(char72(1:4) /= 'ATOM') cycle Lines
        read(char72,1000) nameofatom,resname,vx,vy,vz
        if(nameofatom .eq. '  CA') then
           i = i + 1
           rca(i,1) = vx
           rca(i,2) = vy
           rca(i,3) = vz
           Lsequence(i) = resname

        endif

1000    format(12x,a4,1x,a3,10x,3f8.3)
     end do Lines
     
!############### calculate CM velocity for 2nd unit   
     v_cm_x = 0.0
     v_cm_y = 0.0
     v_cm_z = 0.0
    ! v_cm(1:3) = 0.0
     do j = res_num1 + 1, res_num1 + res_num2
        v_cm_x = v_cm_x + rca(j,1)
        v_cm_y = v_cm_y + rca(j,2)
        v_cm_z = v_cm_z + rca(j,3)
     !   v_cm(1:3) = v_cm(1:3) + rca(j,1:3)
     enddo
     v_cm_x = v_cm_x / real(res_num2)
     v_cm_y = v_cm_y / real(res_num2)
     v_cm_z = v_cm_z / real(res_num2)
     vel_cm = v_cm_x*v_cm_x + v_cm_y*v_cm_y + v_cm_z*v_cm_z  
     vel_cm = sqrt(vel_cm)
     ! v_cm(1:3) = v_cm(1:3) / real(res_num2)
     
     
     do j = res_num1 + 1, res_num1 + res_num2
        
        vr_x = rca(j,1) - v_cm_x
        vr_y = rca(j,2) - v_cm_y
        vr_z = rca(j,3) - v_cm_z
        
        vr = vr_x*vr_x + vr_y*vr_y + vr_z*vr_z
        v_loc(j) = sqrt(vr)
       
     end do
     v_loc(res_num1) = vel_cm
     write(*,'(<res_num2+1>f7.3,2x)')( v_loc(j), j = res_num1 , res_num1 + res_num2)
!    write(*,'(28(i2,1x))')(v_loc(j), j = 1, res_num2) !# use this with gfortran compiler

  end do
999 continue 

end program VelocityCalculation
