program contactmap
  implicit none
  integer,parameter::nmax = 400
  real:: x,y,z,go_nat,dx,dy,dz
  integer:: i,j,ii,m,frame_n,fn,l,mm,ncon11,ncon12,ncon22,jj,i1,ll,ncon
  integer:: res_num1, res_num2, res_num
  character(LEN=72):: char72,char55
  character(LEN=4):: nameofatom
  character(LEN=3):: resname
  integer::nca,icon,imp1,imp2,imp1un,imp2un,iunit1,iunit2
  real,dimension(nmax,3)::rca, rca2, rca1
  integer,dimension(nmax):: q_norm,co,q_local

  real::dist
  integer,dimension(nmax,2)::u
  character (LEN=3),dimension(nmax)::Lsequence
  real,dimension(nmax)::dist_nat,dist_nat11,dist_nat12,dist_nat22
  character(LEN=80):: char500
  character(LEN=7)::nameid
  character(LEN=32)::arg1
   
  
!########## manually given params ###################################
  ncon11 = 271
  ncon12 = 77
  ncon22 = 51
  ncon = ncon11 + ncon12 + ncon22
  res_num1 = 81
  res_num2 = 28

  
!!###################--- Argument passing ---####################################
  do i = 1, iargc()
     call getarg(i, arg1)
     read(arg1,*) frame_n
  enddo
  
  
!!############## Reading native info file ########################################################  
  ll = 0
  open(unit = 71,file = '1KDX.ninfo', status = 'unknown')     
  ninfo_lines:do
     read(71,'(a80)', END = 550) char500          
     if(char500(1:7) /= 'contact') cycle ninfo_lines
     read(char500,'(a7,4x,i3,5x,i2,5x,i2,4x,i3,4x,i3,5x,i2,5x,i2,5x,f10.5)') nameid,icon,iunit1, &
          iunit2,imp1,imp2,imp1un,imp2un,go_nat
     
     ll = ll + 1
     u(ll,1) = imp1
     u(ll,2) = imp2
     dist_nat(ll) = go_nat
     
  end do ninfo_lines
550 continue
!!##################################################################################################  


  
!!## parsing traj file and extracting coords ##################################3  
  do l = 1, frame_n
     i = 0
     Lines: do
        read(5,'(a72)',end=999) char72
        if(char72(1:3)=='END')then 
           exit Lines
        endif
        if(char72(1:4) /= 'ATOM') cycle Lines
        read(char72,1000) nameofatom,resname,x,y,z
        if(nameofatom .eq. ' CA ') then
           i = i + 1
           rca(i,1) = x
           rca(i,2) = y
           rca(i,3) = z
           Lsequence(i) = resname
        endif
1000    format(12x,a4,1x,a3,10x,3f8.3)
     end do Lines


     
!########### calculate distance from frame and check with go_nat ##############################
     
     do i = 1, res_num2
        q_local(i) = 0
     enddo
     
     do i1 = ncon11 + 1, ncon11 + ncon12
        
        ii = u(i1,1)
        jj = u(i1,2)
        dx = rca(ii,1) - rca(jj,1)
        dy = rca(ii,2) - rca(jj,2)
        dz = rca(ii,3) - rca(jj,3)
        
        dist = dx*dx + dy*dy + dz*dz
        dist = sqrt(dist)
               
        if(dist <= 1.500 + dist_nat(i1))then
           
           do i = 1, res_num2
              if(jj - res_num1.eq.i) then
                 q_local(i) = q_local(i) + 1
              endif
           enddo
           
        endif
        
     end do
     write(*,'(<res_num2>(i2,1x))')(q_local(i), i = 1, res_num2)
!     write(*,'(28(i2,1x))')(q_local(i), i = 1, res_num2)

  end do
999 continue 

end program contactmap
