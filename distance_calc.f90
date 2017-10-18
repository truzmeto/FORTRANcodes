!## Talant Ruzmetov 05/03/2013
!## This code calculates ceter of mass distance between pKID and KIX
!## and alignment of helix_B end to end vector with respect to native
!## orientation. It reads movie file from str input.

program dist_calculation
  implicit none
  integer,parameter:: nmax = 200
  real:: x, y, z, x1_ave, y1_ave, z1_ave, x2_ave, y2_ave, z2_ave
  integer:: i, j, l, n, frame_n, N_beads, ires
  integer:: res_num1, res_num2,res_num,res_num_b,res_num_a, n_link
  character(LEN=72):: char72
  character(LEN=4)::nameofatom
  character(LEN=3):: resname
  real,dimension(nmax,3)::rca, rca2, rca1
  real::dist,x1,x2,y1,y2,z1,z2,d_cm,rcm,rcm1,rcm2
  real::x3,y3,z3,d_cm_b,d_cm_a,x3_ave,y3_ave,z3_ave,x4,y4,z4,x4_ave,y4_ave,z4_ave
  character (LEN=3),dimension(nmax)::Lsequence
  character(LEN=32)::arg1

  real:: x_nat95,y_nat95,z_nat95,x_nat109,y_nat109,z_nat109,bx_nat,by_nat,bz_nat,bx,by,bz,theta_b
  real:: x_nat82,y_nat82,z_nat82,x_nat92,y_nat92,z_nat92,ax_nat,ay_nat,az_nat,ax,ay,az,theta_a
  
  !###############################-- some parameters--##############################################
  res_num1 = 81                 !|> number of residues for unit1 (KIX)
  res_num2 = 28                 !|> number of residues for unit2 (pKID)
  res_num_b = 15                !|> number of residues for subunit of unit2 (pKID helix_B)
  res_num_a = 11                !|> number of residues for subunit of unit2 (pKID helix_A)
  n_link = 2                    !|> num res for linker
  N_beads = res_num1 + res_num2 !|> Total number of residues 

  
  !####################--argument passing from term--################################################
  do j = 1, iargc()
     call getarg(j, arg1)
     read(arg1,*) frame_n  !|> total number of frames for trajectory file
  enddo

  
  !#####################-- Native head and tail coordinates of helices A and B --####################    

!# ATOM     95  CA  PRO B  14       0.527   9.018  -1.924
  x_nat95 = 0.527  !|> helix-B head
  y_nat95 = 9.018  
  z_nat95 = -1.924

  
!# ATOM    109  CA  PRO B  28     -14.186  -4.007   7.980  1.00  1.00
  x_nat109 = -14.186  !|> helix-B tail
  y_nat109 = -4.007 
  z_nat109 = 7.980 

!# ATOM     82  CA  THR B   1      -7.965   2.639 -19.356  1.00  1.00
  x_nat82 = -7.965   !|> helix-A head
  y_nat82 =  2.639  
  z_nat82 = -19.356 

!# ATOM     92  CA  SER B  11      -0.113   6.452  -7.244  1.00  1.00
  x_nat92 = -0.113  !|> helix-A tail
  y_nat92 =  6.452
  z_nat92 = -7.244 

  
  !####################-- vector A=ax*i+ay*j+az*k--#################################
  !####################-- vector B=bx*i+by*j+bz*k--#################################
  bx_nat = -x_nat95 + x_nat109
  by_nat = -y_nat95 + y_nat109
  bz_nat = -z_nat95 + z_nat109

  ax_nat = -x_nat82 + x_nat92
  ay_nat = -y_nat82 + y_nat92
  az_nat = -z_nat82 + z_nat92      
  
    
  do l = 1, frame_n
     ires = 0
     Lines: do
        read(5,'(a72)',end=999) char72
        if(char72(1:3)=='END')then 
           exit Lines
        endif
        
        if(char72(1:4) /= 'ATOM') cycle Lines
        read(char72,1000) nameofatom,resname,x,y,z
        if(nameofatom .eq. ' CA ') then
           ires = ires + 1
           rca(ires,1) = x
           rca(ires,2) = y
           rca(ires,3) = z
           Lsequence(ires) = resname
        endif
1000    format(12x,a4,1x,a3,10x,3f8.3)
        
     end do Lines

     x1_ave = 0.0
     y1_ave = 0.0
     z1_ave = 0.0
     x2_ave = 0.0
     y2_ave = 0.0
     z2_ave = 0.0
     x3_ave = 0.0
     y3_ave = 0.0
     z3_ave = 0.0
     x4_ave = 0.0
     y4_ave = 0.0
     z4_ave = 0.0

     
     do i = 1, N_beads
        
        !######## Since KIX is frozen we only need to do this once #############
        !# calc. KIX COM
        if(l == 1) then
           if(i.le.res_num1) then
              x1 = rca(i,1) 
              y1 = rca(i,2)
              z1 = rca(i,3) 
              
              x1_ave = x1_ave + x1
              y1_ave = y1_ave + y1
              z1_ave = z1_ave + z1
           endif
        endif
        
        !# calc. pKID COM 
        if(i.gt.res_num1) then
           x2 = rca(i,1) 
           y2 = rca(i,2)
           z2 = rca(i,3) 
           
           x2_ave = x2_ave + x2
           y2_ave = y2_ave + y2
           z2_ave = z2_ave + z2
        endif

        !# calc. helix-B COM
        if(i.gt.res_num1 + res_num_a + n_link) then         
           x3_ave = x3_ave + rca(i,1)
           y3_ave = y3_ave + rca(i,2)
           z3_ave = z3_ave + rca(i,3)
        endif

        !# calc. helix-A COM
        if(i.gt.res_num1.and.i.le.res_num1 + res_num_a) then         
           x4_ave = x4_ave + rca(i,1)
           y4_ave = y4_ave + rca(i,2)
           z4_ave = z4_ave + rca(i,3)
        endif
         
         
        !#########################--angle calculation--#############################    
        bx = -rca(95,1) + rca(109,1)
        by = -rca(95,2) + rca(109,2)
        bz = -rca(95,3) + rca(109,3)
        
        ax = -rca(82,1) + rca(92,1)
        ay = -rca(82,2) + rca(92,2)
        az = -rca(82,3) + rca(92,3)
        
        !#########################--vector A=ax*i+ay*j+az*k--#######################
        !#########################--vector B=bx*i+by*j+bz*k--#######################
        theta_a = (ax*ax_nat+ay*ay_nat+az*az_nat)/(sqrt((ax*ax+ay*ay+az*az)*(ax_nat*ax_nat+ay_nat*ay_nat+az_nat*az_nat)))
        theta_b = (bx*bx_nat+by*by_nat+bz*bz_nat)/(sqrt((bx*bx+by*by+bz*bz)*(bx_nat*bx_nat+by_nat*by_nat+bz_nat*bz_nat)))
                
     enddo
     
     x1 = x1_ave/real(res_num1)
     y1 = y1_ave/real(res_num1)
     z1 = z1_ave/real(res_num1)
  
     x2 = x2_ave/real(res_num2)
     y2 = y2_ave/real(res_num2)
     z2 = z2_ave/real(res_num2)
     
     x3 = x3_ave/real(res_num_b)
     y3 = y3_ave/real(res_num_b)
     z3 = z3_ave/real(res_num_b)

     x4 = x4_ave/real(res_num_a)
     y4 = y4_ave/real(res_num_a)
     z4 = z4_ave/real(res_num_a)
          

     !####################--D_cm calculation--################################     
     d_cm = sqrt((x1-x2)**2.0 + (y1-y2)**2.0 + (z1-z2)**2.0)
     d_cm_b = sqrt((x1-x3)**2.0 + (y1-y3)**2.0 + (z1-z3)**2.0)
     d_cm_a = sqrt((x1-x4)**2.0 + (y1-y4)**2.0 + (z1-z4)**2.0)
     
     write(6,'(6f7.2)') d_cm, d_cm_a, d_cm_b, theta_a, theta_b 
  end do
999 continue 
end program dist_calculation
