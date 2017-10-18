
!###########################-- Talant Ruzmetov 10/06/2017 --#################################
!##      This code calculates all intermolecular contacts between 2 units.                  #
!##      Spesifically, it stores pairwise total intermolecular contacts for each frame      #
!##      in the form of matrix, where u(i,j)=1 for contact and 0 otherwise,                 #
!##      so we can use its output in order to calculate amount of total distinct            #
!##      contacts per attempt for kinetics                                                  #
!!###########################################################################################

program ContactMatrix 

  implicit none

  integer,parameter:: nmax = 200
  real:: x, y, z, dx, dy, dz
  real:: x1_ave, y1_ave, z1_ave, x2_ave, y2_ave, z2_ave
  real:: x1,x2,y1,y2,z1,z2
  real:: gen_dist_all
  integer:: ii, jj, count, k, m, step_size
  integer:: tot_contacts
  integer:: i, j, l, n, iframe, frame_n, N_beads, ires
  integer:: res_num1, res_num2
  character(LEN=72):: char72
  character(LEN=80):: char500
  character(LEN=4):: nameofatom
  character(LEN=3):: resname
  character(LEN=7):: nameid
  real,dimension(nmax,3):: rca
  integer,dimension(3000,2):: u, u_all
  integer,dimension(120,120):: u_new
  real:: dist, d_cm, go_nat_dist
  character(LEN=3),dimension(nmax):: Lsequence
  character(LEN=32):: arg1
  
!###############################-- some parameters--##############################################
  res_num1 = 81                  !|> number of residues for unit1 (KIX)
  res_num2 = 28                  !|> number of residues for unit2 (pKID)
  N_beads = res_num1 + res_num2  !|> Total number of residues 
  gen_dist_all = 8.3             !|> generic non-native + native contact distance 
  step_size = 100                !|> saved step frequency
  
!####################--argument passing from term--###############################################
  do j = 1, iargc()
     call getarg(j, arg1)
     read(arg1,*) frame_n        !|> total number of frames for trajectory file
  enddo
  
!#######################-- generating all intermolecular contact indecies --######################
  l = 0
  do i = 1, res_num1
     do j = res_num1 + 1, N_beads
        l = l + 1
        u_all(l,1) = i
        u_all(l,2) = j
        u_new(i,j) = 0
     enddo
  enddo
  tot_contacts = l
    
!########################-- looping over movie frames to read coordinates --########################3  
  do iframe = 1, frame_n  !# this is faster than rewinding!
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
     


!!############ dist calculation ####################################################     
     x1_ave = 0.0
     y1_ave = 0.0
     z1_ave = 0.0
     x2_ave = 0.0
     y2_ave = 0.0
     z2_ave = 0.0
    
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
        
     enddo

     x1 = x1_ave/real(res_num1)
     y1 = y1_ave/real(res_num1)
     z1 = z1_ave/real(res_num1)
  
     x2 = x2_ave/real(res_num2)
     y2 = y2_ave/real(res_num2)
     z2 = z2_ave/real(res_num2)

     d_cm = sqrt((x1-x2)**2.0 + (y1-y2)**2.0 + (z1-z2)**2.0)
        
!!#####################################################################################


     
!!############  contact calc ########################     
     count = 0
     do m = 1, tot_contacts
        ii = u_all(m,1)
        jj = u_all(m,2)
        dx = rca(ii,1) - rca(jj,1)
        dy = rca(ii,2) - rca(jj,2)
        dz = rca(ii,3) - rca(jj,3)
        
        dist = dx*dx + dy*dy + dz*dz
        dist = sqrt(dist)
        
        if(dist.le.1.200 * gen_dist_all) then
           u_new(ii,jj) = 1
           count = count + 1                  
        else
           u_new(ii,jj) = 0
        endif

!       write(*,*) dist
     end do
     q_all = count !/ real(tot_contacts)
!!###########################################################################################
     

     
!!########### here I sort out all distinct contacts in the encounter per attempt for unsuccessful trials
     state = 1
     t_cap = 0
     
     dt = dt + step
     if(state.eq.1.and.q_all.ge.q_thresh) then
     !if(state.eq.1.and.r.lt.r_lim) then  ! just to define cap differently
        state = 2
        t_cap = t_cap + dt
        n_cap = n_cap + 1
        dt = 0
     elseif(state.eq.2.and.q.ge.q_b_min) then
        state = 3
        t_evol = t_evol + dt
        n_evol = n_evol + 1
        dt = 0
     elseif(state.eq.2.and.r.gt.r_lim) then
        state = 1
        t_esc = t_esc + dt
        n_esc = n_esc + 1
        dt = 0
     endif
!!#######################################################################################################
     
     
!!################## output results as matrix ##############################################
     write(*,'(a13,1x,i8)') ">>>> iframe =", iframe*step_size
     write(*,'(a25,1x,i4)') ">>>> tot contact amount =", count
     do i = 1, res_num1
        write(*,'(<res_num2>(i2,1x))')(u_new(i,j), j = res_num1 + 1, N_beads)
     enddo
     write(*,'(a50)') " "
     
!!##########################################################################################
     
  end do
999 continue 
end program ContactMatrix
