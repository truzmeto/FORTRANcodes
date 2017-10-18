
!###########################-- Talant Ruzmetov 10/06/2017 --#################################
!##      This code calculates all intermolecular contacts between 2 units.                 ##
!##      Specifically it calculates amount of total distinct contacts made per             ##
!##      attempt from kinetics binding trajectory                                          ##
!############################################################################################

program NumberOfDistinctContactsPerAttempt 

  implicit none

  integer,parameter:: nmax = 200
  real:: x, y, z, dx, dy, dz
  real:: gen_dist_all, r_lim
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
  integer,dimension(120,120):: u_new, collapse_cont
  real:: dist, d_cm, go_nat_dist
  real,dimension(3):: Rcm_pro, Rcm_lig, dr
  integer,dimension(500):: n_dist
  integer state,n_cap,n_esc,n_evol, step, q_thresh,q, tot_dist

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

  r_lim = 40.000 
  n_esc = 0
  n_cap = 0
  n_evol = 0
  step = 100
  state = 1
  q_thresh = 4

  n_dist(1:200) = 0

  
  l = 0
  do i = 1, res_num1
     do j = res_num1 + 1, N_beads
        l = l + 1
        u_all(l,1) = i
        u_all(l,2) = j
        u_new(i,j) = 0
        collapse_cont(i,j) = 0
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

!######---------------- Dcm calculation------------------------------------#################
     Rcm_pro(1:3) = 0.00
     Rcm_lig(1:3) = 0.00
     do i = 1, N_beads
        !######## Since KIX is frozen we only need to do this once #############
        !# calc. KIX COM
        if(l == 1) then
           if(i.le.res_num1) then
              Rcm_pro(1:3) = Rcm_pro(1:3) + rca(i,1:3)
           endif
        endif
        
        !######### calc. pKID COM 
        if(i.gt.res_num1) then
           Rcm_lig(1:3) = Rcm_lig(1:3) + rca(i,1:3)
        endif
     enddo

     Rcm_pro(1:3) = Rcm_pro(1:3) / real(res_num1)
     Rcm_lig(1:3) = Rcm_lig(1:3) / real(res_num2)
     dr(1:3) = Rcm_pro(1:3) - Rcm_lig(1:3)
     
     !####################--D_cm between two units calculation--################################     
     d_cm = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
!###############################################################################################

     
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
     end do
   
     do i = 1, res_num1
        do j = res_num1 + 1, N_beads
           collapse_cont(i,j) = collapse_cont(i,j) + u_new(i,j)
        enddo
     enddo

     !!##---------------- here we calculate # of distinct contacts per attempt-----------------------
     q = count
     if(state.eq.1.and.q.ge.q_thresh) then
        state = 2
        n_cap = n_cap + 1
     elseif(state.eq.2.and.d_cm.gt.r_lim) then
        state = 1
        n_esc = n_esc + 1
      
        do i = 1, res_num1
           do j = res_num1 + 1, N_beads
              if(collapse_cont(i,j).ne.0)  n_dist(n_esc) = n_dist(n_esc) + 1
              collapse_cont(i,j) = 0
           enddo
        enddo
      endif
!!##########################################################################################
     
   end do

   tot_dist = 0
   do i = 1, n_esc
      tot_dist = tot_dist + n_dist(i)
   enddo
   
   write(*,*) n_esc, tot_dist
   
   
999 continue
   
   !write(*,'(a26,1x,i3,2x,i8)') ">>>> #distinct contacts =", i, n_dist(i)
   
   ! do i = 1, res_num1
   !    write(*,'(<res_num2>(i2,1x))')(collapse_cont(i,j), j = res_num1 + 1, N_beads)
   ! enddo
 end program  NumberOfDistinctContactsPerAttempt
