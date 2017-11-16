!###################################################################
!                    T. Ruzmetov  written on 10/10/2017
! This code reads movie + native info file and calculates
! reaction coordinates(q_int, d_cm). Then conditioning on those reaction coordinates
! it splits binding trajectory into sub parts(attempts and final binding). At the end
! it outputs a contact_vector for all attempts sorted by n_esc in the first column.
! contact_vector is a vector of size=n_tot_int_contacts, composed of zero's and ones, where
! 1 means that specific contact is made and vise versa. 
!###################################################################

program ExtractKinSteps
  implicit none
  integer,parameter::nmax = 400
  real:: x, y, z, dx, dy, dz
  real:: go_nat_dist, d_cm, r_lim
  integer:: count, q_thresh
  integer:: q, step_save
  integer:: i, j, ii, l, jj, ll, id, n_cap, n_esc, state
  integer:: frame_n, ncon11, ncon12, ncon22, ncon
  integer:: res_num1, res_num2, res_num, N_beads
  real,dimension(3):: Rcm_pro, Rcm_lig, dr

  character(LEN=72):: char72,char55
  character(LEN=4)::nameofatom
  character(LEN=3):: resname
  integer:: nca,icon,imp1,imp2,imp1un,imp2un,iunit1,iunit2
  real,dimension(nmax,3):: rca
  real:: dist
  integer,dimension(nmax,2):: u
  character (LEN=3),dimension(nmax):: Lsequence
  real,dimension(nmax):: dist_nat 
  character(LEN=80):: char500
  character(LEN=80):: char200
  character(LEN=7):: nameid
  integer,dimension(nmax):: frame_in,k 
  integer,dimension(nmax):: icontact
  character(LEN=32):: arg1
     
!##############################-- some parameters --###########################################
  ncon11 = 271
  ncon12 = 77
  ncon22 = 51
  ncon = ncon11 + ncon12 + ncon22
  step_save = 100

  res_num1 = 81                 !|> number of residues for unit1 (KIX)
  res_num2 = 28                 !|> number of residues for unit2 (pKID)
  N_beads = res_num1 + res_num2 !|> Total number of residues 
  r_lim = 40.0
  q_thresh = 1
!###############################--argument passing n_frames--###################################
  do i = 1, iargc()
     call getarg(i, arg1)
     read(arg1,*) frame_n
  enddo
!==============================================================================================  

    
!####################-- This loop reads native info file, stores contact indecies and calculates
!####################-- n_contacts for helix A and B
  ll = 0
  
  open(unit=71,file = '1KDX.ninfo', status = 'unknown') 
  ninfo_lines:do
     read(71,'(a80)', END = 550) char500
     if(char500(1:7) /= 'contact') cycle ninfo_lines
     
     read(char500,'(a7,4x,i3,5x,i2,5x,i2,4x,i3,4x,i3,5x,i2,5x,i2,5x,f10.5)') nameid,icon,iunit1, &
          iunit2,imp1,imp2,imp1un,imp2un,go_nat_dist
     ll = ll + 1
     u(ll,1) = imp1
     u(ll,2) = imp2
     dist_nat(ll) = go_nat_dist

  end do ninfo_lines
  
550 continue
!==============================================================================================

  
!##### Here we read movie file and extract necessary info ####################################  
  state = 1
  n_esc = 0
  n_cap = 0
  
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


     
!!################ D_cm calculation #######################################################
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
     
!####################--D_cm between two units calculation--##################################################  
     d_cm = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))

     
     
!############--- Calculate all pairwise distances from frame and check with go_nat_dist---###################
     count = 0
     id = 0
     do i = ncon11 + 1, ncon11 + ncon12
        
        id = id + 1
        ii = u(i,1)
        jj = u(i,2)
        dx = rca(ii,1) - rca(jj,1)
        dy = rca(ii,2) - rca(jj,2)
        dz = rca(ii,3) - rca(jj,3)
        
        dist = dx*dx + dy*dy + dz*dz
        dist = sqrt(dist)
        
        if(dist <= 1.500 + dist_nat(i))then
           count = count + 1
           icontact(id) =  1
        endif
     end do


     
!!## Here we partition frames into different kinetic steps (encounters and finale binding) and print it out to analyse further with R 
     q = count
     if(state.eq.1.and.q.ge.q_thresh) then
        state = 2
        n_cap = n_cap + 1
     elseif(state.eq.2.and.d_cm.gt.r_lim) then
        state = 1
        n_esc = n_esc + 1
     endif

     !write(*,'(<ncon12>(i2,1x))')(icontact(i), i = 1, ncon12)
     if(q.ne.0) write(*,'(i3,2x,<ncon12>(i2,1x))') n_esc, (icontact(i), i = 1, ncon12) 
     
  end do

999 continue 

end program ExtractKinSteps
