!###################################################################
!                    Talant Ruzmetov  written on 11/21/2013
!                                     updated on 07/21/2016
! This code reads movie + native info file and calculates
! reaction coordinates such as q (fraction of native contacts for
! folding and binding) 
!###################################################################

program ContactMap
  implicit none
  integer,parameter::nmax = 400
  real:: x, y, z, dx, dy, dz
  real:: go_nat_dist
  real:: count, count11, count12, count22, count_alpha,count_beta
  integer:: q11, q12, q22, q12_alpha, q12_beta, q22_a, q22_b
  integer:: pse_q12, pse_q12_beta, pse_q12_alpha
  integer:: i, j, ii, l, jj, i1, ll
  integer:: frame_n, ncon11, ncon12, ncon22, ncon, ncon_alpha, ncon_beta

  integer:: res_num1, res_num2, res_num, n_link, n_alpha, n_beta, ncon_a, ncon_b
  real:: pse_count, pse_count11, pse_count12, pse_count22, pse_count_alpha, pse_count_beta, count_a, count_b

  character(LEN=72):: char72,char55
  character(LEN=4)::nameofatom
  character(LEN=3):: resname
  integer::nca,icon,imp1,imp2,imp1un,imp2un,iunit1,iunit2
  real,dimension(nmax,3)::rca, rca2, rca1
  real::dist
  integer,dimension(nmax,2)::u
  character (LEN=3),dimension(nmax)::Lsequence
  real,dimension(nmax)::dist_nat,dist_nat11,dist_nat12,dist_nat22
  character(LEN=80):: char500
  character(LEN=80):: char200
  character(LEN=7)::nameid
  integer,dimension(nmax)::frame_in,k
  character(LEN=32)::arg1
     
  !##############################-- some parameters --###########################################
  ncon11 = 271
  ncon12 = 77
  ncon22 = 51
  ncon = ncon11 + ncon12 + ncon22

  n_link = 2    ! linker length
  n_alpha = 11  ! helix-A size 
  n_beta = 15   ! helix-B size
  res_num1 = 81 ! KIX size
  res_num2 = 28 ! pKID size

!###############################--argument passing n_frames--###################################
  do i = 1, iargc()
     call getarg(i, arg1)
     read(arg1,*) frame_n
  enddo
!==============================================================================================  


  
!####################-- This loop reads native info file, stores contact indecies and calculates
!####################-- n_contacts for helix A and B
  ll = 0
  ncon_alpha = 0
  ncon_beta = 0
  ncon_a = 0
  ncon_b = 0
  
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

     !#########################--counting tot # of contacts for binding--#############     
     if(iunit1.eq.1.and.iunit2.eq.2) then
        if(imp2un.le.n_alpha) ncon_alpha = ncon_alpha + 1
        if(imp2un.gt.n_alpha + n_link) ncon_beta = ncon_beta+1
     endif
     
!#########################--counting tot # of contacts for binding--#############     
   !  if(ll.gt.ncon11.and.ll.le.ncon11+ncon12) then
   !     if(imp2.le.res_num1+n_alpha) ncon_alpha = ncon_alpha + 1
   !     if(imp2.ge.res_num1+n_alpha+n_link) ncon_beta = ncon_beta + 1
   !  endif
     
  end do ninfo_lines
  
550 continue
!==============================================================================================

  
!##### Here we read movie file and extract necessary info ####################################  
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
    
   
!############---- Calculate distance from frame and check with go_nat_dist---######################
  count = 0
  count12 = 0
  count22 = 0
  count_alpha = 0
  count_beta = 0
  count_a = 0
  count_b = 0
  pse_count = 0
  pse_count12 = 0
  pse_count_alpha = 0
  pse_count_beta = 0
  
  do i1 = 1, ncon
     
     ii = u(i1,1)
     jj = u(i1,2)
     dx = rca(ii,1) - rca(jj,1)
     dy = rca(ii,2) - rca(jj,2)
     dz = rca(ii,3) - rca(jj,3)
     
     dist = dx*dx + dy*dy + dz*dz
     dist = sqrt(dist)

     !# counting intermolecular contacts
     if(dist <= 1.500+dist_nat(i1))then
        count = count + 1        
        if(i1.gt.ncon11.and.i1.le.ncon11 + ncon12) then
           count12 = count12 + 1
           if(jj.le.res_num1 + n_alpha) count_alpha = count_alpha + 1
           if(jj.gt.res_num1 + n_alpha + n_link) count_beta = count_beta + 1
        endif
     endif

     !# counting intramolecular contacts
     if(dist <= 1.200*dist_nat(i1)) then
        if(i1.gt.ncon11 + ncon12) then 
           count22 = count22 + 1
           if(ii.le.res_num1 + n_alpha) count_a = count_a + 1
           if(jj.gt.res_num1 + n_alpha + n_link) count_b = count_b + 1
        endif
     endif
     
     !# counting intermolecular pseudo contacts
     if(dist.gt.1.500+dist_nat(i1)) then
        if(dist.lt.dist_nat(i1)+4.500) then
           pse_count = pse_count + 1        
           if(i1.gt.ncon11.and.i1.le.ncon11+ncon12) then
              pse_count12 = pse_count12 + 1
              if(jj.le.res_num1 + n_alpha) pse_count_alpha = pse_count_alpha + 1
              if(jj.gt.res_num1 + n_alpha + n_link) pse_count_beta = pse_count_beta + 1
           endif
        endif
     endif
     
  end do
  q12 = count12 !/real(ncon12)
  q22 = count22 !/real(ncon22)
  q12_alpha = count_alpha !/real(ncon_alpha)
  q12_beta = count_beta !/real(ncon_beta)
  pse_q12 = pse_count12 !/real(ncon12)
  pse_q12_alpha = pse_count_alpha !/real(ncon_alpha)
  pse_q12_beta = pse_count_beta !/real(ncon_beta)
  q22_a = count_a !/real(ncon_a)
  q22_b = count_b !/real(ncon_b)
  
  write(*,'(10i7)') l, q12, q22, q12_alpha, q12_beta, pse_q12, pse_q12_alpha, pse_q12_beta, q22_a, q22_b
end do
999 continue 

end program contactmap
