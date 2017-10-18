!## Talant Ruzmetov 05/03/2013
!## This code calculates ceter of mass distance between pKID and KIX
!## and alignment of helix_B end to end vector with respect to native
!## orientation. It reads movie file from str input.

program dist_calculation
  implicit none
  integer,parameter:: nmax = 200
  real:: x, y, z, d_cm 
  integer:: i, j, l, n, frame_n, N_beads, ires
  integer:: res_num1, res_num2, res_num
  character(LEN=72):: char72
  character(LEN=4):: nameofatom
  character(LEN=3):: resname
  real,dimension(nmax,3):: rca
  real,dimension(3):: Rcm_pro, Rcm_lig, dr
  character (LEN=3),dimension(nmax):: Lsequence
  character(LEN=32)::arg1
  
  !###############################-- some parameters--##############################################
  res_num1 = 81                 !|> number of residues for unit1 (KIX)
  res_num2 = 28                 !|> number of residues for unit2 (pKID)
  N_beads = res_num1 + res_num2 !|> Total number of residues 

  
  !####################--argument passing from term--################################################
  do j = 1, iargc()
     call getarg(j, arg1)
     read(arg1,*) frame_n  !|> total number of frames for trajectory file
  enddo
      
  do l = 1, frame_n
     ires = 0
     Lines: do
        read(5,'(a72)',end=20) char72
        if(char72(1:3)=='END')then 
           exit Lines
        endif
        
        if(char72(1:4) /= 'ATOM') cycle Lines
        read(char72,10) nameofatom,resname,x,y,z
        if(nameofatom .eq. ' CA ') then
           ires = ires + 1
           rca(ires,1) = x
           rca(ires,2) = y
           rca(ires,3) = z
           Lsequence(ires) = resname
        endif
10    format(12x,a4,1x,a3,10x,3f8.3)
     end do Lines
     
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

     write(*,*) d_cm
  end do
20 continue 
end program dist_calculation
