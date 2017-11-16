program rates
!###########################################################################################
!## Program to calculate folding rate from long equilibrium simulation based on
!## Q_intra 
!## Reads txt file from std input
!###########################################################################################
  implicit none
  integer:: frame_n, n_fo, n_un,non, n_off
  integer:: t_fo, step
  integer:: t_un 
  integer:: l,i,t,f_state,dt,dtau,count_IF,count_CS
  real::    q_fold,q_unfold
  real::    q11
  character(len=32) :: arg1

  step = 200            !#### step=t_step_save (Cafemol) 
  q_fold = 0.75         !#### Threshold Q value for folding
  q_unfold = 0.30       !#### Threshold Q value for unfolding
  
  l = 0
  dt = 0
  dtau = 0
  
  t_un = 0
  t_fo = 0
  n_un = 0
  n_fo = 0

  f_state = 1
  frame_lines:do
     l=l+1
     read(5,*, END=660) q11       
     
!!###################### Folding rate calculation  #######################!!     
   
        dtau = dtau + step
        if(f_state.eq.1.and.q11.lt.q_unfold) then
           f_state = 2
           t_un = t_un + dtau
           n_un = n_un + 1
           dtau = 0
        elseif(f_state.eq.2.and.q11.gt.q_fold) then
           f_state = 1
           t_fo = t_fo + dtau
           n_fo = n_fo + 1
           dtau = 0
           
        endif
!!***********************************************************************!!
!!***********************************************************************!!
     
  end do frame_lines
660 continue
  
  write(6,'(8i11)') t_fo, t_un, n_fo, n_un  
end program rates
