program rates
!####################################################################################
!## Program to calculate bimolecular rates and fraction of alignment upon successful
!## encounter simulation based on Q_int and dist .....
!## Reads txt file from std input
!####################################################################################
  implicit none
  integer,parameter:: nmax = 1100000
  integer:: frame_n, n_esc, n_cap,n_evol
  integer:: t_cap, t_esc, t_evol, step
  integer:: t_encoun, int_dummy 
  integer:: l,i,j,t,state,dt, indx_save
  real::    d_lim, d_Rcm, theta, d_b, d_cm, d_cm_a, d_cm_b, theta_a, theta_b
  integer:: q12,q22,qa, qb, q_b_min, q_thresh,q22_a, q22_b
  integer:: q12_alpha, q12_beta, pse_q12, pse_q12_alpha, pse_q12_beta, q_non, q_tot
  
  character(len=32) :: arg1
  
  d_lim = 40.0          !#### Threshold distance between units to to define encounter complex    
  step = 100            !#### step=t_step_save (Cafemol) 
  q_thresh = 1          !#### Threshold value for Q_int to define encounter complex
 !q_b_min = 0.70        !#### Value of Q_int at native bound state 

  do i = 1, iargc()
     call getarg(i, arg1) !! first argument should be q bound minima
     read(arg1,*) q_b_min
  enddo
  
  l = 0
  dt = 0
  n_esc = 0
  n_cap = 0
  n_evol = 0
 
  t_esc = 0 
  t_cap = 0
  t_evol = 0

  state = 1
  frame_lines:do
     l=l+1
     read(5,*, END=660) l, q12, q22, q12_alpha, q12_beta, pse_q12, pse_q12_alpha, pse_q12_beta, /
     q22_a, q22_b, d_cm, d_cm_a, d_cm_b, theta_a, theta_b, q_non, q_tot
     dt = dt + step
     if(state.eq.1.and.q12.ge.q_thresh) then
        state = 2
        t_cap = t_cap + dt
        n_cap = n_cap + 1
        dt = 0
        write(*,'(10i6,2x,5f6.2)') l, q12, q22, q12_alpha, q12_beta, pse_q12, pse_q12_alpha, pse_q12_beta,/
        q22_a, q22_b, d_cm, d_cm_a, d_cm_b, theta_a, theta_b, q_non, q_tot
     elseif(state.eq.2.and.q12.ge.q_b_min) then
        state = 3
        t_evol = t_evol + dt
        n_evol = n_evol + 1
        dt = 0
     elseif(state.eq.2.and.q12.lt.q_thresh) then
        state = 1
        t_esc = t_esc + dt
        n_esc = n_esc + 1
        dt = 0
     endif
 
  end do frame_lines
660 continue
  
end program rates
