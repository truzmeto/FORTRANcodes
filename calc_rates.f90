program rates
  implicit none
  integer,parameter::nmax = 2000000
  integer:: frame_n,l,a,t_cap,t_esc,step,n_esc,n_cap,t_evol
  real :: b,c,r_lim,q,r,q_b_min,q_thresh,dum, p_bound, p_unbound
  character(len=50):: input_file_name
  character(len=32):: arg1
!  real,dimension(nmax):: q,r
!  integer,dimension(nmax):: t
  integer tot_time,t,state,dt,n_evol,i
  integer t_encoun,t_bo,t_un
  
  r_lim = 40.000 
  t_evol = 0
  t_esc = 0 
  t_cap = 0
  n_esc = 0
  n_cap = 0
  n_evol = 0
  step = 100
  l = 0
  dt = 0
  state = 1
  q_thresh = 1.0 / 77.0

  do i = 1, iargc()
     call getarg(i, arg1) !! first argument should be q bound minima
     read(arg1,*) q_b_min
  enddo
     
  t_encoun = 0
  t_un = 0
  t_bo = 0
  
  frame_lines:do
     l = l + 1
!     int_dummy, q12, q22, qa, qb, dum, dum, dum, dum, dum, d_Rcm, dum, d_b, dum, theta 
     read(5,*, END = 660) t,q,dum,dum,dum,dum,dum,dum,dum,dum,r,dum,dum !input should be given using <
     
     dt = dt + step
     if(state.eq.1.and.q.ge.q_thresh) then
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
         
!!--------------check pop of bound and unbound in the encounter---------------!!!!
    
        if(q.eq.0.000) then
           t_un = t_un + 1
        else
           t_bo = t_bo + 1
        endif
     
!!-------------------------------------------------------------------------!!!!!!
  end do frame_lines
660 continue

  p_bound = t_bo / real(l)
  p_unbound = t_un / real(l)
  
  write(6,'(6i10,2x,2f5.2)') t_cap, t_esc, t_evol, n_cap, n_esc, n_evol, p_unbound, p_bound
  
end program rates
   
         
   
