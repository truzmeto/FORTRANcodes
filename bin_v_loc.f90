!#########----- Program written by ------###############
!#########----- Talant Ruzmetov----------###############
!#########----- July, 2015---------------###############
!#########----- Benchmarking: 200Mb data frame is evaluated in 5s time

program BinLocalVelocities
  implicit none
  integer,parameter:: nmax = 400
  integer:: i, j, ncon, ncon11, ncon12, l
  integer:: res_num, nbin, ibin, nchain
  real, dimension(nmax)::  v_local
  integer, dimension(nmax):: bin_count
  real:: bin_width, q_int_max, q_int_min, v_cm
  integer:: q_int
  real, dimension(nmax,nmax):: v_res
    
!####################### manual parameters ###############################
  ncon12 = 77
  res_num = 28

 
!##### binnig parameters  
  q_int_max = real(ncon12)
  q_int_min = 0.0
  bin_width = 1.0
  nbin = int((q_int_max - q_int_min)/bin_width) + 1
  nchain = res_num
  
  v_res(1:nbin,1:nchain) = 0.0
  bin_count(1:nbin) = 0

  l = 0
  Lines: do
     !##### reading lines from data frame  
     read(5,*,end = 30) q_int, v_cm, (v_local(i), i = 1, nchain)
     l = l + 1

     !# looping over bins and collecting counts
     do ibin = 1, nbin
        if(q_int.ge.(ibin - 1)*bin_width.and.q_int.lt.ibin*bin_width) then
           bin_count(ibin) = bin_count(ibin) + 1

           !# loop over residues along chain and do summation per bin 
           do i = 1, nchain
              v_res(ibin,i) = v_res(ibin,i) + v_local(i)
           enddo

        endif
     enddo
     
  end do Lines
30 continue
  
  
  !## printing the output to a file
  do i = 1, nchain
     write(6,'(3x,<nbin>(f8.3,1x))')(v_res(ibin,i)/real(bin_count(ibin)), ibin = 1, nbin)
     !write(6,'(3x,78(f8.3,1x))')(v_res(ibin,i)/real(bin_count(ibin)), ibin = 1, nbin) !# for gfortran
  enddo
  
end program BinLocalVelocities
    
