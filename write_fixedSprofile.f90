program write_fixed_S_profile
 use fixedSprofile
 use rho_profile, only:write_softened_profile
 use eos, only:calc_temp_and_ene
 implicit none
 real, allocatable, dimension(:) :: m,r,rho,pres,ene,temp
 character(len=120)              :: outputpath
 integer                         :: ierr,i

 call calc_rho_and_pres(m,r,rho,pres)
 allocate(ene(0:size(m)-1), temp(0:size(m)-1))

 print*,'CALCULATING THERMAL PROFILE'
 ! Calculate thermal profile
 do i = 0,size(m)-1
     call calc_temp_and_ene(rho(i),pres(i),ene(i),temp(i),ierr)
 enddo

 outputpath = 'fixedSprofile.dat'
 call write_softened_profile(outputpath,m,pres,temp,r,rho,ene)
 print*,'PROFILE WRITTEN TO',outputpath

 print*,'May this envelope stay put.'

end program write_fixed_S_profile