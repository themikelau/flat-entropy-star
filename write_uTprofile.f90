program write_uTprofile
 use rho_profile, only:write_softened_profile,read_mesa
 use eos,         only:calc_temp_and_ene,ieos,gmw,X_in,Z_in,init_eos,gamma
 implicit none
 real, allocatable, dimension(:) :: m,r,rho,pres,ene,temp,Xfrac,Yfrac
 real :: Mstar
 integer :: ierr,i
 character(len=120) :: inputpath,outputpath

 !-Settings--------
 inputpath = 'fixedSprofile.dat'
 outputpath = 'fixedSprofile_mesa.dat'
 ieos = 10
 gmw = 0.61821
 gamma = 1.6666666667
 if (ieos == 10) then
    X_in = 0.69843
    Z_in = 0.01426
 endif 
 !--------------------------

 call read_mesa(inputpath,rho,r,pres,m,ene,temp,Xfrac,Yfrac,Mstar,ierr,cgsunits=.true.)
 call init_eos(ieos,ierr)
 
 if (ieos == 2) print*,'Using ieos = ',ieos,'gamma = ',gamma,' gmw = ',gmw
 
 do i=1,size(m)
    call calc_temp_and_ene(rho(i),pres(i),ene(i),temp(i),ierr)
 enddo

 call write_softened_profile(outputpath,m,pres,temp,r,rho,ene)
 print*,'Profile written to',outputpath

end program write_uTprofile