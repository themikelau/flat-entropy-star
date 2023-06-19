program write_uTprofile
 use rho_profile, only:read_mesa
 use eos,         only:calc_temp_and_ene,ieos,gmw,X_in,Z_in,init_eos,gamma,equationofstate
 use units,       only:set_units,unit_density,unit_velocity,unit_ergg
 use physcon,     only:solarm,solarr
 implicit none
 real, allocatable, dimension(:) :: m,r,rho,pres,ene,temp,Xfrac,Yfrac,csound
 real :: Mstar,spsoundi,ponrhoi,eni
 integer :: ierr,i
 character(len=120) :: inputpath,outputpath
 logical :: icsound

 !-Settings--------
 inputpath = 'fixedSprofile.dat'
 outputpath = 'fixedSprofile_ideal_cs.dat'
 icsound = .true. ! Include X, Y, cs columns
 ieos = 2
 gmw = 0.61821
 gamma = 1.6666666667
 if (ieos == 10) then
    X_in = 0.69843
    Z_in = 0.01426
 endif 
 !--------------------------

 call set_units(solarr,solarm,G=1.)
 call read_mesa(inputpath,rho,r,pres,m,ene,temp,Xfrac,Yfrac,Mstar,ierr,cgsunits=.true.)
 if (ierr == 1) print*,'Error: Cannot read file'
 call init_eos(ieos,ierr)
 
 if (ieos == 2) print*,'Using ieos = ',ieos,'gamma = ',gamma,' gmw = ',gmw
 
 do i=1,size(m)
    call calc_temp_and_ene(ieos,rho(i),pres(i),ene(i),temp(i),ierr)
 enddo
 if (icsound) then
    allocate(csound(size(m)))
    do i=1,size(m)
       eni = ene(i)/unit_ergg
       call equationofstate(ieos,ponrhoi,spsoundi,rho(i)/unit_density,1.,1.,1.,eni) ! x,y,z are dummies
       csound(i) = spsoundi * unit_velocity
       Xfrac(i) = X_in
       Yfrac(i) = 1.-X_in-Z_in
    enddo
    call write_mesa(outputpath,m,pres,temp,r,rho,ene,Xfrac,Yfrac,csound)
 else
    call write_mesa(outputpath,m,pres,temp,r,rho,ene)
 endif
 print*,'Profile written to',outputpath

end program write_uTprofile