program calc_binding_energy
 use rho_profile, only:read_mesa
 use eos,         only:calc_temp_and_ene,ieos,gmw,X_in,Z_in,init_eos,gamma
 use kernel,      only:kernel_softening
 use physcon,     only:solarm,solarr,gg
 implicit none
 real, allocatable, dimension(:) :: m,r,rho,pres,ene,temp,Xfrac,Yfrac
 real :: Mstar,hsoft,rcore,mcore,dum,q,Ebind_b,Ebind_g,potencore,phi
 integer :: ierr,i
 character(len=120) :: inputpath

 !-Settings--------
 inputpath = 'fixedSprofile.dat'
 rcore = 18.5 * solarr
 mcore = 3.8405 * solarm
 ieos = 2
 gmw = 0.61821
 gamma = 1.6666666667
 if (ieos == 10) then
    X_in = 0.69843
    Z_in = 0.01426
 endif 
 !--------------------------

 call read_mesa(inputpath,rho,r,pres,m,ene,temp,Xfrac,Yfrac,Mstar,ierr,cgsunits=.true.)
 if (ierr == 1) print*,'Error: Cannot read file'
 call init_eos(ieos,ierr)
 
 if (ieos == 2) print*,'Using ieos = ',ieos,'gamma = ',gamma,' gmw = ',gmw
 
 Ebind_g = 0.
 Ebind_b = 0.
 do i=2,size(m)
    call calc_temp_and_ene(rho(i),pres(i),ene(i),temp(i),ierr)

    ! Get gravitational potential
    hsoft = 0.5 * rcore
    q = r(i) / hsoft
    call kernel_softening(q**2,q,potencore,dum)
    potencore = potencore / hsoft ! Note: gcore is not multiplied by G or mcore yet, and is already negative.
    phi = gg * ( -m(i)/r(i) + mcore*potencore )
    
!    if (r(i) > rcore) then
       Ebind_b = Ebind_b + (phi + ene(i)) * ( m(i) - m(i-1) )
       Ebind_g = Ebind_g + phi * ( m(i) - m(i-1) )
!    endif
 enddo

 print*,'Ebind_b = ', Ebind_b
 print*,'Ebind_g = ', Ebind_g
 

end program calc_binding_energy