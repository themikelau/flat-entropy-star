program calc_binding_energy
 use rho_profile, only:read_mesa
 use eos,         only:calc_temp_and_ene,ieos,gmw,X_in,Z_in,init_eos,gamma
 use kernel,      only:kernel_softening
 use physcon,     only:solarm,solarr,gg,radconst,kb_on_mh
 implicit none
 real, allocatable, dimension(:) :: m,r,rho,pres,ene,temp,Xfrac,Yfrac
 real :: Mstar,hsoft,rcore,mcore,dum,q,bind_int,bind_grav,bind_th,potencore,phi,ethi,mui
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
 
 bind_grav = 0.
 bind_th = 0.
 bind_int = 0.
 do i=2,size(m)
    call calc_temp_and_ene(rho(i),pres(i),ene(i),temp(i),ierr)

    ! Get gravitational potential
    hsoft = 0.5 * rcore
    q = r(i) / hsoft
    call kernel_softening(q**2,q,potencore,dum)
    potencore = potencore / hsoft ! Note: gcore is not multiplied by G or mcore yet, and is already negative.
    phi = gg * ( -m(i)/r(i) + mcore*potencore )
    

    ! Get (specific) thermal energy
    select case (ieos)
    case (2,12) ! Ideal gas or ideal gas plus radiation
       ethi = ene(i) ! Thermal energy is same as internal energy
    case (10) ! We want just the gas + radiation internal energy
       ! Get mu from pres and temp
       mui = rho(i) * kb_on_mh * temp(i) / (pres(i) - radconst * temp(i)**4 / 3.)
       ethi = 1.5 * kb_on_mh * temp(i) / mui + radconst * temp(i)**4 / rho(i) ! 3/2*kT/(mu*mh) + aT^4/rho
    end select

    if (r(i) > rcore) then
       bind_int = bind_int + (phi + ene(i)) * ( m(i) - m(i-1) )
       bind_th = bind_th + (phi + ethi) * ( m(i) - m(i-1) )
       bind_grav = bind_grav + phi * ( m(i) - m(i-1) )
    endif
 enddo

 print*,'bind_grav = ', bind_grav
 print*,'bind_th   = ', bind_th 
 print*,'bind_int  = ', bind_int

end program calc_binding_energy