program calc_binding_energy
 use rho_profile, only:read_mesa
 use eos,         only:calc_temp_and_ene,ieos,gmw,X_in,Z_in,init_eos,gamma
 use kernel,      only:kernel_softening
 use physcon,     only:solarm,solarr,gg,radconst,kb_on_mh
 implicit none
 real, allocatable, dimension(:) :: m,r,rho,pres,ene,temp,Xfrac,Yfrac
 real :: Mstar,hsoft,rcore,mcore,dum,q,potencore
 real :: bind_int,bind_grav,bind_th,Egas,Erad,Erec,phi,ethi,egasi,eradi,ereci,mui,dmi
 integer :: ierr,i
 logical :: iwritefile
 character(len=120) :: inputpath,outpath

 !-Settings--------
 inputpath = 'fixedSprofile.dat'
 iwritefile = .true.
 outpath = 'binding_energy_grav.dat'
 rcore = 18.5 * solarr
 mcore = 3.8405 * solarm
 ieos = 12
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
 Egas = 0.
 Erad = 0.
 Erec = 0.

 do i=2,size(m)
    call calc_temp_and_ene(rho(i),pres(i),ene(i),temp(i),ierr)

    ! Get gravitational potential
    hsoft = 0.5 * rcore
    q = r(i) / hsoft
    call kernel_softening(q**2,q,potencore,dum)
    potencore = potencore / hsoft ! Note: gcore is not multiplied by G or mcore yet, and is already negative.
    phi = gg * ( -m(i)/r(i) + mcore*potencore )

    ! Get (specific) energies depending on EoS
    select case (ieos)
    case (2) ! Ideal gas
       ethi = ene(i) ! Thermal energy is same as internal energy
       egasi = ethi   ! Gas thermal energy makes up entire thermal energy
       eradi = 0.
       ereci = 0.
    case (12) ! Ideal gas plus radiation
       ethi = ene(i) ! Thermal energy is same as internal energy
       egasi = 1.5 * kb_on_mh * temp(i) / gmw
       eradi = radconst * temp(i)**4 / rho(i)
       ereci = 0.
    case (10) ! We want just the gas + radiation internal energy
       ! Get mu from pres and temp
       mui = rho(i) * kb_on_mh * temp(i) / (pres(i) - radconst * temp(i)**4 / 3.)
       egasi = 1.5 * kb_on_mh * temp(i) / mui ! 3/2*kT/(mu*mh)
       eradi = radconst * temp(i)**4 / rho(i) ! aT^4/rho
       ethi = egasi + eradi
       ereci = ene(i) - ethi ! Remaining component is due to ionisation
    end select

   ! if (r(i) > rcore) then
       dmi = m(i) - m(i-1)
       bind_grav = bind_grav + dmi * phi
       bind_th   = bind_th   + dmi * (phi + ethi)
       bind_int  = bind_int  + dmi * (phi + ene(i))
       Egas = Egas + dmi * egasi
       Erad = Erad + dmi * eradi
       Erec = Erec + dmi * ereci
   ! endif

    if (iwritefile) then
       if (i==2) then
          open(unit=42, file=outpath, status='replace')
          write(42,"(a,2x,a,2x,a)") '       r / cm','        m / g','  Ebind / erg'
          write(42,"(es13.6,2x,es13.6,2x,es13.6)") r(i),m(i)+mcore,bind_int
       else
          open(unit=42, file=outpath, position='append')
          write(42,"(es13.6,2x,es13.6,2x,es13.6)") r(i),m(i)+mcore,bind_int
       endif
    endif
 enddo
 close(unit=42)
 

 print*,'bind_grav = ', bind_grav
 print*,'bind_th   = ', bind_th 
 print*,'bind_int  = ', bind_int
 print*,'Egas      = ', Egas
 print*,'Erad      = ', Erad
 print*,'Erec      = ', Erec
 print*,'CONSISTENCY CHECKS:'
 print*,'bind_grav + Egas               = ', bind_grav + Egas
 print*,'bind_grav + Egas + Erad        = ', bind_grav + Egas + Erad
 print*,'bind_grav + Egas + Erad + Eint = ', bind_grav + Egas + Erad + Erec

end program calc_binding_energy