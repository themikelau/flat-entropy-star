program calc_binding_energy
 use readwrite_mesa,only:write_mesa
 use eos,         only:calc_temp_and_ene,ieos,gmw,X_in,Z_in,init_eos,gamma
 use eos_gasradrec, only:irecomb
 use kernel,      only:kernel_softening
 use physcon,     only:solarm,solarr,gg,radconst,kb_on_mh
 implicit none
 character(len=120) :: inputpath,outpath
 integer :: ierr,i,iunit
 logical :: iwritefile
 real, allocatable, dimension(:) :: m,r,rho,pres,ene,temp,Xfrac,Yfrac
 real :: Mstar,hsoft,rcore,mcore,dum,q,potencore,guesseint = 0.,&
         bind_int,bind_grav,bind_th,Egas,Erad,Erec,phi,ethi,egasi,eradi,ereci,mui,dmi

 !-Settings--------
 inputpath = 'fixedS_gasrad.dat'
 iwritefile = .true.
 outpath = 'Ebind.dat'
 rcore = 18.5 * solarr
 mcore = 3.8405 * solarm
 ieos = 20
 gmw = 0.61821  ! Only used for ieos = 2,12
 gamma = 1.6666666667
 irecomb = 3
 if ( (ieos == 10) .or. (ieos == 20) ) then
    X_in = 0.69843
    Z_in = 0.01426
 endif 
 !--------------------------

 ! Only density and pressure are used
 call read_mesa(inputpath,rho,r,pres,m,ene,temp,Xfrac,Yfrac,Mstar,ierr,cgsunits=.true.)
 if (ierr == 1) print*,'Error: Cannot read file'
 call init_eos(ieos,ierr)

 bind_grav = 0.
 bind_th = 0.
 bind_int = 0.
 Egas = 0.
 Erad = 0.
 Erec = 0.
 ethi = 0.
 eradi = 0.
 egasi = 0.
 ereci = 0.

 do i = 2,size(m)
    mui = gmw
    guesseint = 0.
    call calc_temp_and_ene(ieos,rho(i),pres(i),ene(i),temp(i),ierr,guesseint,mui,X_in,Z_in)

    ! Get gravitational potential
    hsoft = 0.5 * rcore
    q = r(i) / hsoft
    call kernel_softening(q**2,q,potencore,dum)
    potencore = potencore / hsoft ! Note: gcore is not multiplied by G or mcore yet, and is already negative.
    phi = gg * ( -m(i)/r(i) + mcore*potencore )

    ! Get (specific) energies depending on EoS
    select case (ieos)
    case(2) ! Ideal gas
       ethi = ene(i) ! Thermal energy is same as internal energy
       egasi = ethi   ! Gas thermal energy makes up entire thermal energy
       eradi = 0.
       ereci = 0.
    case(12) ! Ideal gas plus radiation
       ethi = ene(i) ! Thermal energy is same as internal energy
       egasi = 1.5 * kb_on_mh * temp(i) / gmw
       eradi = radconst * temp(i)**4 / rho(i)
       ereci = 0.
    case(10,20) ! We want just the gas + radiation internal energy
       ! Get mu from pres and temp
       mui = rho(i) * kb_on_mh * temp(i) / (pres(i) - radconst * temp(i)**4 / 3.)
       egasi = 1.5 * kb_on_mh * temp(i) / mui ! 3/2*kT/(mu*mh)
       eradi = radconst * temp(i)**4 / rho(i) ! aT^4/rho
       ethi = egasi + eradi
       ereci = ene(i) - ethi ! Remaining component is due to ionisation
    end select

!    if (r(i) > rcore) then
       dmi = m(i) - m(i-1)
       bind_grav = bind_grav + dmi * phi
       bind_th   = bind_th   + dmi * (phi + ethi)
       bind_int  = bind_int  + dmi * (phi + ene(i))
       Egas = Egas + dmi * egasi
       Erad = Erad + dmi * eradi
       Erec = Erec + dmi * ereci
!    endif

    if (iwritefile) then
       if (i==2) then
          open(newunit=iunit, file=outpath, status='replace')
          write(iunit,"(a,2x,a,2x,a,2x,a,2x,a)") '       r / cm','        m / g','Ebind_g / erg','Ebind_t / erg','Ebind_i / erg'
       else
          open(newunit=iunit, file=outpath, position='append')
       endif
       write(iunit,"(es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6)") r(i),m(i)+mcore,bind_grav,bind_th,bind_int
    endif
    close(unit=iunit)
 enddo
 

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
