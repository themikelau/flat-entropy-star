module fixedSprofile
 implicit none

contains

subroutine calc_rho_and_pres(m,r,rho,pres)
 use physcon, only:solarm,solarr,kb_on_mh
 use eos,only:init_eos,entropy,gmw,ieos
 use rho_profile, only:write_softened_profile
 real, allocatable, dimension(:), intent(out) :: m,r,rho,pres
 real, allocatable, dimension(:) :: dm,dm_c
 real :: Mstar,Rstar,R_old,mcore,rcore,Sc,surfpres,surfpres_old,rhofac,Sfac,rhofac0,tol
 integer :: i,N,ientropy,counter1,counter2,ierr
 logical :: irefined
 character(len=120) :: outputpath
 ! INSTRUCTIONS
 ! * The stellar centre has index 0 while the surface has index N
 ! * r(i), m(i) are cell face values
 ! * pres(i), rho(i) are cell centre values
 ! * dm(i) := m(i) - m(i-1) is the mass width of a cell
 ! * dm_c(i) := 0.5*( dm(i) + dm(i-1) ) is the distance between cell centres m(i-1/2) and m(i+1/2)

 !-----------------------------------------------------------------------------------------
 ! USER SETTINGS
 !-----------------------------------------------------------------------------------------
 !
 ! Choose stellar total mass (including core particle)
 Mstar = 12.0 * solarm
 !
 ! Choose number of zones
 N = 10000
 !
 ! Choose core mass
 mcore = 3.84048 * solarm
 !
 ! Choose core radius
 rcore = 18.9 * solarr
 !
 ! Choose desired stellar radius to shoot for
 Rstar = 4.3061478138863500d13 !500. * solarr !3.416975d13
 !
 ! Desired surface pressure
 surfpres = 300. !2.374667d5
 !
 ! Adjustment factor for shooting
 rhofac0 = 0.01 ! Multiplicative factor for adjusting central density
 Sfac = 0.008     ! Additive factor for adjusting entropy
 tol = 1d-3     ! Relative tolerance for matching surface pressure and radius to desired values
 !
 ! Expression for entropy
 ieos = 12
 ientropy = 2 ! Include both gas and radiation entropy
 gmw = 0.61821
 !-----------------------------------------------------------------------------------------

 print*,'CONSTRUCTING FIXED ENTROPY STAR'
 write(*,'(a30,e13.5,a10)'),'M =',Mstar / solarm,'Msun'
 write(*,'(a30,e13.5,a10)'),'mcore =',mcore / solarm,'Msun'
 write(*,'(a30,e13.5,a10)'),'rcore =',rcore / solarr,'Rsun'
 write(*,'(a30,e13.5,a10)'),'Desired surface pressure =',surfpres,'dyn cm^-2'
 write(*,'(a30,e13.5,a10)'),'Desired radius =',Rstar / solarr,'Rsun'
 call init_eos(ieos,ierr)

 ! Allocate and set up arrays
 allocate( m(0:N), rho(0:N), pres(0:N), r(0:N) )
 allocate( dm(1:N), dm_c(2:N))

 do i = 1,N
    dm(i) = (Mstar - mcore) / N
    dm_c(i) = dm(i)
 enddo

 do i = 0,N
    m(i) = real(i) * dm(1) ! note: this shortcut only works for dm = constant
 enddo

 !-----------------------------------------------------------------------------------------
 ! INITIAL GUESSES
 !-----------------------------------------------------------------------------------------
 !
 ! Central density
 rho(0) = 7.607091d-5
 !
 ! Central pressure
 pres(0) = 2.831573d10       ! Isn't actually used
 !
 ! Constant value of entropy
 Sc = 4569322041.73678
 !-----------------------------------------------------------------------------------------
 
 call set_central_bc(Sc,m,r,rho,pres,ientropy)

 !---------------------------LOOP-OVER-ENTROPY-VALUE---------------------------------------
 ! Vary the value of the fixed entropy to get the desired radius
 r(N) = 1.1 * Rstar
 rhofac = rhofac0
 counter2 = 0

 do
    R_old = r(N) ! Save value of stellar radius in previous shot
    rhofac = rhofac0 ! Reset fractional change in rho(0) to its starting value
    pres(N) = 1.1 * surfpres
    counter1 = 0

    if (irefined) then
       write(*,'(a4,2x,a4,2x,a12,2x,a12,2x,a12,2x,a12,2x,a12,2x,a12)'),'c1','c2','rhofac','S','rho_0','Psurf','R','rho refined' ! DEBUGGING
       irefined = .false.
    else
       write(*,'(a4,2x,a4,2x,a12,2x,a12,2x,a12,2x,a12,2x,a12)'),'c1','c2','rhofac','S','rho_0','Psurf','R' ! DEBUGGING
    endif
    !---------------------------LOOP-OVER-CENTRAL-DENSITY-------------------------------------
    ! Vary central density to get the desired surface pressure
       do
          surfpres_old = pres(N) ! Save value of surface pressure in previous shot
          call one_shot(Sc,mcore,rcore,m,dm,dm_c,r,rho,pres,ientropy,.false.)
          if ( abs(pres(N)-surfpres) < tol * pres(N) ) exit
          write(*,'(i4,2x,i4,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4)'),counter1,counter2,rhofac,Sc,rho(0),pres(N),r(N)/solarr ! DEBUGGING
          if (pres(N) > surfpres) then
             ! Decrease value of central density to decrease surface pressure
             rho(0) = rho(0) * (1 - rhofac)
          else
             ! Increase value of central density to increase surface pressure
             rho(0) = rho(0) * (1 + rhofac)
          endif
          if ( (surfpres - surfpres_old) * (surfpres - pres(N)) < 0. ) then
             rhofac = rhofac * 0.5 ! Refine adjustment to Sc
             irefined = .true.
          endif
          call set_central_bc(Sc,m,r,rho,pres,ientropy) ! Set central boundary conditions for new central density
          counter1 = counter1 + 1
         !  outputpath = 'testprof.dat'
         !  call write_softened_profile(outputpath,m,pres,m,r,rho,m)
       enddo
    !-----------------------------------------------------------------------------------------
    if (abs(Rstar-r(N)) < tol*Rstar) exit
    if (r(n) < Rstar) then
       ! Increase value of constant entropy to increase radius
       Sc = Sc + kb_on_mh * Sfac / gmw
    else
       ! Decrease value of constant entropy to decrease radius
       Sc = Sc - kb_on_mh * Sfac / gmw
    endif
    if ( (Rstar - r(n)) * (Rstar - r(n)) < 0. ) then
       Sfac = 0.5 * Sfac !log(0.5) + log(10.**Sfac + 1.) ! Refine adjustment to Sc
       print*, 'Sfac refined'
    endif
    call set_central_bc(Sc,m,r,rho,pres,ientropy) ! Set central boundary conditions for new Sc
    counter2 = counter2 + 1
 enddo
 !-----------------------------------------------------------------------------------------

 print*,'CONVERGED ONTO DESIRED PROFILE'
 write(*,'(a30,e13.5,a10)'),'Obtained surface pressure =',pres(N),'dyn cm^-2'
 write(*,'(a30,e13.5,a10)'),'Obtained radius =',r(N) / solarr,'Rsun'

end subroutine calc_rho_and_pres


subroutine one_shot(Sc,mcore,rcore,m,dm,dm_c,r,rho,pres,ientropy,idebug)
 use physcon, only:gg,pi,solarm,solarr
 use eos, only:get_rho_from_p_s
 use rho_profile, only:write_softened_profile
 real, intent(in)                               :: Sc,mcore,rcore
 real, allocatable, dimension(:), intent(in)    :: m,dm,dm_c
 real, allocatable, dimension(:), intent(inout) :: r,rho,pres
 integer, intent(in)                            :: ientropy
 integer                                        :: i,N
 logical, intent(in)                            :: idebug
 character(len=120) :: outputpath

 N = size(m)-1
 do i = 1,N-1
    if (idebug) then
      write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4)') i,m(i)/solarm,r(i),rho(i),pres(i)
    endif

    pres(i+1) = pres(i) - dm_c(i) * gg/(4*pi*r(i)**2) * ( m(i)/r(i)**2 + mcore * gcore(r(i),rcore) )
    if (pres(i+1) < 0.) then
       print*,'ERROR: Reached pres < 0 at i = ',i, 'm = ',m(i)/solarm
       write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4)') i,pres(i),pres(i+1),r(i-1),r(i),m(i)/solarm
       !outputpath = 'oneshot.dat'
       !call write_softened_profile(outputpath,m(0:i),pres(0:i),m(0:i),r(0:i),rho(0:i),m(0:i))
       pres(N) = -1. ! Force shooter to increase central density again and refine adjustment to central pressure
       return
    endif
    if (pres(i+1) > pres(i)) then
       print*,'ERROR: Pressure inversion at i = ',i, 'm = ',m(i)/solarm
       write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4)') i,pres(i),pres(i+1),r(i-1),r(i),m(i-1)/solarm,m(i)/solarm
       stop
    endif
    call get_rho_from_p_s(pres(i+1),Sc,rho(i+1),rho(i),ientropy)
    if (rho(i+1) < 0.) then
       print*,'ERROR: Reached rho < 0 at i = ',i, 'm = ',m(i)/solarm
       write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4)') i,pres(i+1),rho(i),Sc
       stop
    endif
    if (rho(i+1) > rho(i)) then
       print*,'ERROR: Density inversion at i = ',i, 'm = ',m(i)/solarm
       write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4)') i,pres(i+1),rho(i),Sc
       stop
    endif
    r(i+1) = (3*dm(i+1)/(4*pi*rho(i+1)) + r(i)**3)**(1./3.)
    if (r(i+1) < 0.) then
       print*,'ERROR: Reached r < 0 at i = ',i, 'm = ',m(i)/solarm
       write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4)') i,rho(i+1),r(i),r(i+1)
       stop
    endif
    if (r(i+1) < r(i)) then
       print*,'ERROR: Radius inversion at i = ',i, 'm = ',m(i)/solarm
       write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4)') i,rho(i+1),r(i),r(i+1)
       stop
    endif
 enddo

end subroutine one_shot


!-----------------------------------------------------------------------
!+
!  Set central boundary conditions
!+
!-----------------------------------------------------------------------
subroutine set_central_bc(Sc,m,r,rho,pres,ientropy)
 use eos, only:get_pres_from_rho_s,get_rho_from_p_s
 use physcon, only:gg,pi
 real, intent(in) :: Sc
 real, allocatable, dimension(:), intent(in) :: m
 real, allocatable, dimension(:), intent(inout) :: r,rho,pres
 integer, intent(in) :: ientropy
 real :: presguess

 presguess = pres(0)
 call get_pres_from_rho_s(rho(0),Sc,pres(0),presguess,ientropy) ! Get pres(0)
 r(1) = ( 3.*m(1) / (4.*pi*rho(0)) )**(1./3.) ! This is the condition of m(0) = m(1/2), i.e. drho/dm = 0 at m(0)
 pres(1) = pres(0) - 3.*gg/(8.*pi) * (4.*pi/3 * rho(0))**(4./3.) * m(1)**(2./3.) ! Ref: Kippenhahn & Weigert
 call get_rho_from_p_s(pres(1),Sc,rho(1),rho(0),ientropy) ! Get rho(1) from pres(1)

end subroutine set_central_bc


!-----------------------------------------------------------------------
!+
!  Acceleration from softened potential of stellar core
!+
!-----------------------------------------------------------------------
function gcore(r,rcore)
 use kernel, only:kernel_softening
 real, intent(in) :: r,rcore
 real             :: gcore,dum
 real             :: hsoft,q

 hsoft = 0.5 * rcore
 q = r / hsoft
 call kernel_softening(q**2,q,dum,gcore)
 gcore = gcore / hsoft**2 ! Note: gcore is not multiplied by G or mcore yet.

end function gcore


end module fixedSprofile