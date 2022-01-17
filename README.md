## Getting the code
This code requires several modules from Phantom, which are incorporated as a git submodule to this repository. This means that this repository stores a copy of the Phantom repository, but the commit history of Phantom is independent from this code's. To clone this repository including the Phantom submodule, you must include the `--recurse-submodules` flag:
```
git clone --recurse-submodules git@github.com:themikelau/flat-entropy-star.git
```
You should see that the `flat-entropy-star` repository contains a copy of `phantom`.

## Setting up your calculation
Before compile the code, you should set up your calculation by editing `fixed_S_profile.f90`. Calculation settings are specified by editing part of the `calc_rho_and_pres` subroutine:
```
 !-----------------------------------------------------------------------------------------
 ! USER SETTINGS
 !-----------------------------------------------------------------------------------------
 !
 ! Choose total stellar mass (including core particle)
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
 Rstar = 4.3061478138863500d13 !500. * solarr
 !
 ! Desired surface pressure
 surfpres = 300. !2.374667d5
 !
 ! Adjustment factor for shooting
 rhofac0 = 0.01   ! Multiplicative factor for adjusting central density
 Sfac = 0.008     ! Additive factor for adjusting entropy
 tol = 1d-3       ! Relative tolerance for matching surface pressure and radius to desired values
 !
 ! Expression for entropy
 ieos = 12
 ientropy = 2 ! Include both gas and radiation entropy
 gmw = 0.61821 ! Assumed mean molecular weight
 !-----------------------------------------------------------------------------------------
```
In a nutshell, specify the stellar mass and radius and the core mass and radius. Choose an appropriate value for the star's surface pressure, which should not be "too high" to avoid expansion. You should also specify EoS-related settings: `ieos` (equation of state type, following what Phantom uses), `ientropy` (see the `entropy` function in Phantom's `eos` module), and `gmw` (mean molecular weight, assumed to be constant throughout star).

If may also adjust initial guesses for the central density and fixed value of specific entropy. Play with these values if the code is struggling to converge onto a solution:
```
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
```

## Compilation and execution
To compile the code, simply type
```
make fixedSprofile
```
Hopefully this is successful, and you should now have the executable `fixedSprofile`, which can be run by typing
```
./fixedSprofile
```
You must delete `fixedSprofile` and recompile every time you tweak the settings in `fixed_S_profile.f90`.