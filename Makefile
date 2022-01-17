PHANTOMDIR=./phantom
EDITOR=vi
ifndef SETUP
SETUP=star
endif

again:
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME=; cd -; cp ${PHANTOMDIR}/bin/phantom .; cp ${PHANTOMDIR}/bin/phantom_version .

all: again setup infile

.PHONY: phantom phantomsetup phantom2power phantom2grid phantomanalysis phantomevcompare libphantom mflow
phantom         : again
phantomsetup    : setup
phantom2power   : power
phantom2grid    : grid
phantomanalysis : analysis
phantomevcompare: evcompare
phantomsinks    : sinks
libphantom      : phantomlib
mflow           : mflow

clean:
	cd ${PHANTOMDIR}; make clean KROME=krome
setup:
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME= setup; cd -; cp ${PHANTOMDIR}/bin/phantomsetup .
fixedSprofile:
	cp Makefile_phantom ${PHANTOMDIR}/build/Makefile; cp {fixed_S_profile.f90,write_fixedSprofile.f90} ${PHANTOMDIR}/src/setup; cd ${PHANTOMDIR}; make fixedSprofile SETUP=${SETUP} RUNDIR=${PWD}; cd -; cp ${PHANTOMDIR}/bin/fixedSprofile .
mike2:
	cd ${PHANTOMDIR}; make mike2 SETUP=${SETUP} RUNDIR=${PWD}; cd -; cp ${PHANTOMDIR}/bin/write_uTprofile .
mike3:
	cd ${PHANTOMDIR}; make mike3 SETUP=${SETUP} RUNDIR=${PWD}; cd -; cp ${PHANTOMDIR}/bin/calc_binding_energy .
moddump:
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME= moddump; cd -; cp ${PHANTOMDIR}/bin/phantommoddump .
analysis:
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME= phantomanalysis; cd -; cp ${PHANTOMDIR}/bin/phantomanalysis .
phantomlib:
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME= libphantom; cd -; cp ${PHANTOMDIR}/bin/libphantom.so .
pyanalysis:
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME= pyanalysis; cd -; cp ${PHANTOMDIR}/bin/libphantom.so .; ${PHANTOMDIR}/scripts/pyphantom/writepyanalysis.sh 
power:
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME= phantom2power; cd -; cp ${PHANTOMDIR}/bin/phantom2power .
grid:
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME= phantom2grid; cd -; cp ${PHANTOMDIR}/bin/phantom2grid .
evcompare:
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME= phantomevcompare; cd -; cp ${PHANTOMDIR}/bin/phantomevcompare .
sinks:
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME= phantomsinks; cd -; cp ${PHANTOMDIR}/bin/phantomsinks .
mflow:
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME= mflow; cd -; cp ${PHANTOMDIR}/bin/mflow .; cp ${PHANTOMDIR}/bin/ev2mdot .; cp ${PHANTOMDIR}/bin/lombperiod .
make:
	cd ${PHANTOMDIR}; ${EDITOR} build/Makefile &
qscript:
	@cd ${PHANTOMDIR}; make --quiet SETUP=${SETUP} RUNDIR=${PWD} KROME= qscript
%::
	cd ${PHANTOMDIR}; make SETUP=${SETUP} RUNDIR=${PWD} KROME= ${MAKECMDGOALS}; if [ -e ${PHANTOMDIR}/bin/$@ ]; then cd -; cp ${PHANTOMDIR}/bin/$@ .; fi