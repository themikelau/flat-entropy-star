PHANTOMDIR=./phantom
SETUP=star

fixedSprofile:
	cp Makefile_phantom ${PHANTOMDIR}/build; cp {fixed_S_profile.f90,write_fixedSprofile.f90} ${PHANTOMDIR}/src/setup; cd ${PHANTOMDIR}/build; make -f Makefile_phantom fixedSprofile SETUP=${SETUP} RUNDIR=${PWD}; cd -; cp ${PHANTOMDIR}/bin/fixedSprofile .
writeuT:
	cp Makefile_phantom ${PHANTOMDIR}/build; cp write_uTprofile.f90 ${PHANTOMDIR}/src/setup; cd ${PHANTOMDIR}/build; make -f Makefile_phantom writeuT SETUP=${SETUP} RUNDIR=${PWD}; cd -; cp ${PHANTOMDIR}/bin/write_uTprofile .
calcbinding:
	cp Makefile_phantom ${PHANTOMDIR}/build; cp calc_binding_energy.f90 ${PHANTOMDIR}/src/setup; cd ${PHANTOMDIR}/build; make -f Makefile_phantom calcbinding SETUP=${SETUP} RUNDIR=${PWD}; cd -; cp ${PHANTOMDIR}/bin/calc_binding_energy .