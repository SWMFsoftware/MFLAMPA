# Copyright (C) 2002 Regents of the University of Michigan, 
# portions used with permission 
# For more information, see http://csem.engin.umich.edu/tools/swmf

DEFAULT_TARGET = MFLAMPA
DEFAULT_EXE    = ${DEFAULT_TARGET}.exe

default : ${DEFAULT_TARGET}

include Makefile.def
include Makefile.conf

# Menu of make options
help:
	@echo ' '
	@echo '  You can "make" the following:'
	@echo ' '
	@echo '    <default> ${DEFAULT_TARGET} in stand alone mode, help in SWMF'
	@echo ' '
	@echo '    help          (makefile option list)'
	@echo '    install       (install ${DEFAULT_TARGET})'
	@echo ' '
	@echo '    LIB           (Component library libSP for SWMF)'
	@echo '    ${DEFAULT_TARGET}       (Multiple-Field-Line Advection Model for Particle Acceleration)'
	@echo '    NOMPI         (NOMPI library for compilation without MPI)'
	@echo ' '
	@echo '    rundir        (create run directory for standalone or SWMF)'
	@echo '    rundir RUNDIR=run_test (create run directory run_test)'
	@echo ' '
	@echo '    test          (run all ${DEFAULT_TARGET} tests)'
	@echo ' '
	@echo '    clean         (remove temp files like: *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    distclean     (equivalent to ./Config.pl -uninstall)'


install: src/ModSize.f90

src/ModSize.f90: src/ModSize_orig.f90
	cp -f src/ModSize_orig.f90 src/ModSize.f90

LIB:    install
	cd src;          make LIB
	cd srcInterface; make LIB

${DEFAULT_TARGET}:
	cd ${SHAREDIR}; ${MAKE} LIB
	cd ${TIMINGDIR}; ${MAKE} LIB
	cd src; ${MAKE} LIB
	cd src; ${MAKE} ${DEFAULT_TARGET}

NOMPI:
	cd util/NOMPI/src; make LIB

COMPONENT = SP

rundir:
	mkdir -p ${RUNDIR}/SP 
	cd ${RUNDIR}/SP; \
		mkdir restartIN restartOUT IO2; \
		ln -s ${SPDIR}/Param .
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		touch ${DIR}/share/JobScripts/job._TMP_${MACHINE}; \
		touch ${DIR}/share/JobScripts/_TMP_.${MACHINE}.pl; \
		cp ${DIR}/share/JobScripts/job.*${MACHINE}* ${RUNDIR}/; \
		cp ${DIR}/share/JobScripts/*.${MACHINE}.pl ${RUNDIR}/; \
		rm -f ${RUNDIR}/*_TMP_* ${DIR}/share/JobScripts/*_TMP_*; \
		cp -f Param/PARAM.test ${RUNDIR}/PARAM.in; \
		touch ${RUNDIR}/core; chmod 444 ${RUNDIR}/core; \
		cd ${RUNDIR}; ln -s ${BINDIR}/${DEFAULT_EXE} .; \
		ln -s SP/* .; \
	fi);

clean:  install
	@(if [ -r "Makefile.conf" ]; then \
		cd src; make clean; \
		cd ../srcInterface; make clean; \
	fi)

distclean: 
	./Config.pl -uninstall

allclean:
	cd src; $(MAKE) distclean
	cd srcInterface; $(MAKE) distclean
	cd Doc/Tex; $(MAKE) cleanpdf

# Testing

TESTDIR = run_test

test: test_mflampa

test_mflampa:
	@echo "test_mflampa_compile..." > test_mflampa.diff
	${MAKE} test_mflampa_compile
	@echo "test_mflampa_rundir..." >> test_mflampa.diff
	${MAKE} test_mflampa_rundir
	@echo "test_mflampa_run..." >> test_mflampa.diff
	${MAKE} test_mflampa_run
	@echo "test_mflampa_check..." >> test_mflampa.diff
	${MAKE} test_mflampa_check

test_mflampa_compile:
	./Config.pl -g=20000
	${MAKE} 

test_mflampa_rundir: 
	rm -rf ${TESTDIR}
	${MAKE} rundir RUNDIR=${TESTDIR} STANDALONE=YES SPDIR=`pwd`
	cd ${TESTDIR}; cp -f Param/PARAM.in.test PARAM.in
	cd ${TESTDIR}; tar xzf ../data/input/test_mflampa/MH_data_e20120123.tgz

test_mflampa_run:
	cd ${TESTDIR}; ${MPIRUN} ./MFLAMPA.exe | tee runlog

test_mflampa_check:
	cat ${TESTDIR}/SP/IO2/MH_data_{*???_???,*n000006}.out \
		> ${TESTDIR}/SP/IO2/MH_data.outs
	${SCRIPTDIR}/DiffNum.pl -t -r=1e-6 -a=1e-6 \
		${TESTDIR}/SP/IO2/MH_data.outs \
		data/output/test_mflampa/MH_data.ref.gz \
		> test_mflampa.diff
	ls -l test_mflampa.diff
