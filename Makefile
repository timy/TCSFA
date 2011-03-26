include ./Makefile.cfg

ifeq ($(PARALLEL), 1)
DIR_ENTRY	=	p_entry
DIR_MMFF	=	./src/mmff
else
DIR_ENTRY	=	s_entry
DIR_MMFF	=	
endif

SUBDIRS = $(DIR_MMFF) $(DIR_MAIN)/$(DIR_ENTRY) 


all: $(SUBDIRS)
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done 

run:
ifeq ($(PARALLEL), 1)
	qsub run.sh -q opt.q
else
	./$(EXE_NAME)
endif

clean:
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done
	$(RM) *~ 
	$(RM) ccsfa_traj.* 

cleanall: clean
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done
	$(RM) $(EXE_NAME)
	$(RM) dat/*.dat
	$(RM) LOG

commit:
	$(RM) ana/data/*.dat
	$(RM) ana/plot/*~
	$(RM) ana/proc/app
	git add .
	git commit -a -m 'p_entry/console_slave.f90: set_pulse function have 5 parameters intead of 4'
	git push origin master