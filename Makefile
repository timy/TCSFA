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
	$(RM) ana/data/plot/*~
	$(RM) ana/data/proc/app
	git add .
	git commit -a -m '1. options for representation-dependent problem: new W-representation (Lagranian-based Coulomb correction) is added for discussion; 2. information output is explicitly controlled in the configure file; 3. too-close-to-core situation is raised, which may cause high-order term which spoil the spectra; one can control the radius-threshold to discard trajectories which are too close to the core, to somewhat alleviate the rampant artificial higher-order caustics, and to significantly speed up the computation'
	git push origin master