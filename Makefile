include ./Makefile.cfg

ifeq ($(PARALLEL),1)
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
	./$(EXE_NAME)
	#qsub run.sh  -q opt.q

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
	git commit -a -m '1. New files: analyze_energy_spec_2D.f90 - angular-resolved energy spectra; 2. CEP is included in the pulse.'
	git push origin master