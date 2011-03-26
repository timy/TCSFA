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
	$(RM) ana/data/*~
	$(RM) ana/plot/*~
	$(RM) ana/plot/*.pyc
	$(RM) ana/proc/app
	$(RM) doc/*~

commit: cleanall
	$(RM) ana/data/*.dat
	git add .
	git commit -a -m 'bug fix'
	git push origin master