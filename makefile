include makefile.inc

all:
	for i in $(SUBDIRS); do \
	echo "make all in $$i.."; \
	(cd $$i; $(MAKE) all); \
	$(CC) $(LLFLAGS) $(MAINFILE) $(OBJSDIR)$(O) $(INCSDIR) $(OUT) $(EXEC) $(LIBS); done

clean:
	for i in $(SUBDIRS); do \
	echo "clean all in $$i.."; \
	(cd $$i; $(MAKE) clean); done
	$(RM) $(EXEC)
	
moreclean:
	for i in $(TESTSUBDIRS); do \
	echo "clean all in $$i.."; \
	(cd $$i; $(MAKE) clean); done
	$(RM) $(EXEC)

run: all
	./$(EXEC)
	
runtests: clean
	for i in $(TESTSUBDIRS); do \
	echo "make all in $$i.."; \
	(cd $$i; $(MAKE) tests); done