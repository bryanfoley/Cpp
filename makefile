include makefile.inc

all:
	for i in $(SUBDIRS); do \
	echo "make all in $$i.."; \
	(cd $$i; $(MAKE) all); \
	$(CC) $(LLFLAGS) $(MAINFILE) $(OBJSDIR)$(O) $(INCSDIR) $(OUT) $(EXEC); done

clean:
	for i in $(SUBDIRS); do \
	echo "clean all in $$i.."; \
	(cd $$i; $(MAKE) clean); done
	$(RM) $(EXEC)
	
run: all
	./$(EXEC)
	
tests: $(TESTEXEC)
	echo 'This is where the tests run'