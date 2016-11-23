include makefile.inc

all: freeze_thaw.cpp $(OBJS)
	$(CC) -o freezethaw freeze_thaw.cpp $(OBJS)

.cpp.o:
	$(CC) -c $< -o $@ $(INCS)

clean:
	rm *.o,$(EXEC)
	
run: all
	./$(EXEC)