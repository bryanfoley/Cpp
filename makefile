include makefile.inc

all: freeze_thaw.cpp disk.cpp
	$(CC) -o freezethaw freeze_thaw.cpp disk.cpp

clean:
	rm -f *.o
	
run: all
	./freezethaw