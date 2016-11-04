include makefile.inc

all: freeze_thaw.o disk.o display.o 
	$(CC) -o freezethaw freeze_thaw.o disk.o display.o -I/usr/include/GL -I/usr/include/SDL2 -D_REENTRANT -L/usr/lib/x86_64-linux-gnu -lSDL2 -L/usr/lib/x86_64-linux-gnu -lGLEW -L/usr/lib/x86_64-linux-gnu/mesa -lGL

freeze_thaw.o: freeze_thaw.cpp disk.h display.h
	$(CC) -c freeze_thaw.cpp

disk.o: disk.cpp disk.h
	$(CC) -c disk.cpp
	
display.o: display.cpp display.h
	$(CC) -c display.cpp

clean:
	rm -f freezethaw *.o *~
	
run: all
	./freezethaw