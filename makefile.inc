# Common definitions are placed in here
CC = g++
LD = g++
RM = rm
SHELL = /bin/sh
ECHO = echo

SUBDIRS =  source
OBJSDIR = source/
INCSDIR = -I./header/ -I/usr/include/GL -I/usr/include/SDL2
LIBS = -L/usr/lib/x86_64-linux-gnu -lSDL2 -L/usr/lib/x86_64-linux-gnu -lGLEW -L/usr/lib/x86_64-linux-gnu/mesa -lGL
O = *.o
H = *.h
OUT = -o

TESTSUBDIRS = source test/unitTests

MAINFILE = freeze_thaw.cpp
EXEC = freezethaw
CCFLAGS = -c -Wall
EXECFLAGS = -o -Wall $(EXEC)
LLFLAGS = -Wall
INCS = -I$(SUBDIRS)
OBJS = Vector.o disk.o

.SILENT :