# Common definitions are placed in here
CC = g++
LD = g++
RM = rm
SHELL = /bin/sh
ECHO = echo

SUBDIRS =  source
OBJSDIR = source/
INCSDIR = -I./header/
O = *.o
H = *.h
OUT = -o

TESTDIRS = test/unitTests

MAINFILE = freeze_thaw.cpp
EXEC = freezethaw
TESTEXEC = testall
CCFLAGS = -c -Wall
EXECFLAGS = -o -Wall $(EXEC)
LLFLAGS = -Wall
INCS = -I$(SUBDIRS)
OBJS = Vector.o disk.o

.SILENT :