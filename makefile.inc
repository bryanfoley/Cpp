# Common definitions are placed in here
CC = g++
LD = g++
RM = rm
SHELL = /bin/sh
ECHO = echo

CFLAGS = -c -Wall
LFLAGS = -Wall
INCS = -I./
OBJS = disk.o Vector.o
EXEC = freezethaw

.SILENT :