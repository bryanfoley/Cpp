CC = g++
CCFLAGS = -c -Wall
INCS = -I../header
OBJS = Vector.o disk.o display.o
TESTOBJS = Vector.o
RM = rm
O = *.o
TESTEXEC = alltests

tests: $(TESTOBJS)
	echo 'making test objects...'

all: $(OBJS)
	echo 'making objects...'	
	
Vector.o:
	$(CC) $(CCFLAGS) Vector.cpp $(INCS)
	
disk.o:
	 $(CC) $(CCFLAGS) disk.cpp $(INCS)
	 
display.o:
	$(CC) $(CCFLAGS) display.cpp $(INCS)
	 
clean:
	$(RM) $(O)