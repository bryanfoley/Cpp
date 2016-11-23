CC = g++
CCFLAGS = -c -Wall
INCS = -I../header
OBJS = Vector.o disk.o
RM = rm
O = *.o

all: $(OBJS)
	echo 'making objects...'
	
Vector.o:
	$(CC) $(CCFLAGS) Vector.cpp $(INCS)
	
disk.o:
	 $(CC) $(CCFLAGS) disk.cpp $(INCS)
	 
clean:
	$(RM) $(O)