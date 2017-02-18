CC = g++ -std=c++11
CCFLAGS = -c -Wall
LLFLAGS = -Wall
INCSDIR = -I../../header
OBJS = VectorTests.o
OBJSDIR = ../../source/
GTESTLIB = -L /usr/src/gtest
LIBS = -lgtest -lgtest_main -lpthread
TESTFILE = AllTests.cpp
TESTEXEC = AllTests
RM = rm
O = *.o
H = *.h

tests: $(OBJS)
	echo 'making test objects...'
	$(CC) $(LLFLAGS) $(TESTFILE) $(OBJSDIR)$(O) $(INCSDIR) -o $(TESTEXEC) $(GTESTLIB) $(LIBS)
	./$(TESTEXEC)
	
VectorTests.o:
	$(CC) $(CCFLAGS) VectorTests.cpp $(INCSDIR)

clean:
	$(RM) $(O) $(TESTEXEC)