PROGRAM=bem2d.exe
CXX=icpc --std=c++11
OBJS=bem2d.o bessel.o

.SUFFIXES: .c .o

$(PROGRAM): $(OBJS)
	$(CXX) $(OBJS) -o $(PROGRAM)

.c.o:
	$(CXX) $(CFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(PROGRAM) $(OBJS)
