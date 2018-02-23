CC=g++
CFLAGS= -std=c++11 -lgsl -lgslcblas -lm -lboost_system -lboost_filesystem -fPIC
PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
OBJ = tandem.o
DEPS = conversions.h tandem.h cholesky.cpp tandem_function_interface.h
# location of the Boost Python include files and library
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib
# compile mesh classes
TARGET = PyTandem
PWD=/nfs/users3/emiljoha/venv-e/ex-job/ANNAN/tandem
# LD_LIBRARY_PATH=$(PWD)
# The mess known as compiling shared library for python bindings.
$(TARGET).so: $(TARGET).o $(DEPS) tandem.o
	g++ -shared -Wl,--export-dynamic,-rpath=$(PWD) $(TARGET).o tandem.o -L$(PWD) -L$(BOOST_LIB) -lboost_python \
-L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION) -o $(TARGET).so $(CFLAGS) \


$(TARGET).o: $(TARGET).cpp $(DEPS)
	g++ -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c $(TARGET).cpp $(CFLAGS)

libtandem.so : tandem.o
	gcc -shared -o libtandem.so tandem.o $(CFLAGS)

tandem.o: tandem.cpp $(DEPS)
	g++ -c -o $@ $< $(CFLAGS)

cli : tandem_cli.cpp tandem_function_interface.h tandem.o
	$(CC) -I. -c tandem_cli.cpp --std=c++11
	$(CC) -o tandem tandem_cli.o tandem.o -lgsl -lgslcblas -lm -lboost_filesystem -lboost_system

debug : tandem_cli.cpp tandem.h tandem.o
	$(CC) -I. -c tandem_cli.cpp --std=c++11	-g
	$(CC) -o tandem tandem_cli.o tandem.o -lgsl -lgslcblas -lm -lboost_filesystem -lboost_system -g

cholesky: cholesky.cpp cholesky.h
	g++ -I. -Wall -c -g cholesky.cpp --std=c++11
	g++ -o cholesky.out cholesky.o -lgsl -lgslcblas -lm
	./cholesky.out
test:
	make -C ./tests
clean:
	rm *.o
	rm *.so
	rm *~

