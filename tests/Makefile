CC=g++
CFLAGS=-I. -std=c++11 -g
# DEPS = ../tandem.h ../cholesky.h test.h
#../tandem_function_interface.h, ../pairs.h
DEPS=../*.h
OBJ = tandem.o test.o
LINKFLAGS= -lgsl -lgslcblas -lm -lboost_filesystem -lboost_system
# TESTS = zero_spin.out timer.out
TESTS = zero_spin.out hf_test.out two_way_test.out energy_test.out tandem_test.out hf_energy_test.out \
	random_test_cases.out test_cholesky.out timer.out

all: $(TESTS) HF.wf one_mat.m

run: $(TESTS) HF.wf one_mat.m
	./zero_spin.out
	./hf_test.out
	./two_way_test.out
	./energy_test.out
	./tandem_test.out
	./hf_energy_test.out
	./random_test_cases.out
	./test_cholesky.out
	./timer.out

%.out: %.o $(OBJ)
	g++ -o $@ $^ $(LINKFLAGS)

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

tandem.o: ../tandem.cpp $(DEPS)
	g++ -c -o $@ $< $(CFLAGS)

HF.wf one_mat.m : two_way_test.py
	python two_way_test.py

timer.out : timer.cpp ../tandem.h
	$(CC) -I. -c timer.cpp --std=c++11 -g
	$(CC) -o timer.out timer.o -lgsl -lgslcblas -lm -lboost_filesystem -lboost_system tandem.o

clean:
	rm HF.wf one_mat.m
	rm *.out
	rm *~
