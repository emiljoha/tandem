CC=g++
CFLAGS=-I. -std=c++11 -g
# DEPS = ../tandem.h ../cholesky.h test.h
#../tandem_function_interface.h, ../pairs.h
DEPS=../*.h
OBJ = tandem.o test.o
LINKFLAGS= -lgsl -lgslcblas -lm -lboost_filesystem -lboost_system
# TESTS = test_case_0.out test_case_2.out test_case_4.out \
# random_test_cases.out test_cholesky.out	timer.out zero_spin.out
TESTS = zero_spin.out hf_test.out two_way_test.out energy_test.out tandem_test.out

all: $(TESTS) HF.wf one_mat.m

run: $(TESTS) HF.wf one_mat.m
	./zero_spin.out
	./hf_test.out
	./two_way_test.out
	./energy_test.out
	./tandem_test.out
# ./test_case_0.out
# ./test_case_2.out
# ./test_case_4.out
# ./random_test_cases.out
# ./test_cholesky.out
# ./timer.out -p 4 -e 200 -f test.cho -d uniform

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
