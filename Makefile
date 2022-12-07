
ifdef CLANG
CC = clang++
else
CC = g++
endif

CPPFLAGS = -std=c++17 -O3 -pthread -Wall -Wextra

EDIT_DISTANCE = edit_distance.h edit_distance_sequential.h edit_distance_dp.h
ALL = suffix_array_test test_framework

INCLUDE_DIR = parlaylib/examples/

all : $(ALL)

test_framework: test_framework.cpp edit_distance_dp.o edit_distance_parallel.o dac_mm.o
	$(CC) $(CPPFLAGS) test_framework.cpp -I$(INCLUDE_DIR) -o $@ edit_distance_dp.o edit_distance_parallel.o dac_mm.o


suffix_array_test : suffix_array_test.o suffix_array_sequential.o
	$(CC) $(CPPFLAGS) -o $@ suffix_array_test.o suffix_array_sequential.o

# ------

suffix_array_test.o: suffix_array_test.cpp suffix_array_parallel.h range_min.h longest_common_prefix.h
	$(CC) $(CPPFLAGS) -c suffix_array_test.cpp

suffix_array_sequential.o: suffix_array_sequential.h suffix_array_sequential.cpp
	$(CC) $(CPPFLAGS) -c suffix_array_sequential.cpp


edit_distance_dp.o: edit_distance_dp.h edit_distance_dp.cpp
	$(CC) $(CPPFLAGS) -c edit_distance_dp.cpp

edit_distance_sequential.o: edit_distance_sequential.h edit_distance_sequential.cpp
	$(CC) $(CPPFLAGS) -c edit_distance_sequential.cpp

edit_distance_parallel.o: edit_distance_parallel.h edit_distance_parallel.cpp
	$(CC) $(CPPFLAGS) -c edit_distance_parallel.cpp

dac_mm.o: dac_mm.h dac_mm.cpp
	$(CC) $(CPPFLAGS) -c dac_mm.cpp

clean :
	rm -f *.o $(ALL)
