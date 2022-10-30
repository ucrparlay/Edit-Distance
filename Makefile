
ifdef CLANG
CC = clang++
else
CC = g++
endif

CPPFLAGS = -std=c++17 -O3 -pthread -Wall -Wextra -Werror

EDIT_DISTANCE = edit_distance.h edit_distance_sequential.h edit_distance_dp.h
ALL = suffix_array_test edit_distance_test

all : $(ALL)


suffix_array_test : suffix_array_test.o suffix_array_sequential.o
	$(CC) $(CPPFLAGS) -o $@ suffix_array_test.o suffix_array_sequential.o

edit_distance_test: edit_distance_test.o edit_distance_sequential.o edit_distance_dp.o
	$(CC) $(CPPFLAGS) -o $@ edit_distance_test.o edit_distance_sequential.o edit_distance_dp.o

# ------

suffix_array_test.o: get_time.hpp suffix_array_test.cpp
	$(CC) $(CPPFLAGS) -c suffix_array_test.cpp

suffix_array_sequential.o: suffix_array_sequential.h suffix_array_sequential.cpp
	$(CC) $(CPPFLAGS) -c suffix_array_sequential.cpp

edit_distance_test.o: get_time.hpp edit_distance_test.cpp
	$(CC) $(CPPFLAGS) -c edit_distance_test.cpp

edit_distance_dp.o: edit_distance_dp.h edit_distance_dp.cpp
	$(CC) $(CPPFLAGS) -c edit_distance_dp.cpp

edit_distance_sequential.o: edit_distance_sequential.h edit_distance_sequential.cpp
	$(CC) $(CPPFLAGS) -c edit_distance_sequential.cpp


clean :
	rm -f *.o $(ALL)
