
ifdef CLANG
CC = clang++
else
CC = g++
endif

CPPFLAGS = -std=c++17 -Wall -Wextra

ifdef CILKPLUS
CPPFLAGS += -DPARLAY_CILKPLUS -DCILK -fcilkplus
else ifdef OPENCILK
CPPFLAGS += -DPARLAY_OPENCILK -DCILK -fopencilk
else ifdef SERIAL
CPPFLAGS += -DPARLAY_SEQUENTIAL
else
CPPFLAGS += -pthread
endif

ifdef DEBUG
CPPFLAGS += -Og -mcx16 -DDEBUG -g
else ifdef PERF
CPPFLAGS += -Og -mcx16 -march=native -g
else ifdef MEMCHECK
CPPFLAGS += -Og -mcx16 -DPARLAY_SEQUENTIAL
else
CPPFLAGS += -O3 -mcx16 -march=native
endif

ifdef STDALLOC
CPPFLAGS += -DPARLAY_USE_STD_ALLOC
endif

EDIT_DISTANCE = edit_distance.h edit_distance_sequential.h edit_distance_dp.h
ALL = suffix_array_test test_framework

INCLUDE_DIR = parlaylib/examples/

all : $(ALL)

test_framework: test_framework.cpp edit_distance_dp.o edit_distance_parallel.o minimum_edit_distance.h dac_mm_k.h dac_mm.h
	$(CC) $(CPPFLAGS) test_framework.cpp -I$(INCLUDE_DIR) -o $@ edit_distance_dp.o edit_distance_parallel.o


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

clean :
	rm -f *.o $(ALL)
