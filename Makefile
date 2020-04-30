TARGET		= asicp

CC		= g++
CFLAGS		:= -I/usr/include/eigen3/ -I/usr/include/nanoflann -std=c++17 -fopenmp
RCFLAGS		:= -O3 -ftree-vectorize -finline-functions -march=native
DCFLAGS		:= -O0 -g3 -ggdb

CFLAGS		+= $(RCFLAGS)

SRC		= $(wildcard *.cxx)
INC		= $(wildcard *.hxx)

OBJ		= $(SRC:.cxx=.o)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.cxx
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	rm -f $(TARGET) $(OBJ)
