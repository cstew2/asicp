TARGET		= asicp

CC		= g++
CFLAGS		= -I/usr/include/eigen3 -I/usr/include/ -std=c++17 -O0 -g3 -ggdb

SRC		= $(wildcard *.cxx)
INC		= $(wildcard *.hxx)

OBJ		= $(SRC:.cxx=.o)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.cxx
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	rm -f $(PROG) $(OBJ)
