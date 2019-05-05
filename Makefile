PROG		= asicp

CC		= g++
CFLAGS		= -I/usr/include/eigen3 -Wpedantic -Wall -g3 -ggdb

SRC		= $(wildcard *.cxx)
INC		= $(wildcard *.h)

OBJ		= $(SRC:.cxx=.o)

$(PROG): $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.cxx
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	rm -f $(PROG) $(OBJ)
