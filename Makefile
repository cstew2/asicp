TARGET		= asicp

CC		= g++
CFLAGS		:= -I/usr/include/eigen3/ -I/usr/include/nanoflann -std=c++17 -fopenmp
RCFLAGS		:= -Ofast -ftree-vectorize -finline-functions -march=native -flto
DCFLAGS		:= -O0 -g3 -ggdb

SRC		= $(wildcard src/*.cxx)
INC		= $(wildcard src/*.hxx)

LIBS		= -lgomp

$(shell mkdir -p obj)

OBJ		= $(subst src,obj,$(SRC:.cxx=.o))

.PHONY: release
release: CFLAGS		+= $(RCFLAGS)
release: $(TARGET)

.PHONY: debug
debug: CFLAGS		+= $(DCFLAGS)
debug: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

obj/%.o: src/%.cxx src/%.hxx
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS)

clean:
	rm -f $(TARGET) $(OBJ)
	rmdir obj
