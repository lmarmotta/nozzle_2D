#
# Name of the exec file.
TARGET=nozzle_2d

#
# Link flags abreviation.
LIBS=-lm -lcgns

#
# Path to hdf5 parallel compiler (needs compilation first).
CC=gcc

# 
# Optimization or debug flags.
CFLAGS=-O0 -g -Wall -Werror

#
# CGNS include path.
INCLUDE=/home/leomm/opt/cgns/include

#
# CGNS lib path.
LDFLAGS=/home/leomm/opt/cgns/lib64

#
# ----------- DO NOT MODIFY -----------  
#

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS=$(patsubst %.c, %.o, $(wildcard *.c))
HEADERS=$(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(INCLUDE) -c $< -o $@ -L$(LDFLAGS) $(LIBS)

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -I$(INCLUDE) -o $@ -L$(LDFLAGS) $(LIBS)

clean:
	-rm -f *.o
	-rm -f $(TARGET)


