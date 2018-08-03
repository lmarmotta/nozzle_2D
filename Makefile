#
# Name of the exec file.
TARGET=code

#
# Link flags abreviation.
LIBS=-lm -lcgns -lhdf5

#
# Path to hdf5 parallel compiler (needs compilation first).
CC=/home/leonardo/opt/hdf5/bin/h5pcc

# 
# Optimization or debug flags.
CFLAGS=-O0 -g -Wall

#
# CGNS include path.
INCLUDE=/home/leonardo/opt/cgns/include

#
# CGNS lib path.
LDFLAGS=/home/leonardo/opt/cgns/lib64

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
	-rm -f *.o *.dat
	-rm -f $(TARGET)


