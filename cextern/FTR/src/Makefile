#
#  Makefile
#  FTR
#
#  Created by Alexander Rudy on 2015-08-13.
#  Copyright 2015 Alexander Rudy. All rights reserved.
#

LDFLAGS = -g -O2 -lfftw3 -L /opt/local/lib
OBJECTS = ftr.o
TARGET = ftr.so

%.o : %.c %.h 
	$(CC) $(CFLAGS) -c $< -o $@

%.so : %.o
	$(CC) $(LDFLAGS) -shared -o $@ $<

all: $(TARGET)
