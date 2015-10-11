#
#  Makefile
#  FTR
#
#  Created by Alexander Rudy on 2015-08-13.
#  Copyright 2015 Alexander Rudy. All rights reserved.
#

LDFLAGS = -g -O2 -lfftw3 -lm -L /opt/local/lib
BUILDDIR = build/
SRCDIR = src/
LIBDIR = lib/
EXAMPLE = example/
INCLUDEDIR = include/
INCLUDE = $(INCLUDEDIR)ftr.h $(INCLUDEDIR)slopemanage.h
TARGET = $(LIBDIR)libftr.so ftr-example $(INCLUDE)
OBJECTS = $(BUILDDIR)ftr.o $(BUILDDIR)slopemanage.o

$(BUILDDIR)%.o : $(SRCDIR)%.c $(INCLUDE)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -I$(INCLUDEDIR) -c $< -o $@

$(EXAMPLE)%.o : $(EXAMPLE)%.c
	$(CC) $(CFLAGS) -I$(INCLUDEDIR) -c $< -o $@

$(LIBDIR)lib%.so : $(OBJECTS)
	@mkdir -p $(LIBDIR)
	$(CC) $(LDFLAGS) -shared -o $@ $^
	
ftr-example : $(EXAMPLE)ftr_example.o $(LIBDIR)libftr.so
	$(CC) $(LDFLAGS) -L$(LIBDIR) -I$(INCLUDEDIR) -lftr -o $@ $<

$(INCLUDEDIR)%.h : $(SRCDIR)%.h
	@mkdir -p $(INCLUDEDIR)
	cp $< $(INCLUDEDIR)

all: $(TARGET)

clean:
	rm -f $(TARGET)
	rm -rf $(LIBDIR)
	rm -rf $(BUILDDIR)
	rm -rf $(INCLUDEDIR)
