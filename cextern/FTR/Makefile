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
TARGET = $(LIBDIR)libftr.so ftr-example
OBJECTS = $(BUILDDIR)ftr.o $(BUILDDIR)slopemanage.o

$(BUILDDIR)%.o : $(SRCDIR)%.c $(SRCDIR)%.h
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(EXAMPLE)%.o : $(EXAMPLE)%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBDIR)lib%.so : $(OBJECTS)
	@mkdir -p $(LIBDIR)
	$(CC) $(LDFLAGS) -shared -o $@ $<
	
ftr-example : $(EXAMPLE)ftr_example.o $(LIBDIR)libftr.so
	$(CC) $(LDFLAGS) -L$(LIBDIR) -lftr -o $@ $<

all: $(TARGET)

clean:
	rm -r $(BUILDDIR)
	rm $(TARGET)
