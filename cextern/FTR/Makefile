#
#  Makefile
#  FTR
#
#  Created by Alexander Rudy on 2015-08-13.
#  Copyright 2015 Alexander Rudy. All rights reserved.
#

LDFLAGS = -g -O2 -lfftw3 -lm -L /opt/local/lib
OBJECTS = ftr.o ftr_example.o
TARGET = libftr.so ftr-example
BUILDDIR = build/
SRCDIR = src/
EXAMPLE = example/

$(BUILDDIR)%.o : $(SRCDIR)%.c $(SRCDIR)%.h
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILDDIR)%.o : $(EXAMPLE)%.c
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

lib%.so : $(BUILDDIR)%.o
	$(CC) $(LDFLAGS) -shared -o $@ $<
	
ftr-example : $(BUILDDIR)ftr_example.o libftr.so
	$(CC) $(LDFLAGS) -L. -lftr -o $@ $<

all: $(TARGET)

clean:
	rm -r $(BUILDDIR)
	rm $(TARGET)
