#
#  Makefile
#  FTR
#
#  Created by Alexander Rudy on 2015-08-13.
#  Copyright 2015 Alexander Rudy. All rights reserved.
#

LDFLAGS = -g -O3 -Wall -lfftw3_threads -lfftw3 -lm -lpthread -L /opt/local/lib
CCFLAGS = -g -O3 -Wall
BUILDDIR = build/
SRCDIR = src/
LIBDIR = lib/
EXAMPLE = example/
INCLUDEDIR = include/
HEADERS = $(wildcard $(SRCDIR)*.h)
INCLUDE = $(patsubst $(SRCDIR)%.h,$(INCLUDEDIR)%.h, $(HEADERS))
PROGS_C = $(wildcard $(EXAMPLE)*.c)
PROGS = $(subst _,-, $(patsubst $(EXAMPLE)%.c,%, $(PROGS_C)))
TARGET = $(LIBDIR)libftr.so $(INCLUDE) $(PROGS)
OBJECTS = $(BUILDDIR)ftr.o $(BUILDDIR)slopemanage.o

all: $(TARGET)

$(BUILDDIR)%.o : $(SRCDIR)%.c $(INCLUDE)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(EXAMPLE)%.o : $(EXAMPLE)%.c $(INCLUDE)
	$(CC) $(CFLAGS) -I$(INCLUDEDIR) -c $< -o $@

$(LIBDIR)lib%.so : $(OBJECTS)
	@mkdir -p $(LIBDIR)
	$(CC) $(LDFLAGS) -shared -o $@ $^
	
ftr-example : $(EXAMPLE)ftr_example.o $(EXAMPLE)clock.o $(LIBDIR)libftr.so
	$(CC) $(LDFLAGS) -L$(LIBDIR) -I$(INCLUDEDIR) -lftr -o $@ $< $(EXAMPLE)clock.o

recon-timing : $(EXAMPLE)recon_timing.o $(EXAMPLE)clock.o $(LIBDIR)libftr.so
	$(CC) $(LDFLAGS) -L$(LIBDIR) -I$(INCLUDEDIR) -lftr -lopenblas -o $@ $< $(EXAMPLE)clock.o

aperture-test : $(EXAMPLE)aperture_test.o $(LIBDIR)libftr.so
	$(CC) $(LDFLAGS) -L$(LIBDIR) -I$(INCLUDEDIR) -lftr -lopenblas -o $@ $<

$(INCLUDEDIR)%.h : $(SRCDIR)%.h
	@mkdir -p $(INCLUDEDIR)
	cp $< $(INCLUDEDIR)
	
.PHONY: clock
clock:
	

clean:
	rm -f $(TARGET)
	rm -rf $(LIBDIR)
	rm -rf $(BUILDDIR)
	rm -rf $(INCLUDEDIR)

show:
	@echo $(PROGS)