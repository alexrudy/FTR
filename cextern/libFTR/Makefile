#
#  Makefile
#  FTR
#
#  Created by Alexander Rudy on 2015-08-13.
#  Copyright 2015 Alexander Rudy. All rights reserved.
#

LDFLAGS = -g -O3 -Wall -lfftw3_threads -lfftw3 -lm -lpthread -L /opt/local/lib
CFLAGS = -g -O3 -Wall -I$(INCLUDEDIR)
LIBFLAGS = $(LDFLAGS) -shared -fPIC
BUILDDIR = build/
SRCDIR = src/
LIBDIR = lib/
EXAMPLE = example/
INCLUDEDIR = include/
HEADERS = $(wildcard $(SRCDIR)*.h)
INCLUDE = $(patsubst $(SRCDIR)%.h,$(INCLUDEDIR)%.h, $(HEADERS))

TARGET = $(LIBDIR)libftr.so $(INCLUDE)
OBJECTS = $(patsubst $(SRCDIR)%.c,$(BUILDDIR)%.o, $(wildcard $(SRCDIR)*.c))

PROGS_C = $(wildcard $(EXAMPLE)*.c)
PROGS = $(subst _,-, $(patsubst $(EXAMPLE)%.c,%,$(PROGS_C)))

all: $(TARGET)
	
tests: $(PROGS)

$(BUILDDIR)%.o: $(SRCDIR)%.c $(INCLUDE)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(EXAMPLE)%.o: $(EXAMPLE)%.c $(INCLUDE)

$(LIBDIR)lib%.so: $(OBJECTS)
	@mkdir -p $(LIBDIR)
	$(CC) $(LIBFLAGS) -o $@ $^


ftr-example: LDFLAGS += -L$(LIBDIR) -I$(INCLUDEDIR) -lftr
ftr-example: $(EXAMPLE)ftr_example.o $(LIBDIR)libftr.so
	$(CC) $(LDFLAGS) -o $@ $<

recon-timing:  LDFLAGS += -L$(LIBDIR) -I$(INCLUDEDIR) -lftr -lopenblas
recon-timing: $(EXAMPLE)recon_timing.o $(LIBDIR)libftr.so
	$(CC) $(LDFLAGS) -o $@ $<

aperture-test:  LDFLAGS += -L$(LIBDIR) -I$(INCLUDEDIR) -lftr
aperture-test: $(EXAMPLE)aperture_test.o $(LIBDIR)libftr.so
	$(CC) $(LDFLAGS) -o $@ $<

$(INCLUDEDIR)%.h: $(SRCDIR)%.h
	@mkdir -p $(INCLUDEDIR)
	cp $< $(INCLUDEDIR)

.PHONY: clean all
	
clean:
	rm -f $(TARGET)
	rm -rf $(LIBDIR)
	rm -rf $(BUILDDIR)
	rm -rf $(INCLUDEDIR)