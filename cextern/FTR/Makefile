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
BUILDDIR = build/
SRCDIR = src/

$(BUILDDIR)%.o : $(SRCDIR)%.c $(SRCDIR)%.h
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

%.so : $(BUILDDIR)%.o
	$(CC) $(LDFLAGS) -shared -o $@ $<

all: $(TARGET)

clean:
	rm -r $(BUILDDIR)
	rm $(TARGET)
