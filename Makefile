CC  	:= gcc -arch x86_64

CFLAGS 	:= -Wall -O2 -dynamic
LIBCFLAGS := -Wall -O3 \
		-current_version 1.0 \
		-compatibility_version 1.0 \
		-fvisibility=hidden \
		-dynamiclib \
		-std=gnu99

INCINSTALLDIR := /usr/local/include/snopt
LIBINSTALLDIR := /usr/local/lib
VERSION := 1

SNOPT_INCLUDE   := /usr/local/include/snopt
SNOPT_LIBDIR	:= /usr/local/lib
SNOPT_LIBS  	:= -lsnopt_c -lsnprint_c -lblas_c
SNOPT_AR     	:= libsnopt_c libsnprint_c libblas_c
SNOPT_AR_LIBS 	:= $(SNOPT_AR:%=$(SNOPT_LIBDIR)/%.a)
SNOPT_EXTRAS	:= libsnoextras
SNOPT_EX_LIBS	:= $(SNOPT_EXTRAS:%=$(SNOPT_LIBDIR)/%.a)

F2C_INCLUDE     := /usr/local/include/snopt
F2C_LIBDIR	:= /usr/local/lib

BLASLIBS	:= -framework Accelerate

default:

	# compiling object files
	$(CC) -c snoextras.c -o snoextras.o $(CFLAGS) \
		-I$(SNOPT_INCLUDE) \
		# -L/usr/local/lib $(BLASLIBS) $(SNOPT_AR_LIBS) $(F2C_LIBDIR)/libf2c.a -lm

	# creating static library
	libtool snoextras.o -o libsnoextras.a -static -s -v
	
	# creating dynamic library
	libtool snoextras.o -o libsnoextras.dylib \
		-L/usr/local/lib $(BLASLIBS) $(SNOPT_AR_LIBS) $(F2C_LIBDIR)/libf2c.a -lm \
		-dynamic \
		# -install_name $(LIBINSTALLDIR)/ \
		-current_version $(VERSION) \
		-v

tests: 
	
	$(CC) $(CFLAGS) tests.c -o tests \
		$(INCLUDE) -I$(SNOPT_INCLUDE) \
		-L/usr/local/lib $(BLASLIBS) $(SNOPT_AR_LIBS) $(SNOPT_EX_LIBS) $(F2C_LIBDIR)/libf2c.a -lm

install:
	# copying files to their installed location
	cp snoextras.h $(INCINSTALLDIR)/
	cp libsnoextras* $(LIBINSTALLDIR)/

uninstall:

	rm $(INCINSTALLDIR)/snoextras.h
	rm $(LIBINSTALLDIR)/libsnoextras*

clean: 
	rm snoextras.o
	rm libsnoextras*
	rm tests

# Fake target to remind people to set the F2C environment variable with SNOPT

$(F2C_INCLUDE)/f2c.h:
	@echo "Could not find the f2c distribution."
	@echo "Set the following environment variables:"
	@echo "  F2CINCLUDE should be the path to f2c.h"
	@false