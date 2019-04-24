# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile,v 1.9 2003/10/21 17:17:51 hulsberg Exp $
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES := 
LIBTMPLFILES := 
LIBDDLORDERED := 

#################################################################
#++ extra binary products	[not in production, build it with extrabin]

EXTRABINS :=

$(addprefix $(bindir),$(EXTRABINS)): $(bindir)% : %.o

#################################################################
#++ binary products				[build it with 'bin']
BINS := HelixTestApp
BINCCFILES := AppHelixTest.cc

ifdef VTK_BOUNDSCHECKING
  CXXFLAGS += -DVTK_BOUNDSCHECKING
endif
#++ Binary rules		 [in production, build it with 'bin']

$(bindir)HelixTestApp: AppHelixTest.o

#++ shell script products.. 			[build it with 'bin']
BINSCRIPTS := 

#################################################################
#++ regression test scripts			[build it with 'test']

#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
