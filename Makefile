##############################################
#  NB: This is the main Makefile for PRISM.  #
#      It calls all the other Makefiles in   #
#      subdirectories, passing in all the    #
#      options configured here.              #
##############################################

####################
# Operating system #
####################

# OSTYPE needs to be one of: linux, solaris, cygwin, darwin
# This makefile will try to detect which one of these is appropriate.
# If this detection does not work, or you wish to override it,
# either uncomment one of the lines directly below
# or pass a value to make directly, e.g.: make OSTYPE=linux

#OSTYPE = linux
#OSTYPE = solaris
#OSTYPE = cygwin
#OSTYPE = darwin

ifdef OSTYPE
	# Look for common variants, e.g. gnu-linux -> linux
	ifneq (,$(findstring linux, $(OSTYPE)))
	  OSTYPE = linux
	endif
	ifneq (,$(findstring solaris, $(OSTYPE)))
	  OSTYPE = solaris
	endif
	ifneq (,$(findstring cygwin, $(OSTYPE)))
	  OSTYPE = cygwin
	endif
	ifneq (,$(findstring darwin, $(OSTYPE)))
	  OSTYPE = darwin
	endif
else
	# If OSTYPE is not defined/available, try uname
	ifneq (,$(findstring Linux, $(shell uname -s)))
		OSTYPE = linux
	endif
	ifneq (,$(findstring SunOS, $(shell uname -s)))
		OSTYPE = solaris
	endif
	ifneq (,$(findstring CYGWIN, $(shell uname -s)))
		OSTYPE = cygwin
	endif
	ifneq (,$(findstring Darwin, $(shell uname -s)))
		OSTYPE = darwin
	endif
endif

########
# Java #
########

# JAVA_DIR needs to be set to the location of your Java installation.
# This makefile will try to detect this automatically based on the location of the javac command.
# If this detection does not work, or you wish to override it,
# either set the variable yourself by uncommenting and/or modifying one of the lines directly below
# or pass a value to make directly, e.g.: make JAVA_DIR=/usr/java

#JAVA_DIR =	/usr/java
#JAVA_DIR =	/usr/java/j2sdk1.4.2
#JAVA_DIR =	/bham/java/packages/j2sdk1.4.2
#JAVA_DIR =	/cygdrive/c/java/j2sdk1.4.2
#JAVA_DIR =	/System/Library/Frameworks/JavaVM.framework

JAVA_JAVAC = $(shell JAVA_JAVAC=`which javac`; while [ -h $$JAVA_JAVAC ]; do JAVA_JAVAC=`/bin/ls -l $$JAVA_JAVAC | sed 's/.* -> //'`; done; echo $$JAVA_JAVAC)
ifneq (darwin,$(OSTYPE))
	JAVA_DIR =	$(shell dirname $(JAVA_JAVAC) | sed 's/\/bin//')
else
	JAVA_DIR =	$(shell dirname $(JAVA_JAVAC) | sed 's/\/Commands//')
endif

##################
# Compilers etc. #
##################

C = gcc
CPP = g++
LD = $(CPP)
JAVAC = javac
JAVAH = javah

##############
# Flags etc. #
##############

DEBUG = 
#DEBUG = -g

OPTIMISE = -O3
#OPTIMISE =

# Flags for compilation/linking
# Flags to generate shared libraries
# Executable/library naming conventions
# Suffix for binary distribution directory
# (requires GNU make for conditional evaluation)

# Linux
ifeq ($(OSTYPE),linux)
	CFLAGS = $(DEBUG) $(OPTIMISE)
	CPPFLAGS = $(DEBUG) $(OPTIMISE)
	LDFLAGS = $(DEBUG) $(OPTIMISE)
	SHARED = -shared
	#SHARED = -G
	EXE =
	LIBPREFIX = lib
	LIBSUFFIX = .so
	BINDISTSUFFIX = linux
endif
# Solaris
ifeq ($(OSTYPE),solaris)
	CFLAGS = $(DEBUG) $(OPTIMISE)
	CPPFLAGS = $(DEBUG) $(OPTIMISE)
	LDFLAGS = $(DEBUG) $(OPTIMISE)
	SHARED = -shared -mimpure-text
	EXE =
	LIBPREFIX = lib
	LIBSUFFIX = .so
	BINDISTSUFFIX = solaris
endif
# Cygwin
ifeq ($(OSTYPE),cygwin)
	CFLAGS = -mno-cygwin $(DEBUG) $(OPTIMISE)
	CPPFLAGS = -mno-cygwin $(DEBUG) $(OPTIMISE)
	LDFLAGS = -mno-cygwin -Wl,--add-stdcall $(DEBUG) $(OPTIMISE)
	SHARED = -shared
	#SHARED = -G
	EXE = .exe
	LIBPREFIX =
	LIBSUFFIX = .dll
	BINDISTSUFFIX = win
endif
# Darwin
ifeq ($(OSTYPE),darwin)
	CFLAGS = $(DEBUG) $(OPTIMISE)
	CPPFLAGS = $(DEBUG) $(OPTIMISE)
	LDFLAGS = $(DEBUG) $(OPTIMISE)
	SHARED = -dynamiclib
	EXE =
	LIBPREFIX = lib
	LIBSUFFIX = .dylib
	BINDISTSUFFIX = osx
endif

#APMC = no
#APMC = yes

###############
# Directories #
###############

# Note that these are all relative to the CUDD directory
# to make the distribution more 'portable'.
# If this is a problem, the best solution is to create symlinks.

CUDD_DIR =		cudd
SRC_DIR =		src
CLASSES_DIR =	classes
OBJ_DIR =		obj
LIB_DIR =		lib
INCLUDE_DIR =	include

# Now we locate the JNI header files jni.h and jni_md.h
# (this is the only reason we need JAVA_DIR).
# If this doesn't work for some reason, locate the two files
# and set the JAVA_INCLUDES directory manually

# Choose name for include directories
ifneq (darwin,$(OSTYPE))
	OSTYPE_INCLUDE = include
else
	OSTYPE_INCLUDE = Headers
endif
# Look in Java directory
JAVA_JNI_H = $(shell ls $(JAVA_DIR)/$(OSTYPE_INCLUDE)/jni.h 2>/dev/null)
JAVA_JNI_MD_H = $(shell ls $(JAVA_DIR)/$(OSTYPE_INCLUDE)/jni_md.h 2>/dev/null)
# Then try subdirectories (for OS-specific files)
ifeq (,$(JAVA_JNI_MD_H))
	JAVA_JNI_MD_H = $(shell ls $(JAVA_DIR)/$(OSTYPE_INCLUDE)/*/jni_md.h 2>/dev/null)
endif
# Now strip off filename to leave directories
JAVA_JNI_H_DIR = $(shell echo $(JAVA_JNI_H) | sed 's/\/jni.h//')
JAVA_JNI_MD_H_DIR = $(shell echo $(JAVA_JNI_MD_H) | sed 's/\/jni_md.h//')
# Store result in JAVA_INCLUDES variable
JAVA_INCLUDES = -I $(JAVA_JNI_H_DIR) -I $(JAVA_JNI_MD_H_DIR)

#########################
# Main part of Makefile #
#########################

MAKE_DIRS = dd jdd odd dv prism mtbdd sparse hybrid parser settings chart userinterface pepa/compiler apmc simulator

default: all

all: checks cuddpackage prism

cuddpackage:
	echo Making cudd ...; \
	cd cudd && \
	/bin/cp Makefile.$(OSTYPE) Makefile && \
	$(MAKE)

checks:
	@(if [ "$(OSTYPE)" != "linux" -a "$(OSTYPE)" != "solaris" -a "$(OSTYPE)" != "cygwin" -a "$(OSTYPE)" != "darwin" ]; then \
	  cat .ostype.txt; \
	  exit 1; \
	fi; \
	if [ "$(JAVA_DIR)" = "" ]; then \
	  cat .java_dir.txt; \
	  exit 1; \
	fi; \
	if [ ! -f $(JAVA_JNI_H_DIR)/jni.h ]; then \
	  echo "Failed to locate JNI headers jni.h and jni_md.h. Are you sure Java is installed?"; \
	  exit 1; \
	fi)

prism: checks sortplugins make_dirs post_make

sortplugins:
# 	@(rm -f $(SRC_DIR)/apmc; \
# 	if [ "$(APMC)" = "yes" ]; then \
# 	  (ln -s plugins/apmc $(SRC_DIR)/apmc) \
# 	else \
# 	  (ln -s stubs/apmc $(SRC_DIR)/apmc) \
# 	fi)

make_dirs:
	@for dir in $(MAKE_DIRS); do \
	  echo Making src/$$dir ...; \
	  (cd src/$$dir && \
	  $(MAKE) \
	  CUDD_DIR="$(CUDD_DIR)" \
	  SRC_DIR="$(SRC_DIR)" \
	  CLASSES_DIR="$(CLASSES_DIR)" \
	  OBJ_DIR="$(OBJ_DIR)" \
	  LIB_DIR="$(LIB_DIR)" \
	  INCLUDE_DIR="$(INCLUDE_DIR)" \
	  JAVA_INCLUDES="$(JAVA_INCLUDES)" \
	  C="$(C)" \
	  CPP="$(CPP)" \
	  LD="$(LD)" \
	  JAVAC="$(JAVAC)" \
	  JAVAH="$(JAVAH)" \
	  CFLAGS="$(CFLAGS)" \
	  CPPFLAGS="$(CPPFLAGS)" \
	  LDFLAGS="$(LDFLAGS)" \
	  SHARED="$(SHARED)" \
	  EXE="$(EXE)" \
	  LIBPREFIX="$(LIBPREFIX)" \
	  LIBSUFFIX="$(LIBSUFFIX)") \
	  || exit 1; \
	done

post_make:
	@(if [ "$(OSTYPE)" = "darwin" ]; then \
	  echo Creating shared library symlinks...; \
	  (cd $(LIB_DIR) && \
	  for lib in `ls *$(LIBSUFFIX)`; do ln -fs $$lib `echo $$lib | sed s/$(LIBSUFFIX)/.jnilib/`; done;); \
	fi; \
	echo Fixing startup scripts...; \
	(cd bin && \
	if [ "$(OSTYPE)" = "darwin" ]; then \
	  sed 's/\(DY\)*LD_LIBRARY_PATH/DYLD_LIBRARY_PATH/g' prism > prism.tmp; \
	else \
	  sed 's/\(DY\)*LD_LIBRARY_PATH/LD_LIBRARY_PATH/g' prism > prism.tmp; \
	fi; \
	mv prism.tmp prism; \
	chmod 755 prism ))

dist: dist_copy clean_all dist_files dist_tidy

dist_bin: binary dist_tidy dist_bin_copy

dist_copy:
	@echo Adding cudd...; test -d ~dxp/public/cudd-linux && (cd ~dxp/public && tar cf - cudd-linux) | tar xf - && rm -rf cudd && mv cudd-linux cudd
	@echo Adding doc...; test -f ~dxp/prism-dev/doc/manual.pdf && rm -rf doc && mkdir doc && cp ~dxp/prism-dev/doc/manual.pdf doc
	@echo Adding examples...; test -d ~dxp/prism-examples && rm -rf examples prism-examples && cp -r ~dxp/prism-examples . && mv prism-examples examples

dist_bin_copy:
	@BIN_DIST_DIR=`/bin/pwd | sed 's/-src$$//'`"-$(BINDISTSUFFIX)" && \
	echo Creating binary distribution in $$BIN_DIST_DIR... && \
	mkdir $$BIN_DIST_DIR && \
	tar cf - * --exclude classes --exclude obj --exclude cudd --exclude src --exclude include --exclude Makefile | ( cd $$BIN_DIST_DIR; tar xfp -) && \
	if [ "$(BINDISTSUFFIX)" = "win" ]; then \
		zip -rq $$BIN_DIST_DIR.zip $$BIN_DIST_DIR; \
	else \
		tar cfz $$BIN_DIST_DIR.tar.gz $$BIN_DIST_DIR; \
	fi

dist_files:
	@echo Detecting unwanted files...
	@find . \( -name '*.o' -o -name '*.so' -o -name '*.dll' -o -name '*.exe' \)
	@find . \( -name '*~*' -o -name '*bak*' -o -name '*CVS*' \)
	@find . \( -name '*NEW*' -o -name '*new*' -o -name '*old*' \) | grep -v old\.gif || test 1
	@find . -name '*tmp*' | grep -v cudd/util/tmpfile.c || test 1
	@find . -name 'log*' | grep -v userinterface/log || test 1
	@find . -name '.nbattrs'
	@find . -name '.xv*'
	@find . -name '*NOTES*' | grep -v src/parser/NOTES | grep -v cudd/RELEASE.NOTES || test 1

dist_tidy:
	@echo Processing text files...
	@unix2dos *.txt 2> /dev/null
	@find examples -type f -not -name auto -exec unix2dos {} \; 2> /dev/null
	@echo Processing file permissions...
	@find . -type f -exec chmod 644 {} \;
	@find . \( -type d -o -type s \) -exec chmod 755 {} \;
	@find . -type f  \( -name '*.sh' -o -name '*.so' -o -name '*.dll' \) -exec chmod 755 {} \;
	@find examples -type f -name 'auto' -exec chmod 755 {} \;
	@chmod 755 bin/*

binary:
	@echo Generating jar file...
	@cd $(CLASSES_DIR) && jar cf ../lib/prism.jar *

undist:
	@rm -rf cudd && ln -s ../cudd cudd
	@rm -rf doc
	@rm -rf examples && ln -s ../prism-examples examples

tarcf:
	@TARCF_DIR=`/bin/pwd | sed 's/.\+\///'` && \
	if [ $$TARCF_DIR = "." ]; then exit 1; fi && \
	echo Building tar file "../"$$TARCF_DIR".tar.gz" && \
	(cd ..; tar cfz $$TARCF_DIR".tar.gz" $$TARCF_DIR)

javadoc:
	javadoc -classpath $(SRC_DIR) -d ../prism-dev/javadoc dd jdd odd dv mtbdd sparse hybrid parser prism chart userinterface apmc

clean: checks
	@(for dir in $(MAKE_DIRS); do \
	  echo Cleaning src/$$dir ...; \
	  (cd src/$$dir && \
	  $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean) \
	  || exit 1; \
	done; \
	find $(CLASSES_DIR) -name '*.class' -exec rm {} \; ; \
	rm -f lib/*jnilib)

celan: clean

clean_all: checks clean_cudd clean

clean_cudd:
	@(if [ ! -h cudd ]; then \
	  cd cudd && $(MAKE) distclean; \
	fi)

clean_dd: checks
	@(cd src/dd && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_jdd: checks
	@(cd src/jdd && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_odd: checks
	@(cd src/odd && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_dv: checks
	@(cd src/dv && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_prism: checks
	@(cd src/prism && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_mtbdd: checks
	@(cd src/mtbdd && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_sparse: checks
	@(cd src/sparse && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_hybrid: checks
	@(cd src/hybrid && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_parser: checks
	@(cd src/parser && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_chart: checks
	@(cd src/chart && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_userinterface: checks
	@(cd src/userinterface && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_apmc: checks
	@(cd src/apmc && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)
clean_simulator: checks
	@(cd src/simulator && $(MAKE) -s SRC_DIR="$(SRC_DIR)" CLASSES_DIR="$(CLASSES_DIR)" OBJ_DIR="$(OBJ_DIR)" LIB_DIR="$(LIB_DIR)" EXE="$(EXE)" LIBPREFIX="$(LIBPREFIX)" LIBSUFFIX="$(LIBSUFFIX)" clean)


#################################################
