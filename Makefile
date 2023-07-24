# Directories
OBJDIR := obj
SRCDIR := src
OUTDIR := output
INCDIR := headers 
# Source directories
SRCDIR_TRACER := $(SRCDIR)/tracer
SRCDIR_INTERP := $(SRCDIR)/interpolation
SRCDIR_DSTR   := $(SRCDIR)/distribution
# Object directories
OBJDIR_TRACER := $(OBJDIR)/tracer
OBJDIR_INTERP := $(OBJDIR)/interpolation
OBJDIR_DSTR   := $(OBJDIR)/distribution
# Output directories
OFILESDIR     := $(OUTDIR)/files
PLOTSDIR 	  := $(OUTDIR)/plots

# Targets
BIN := tracer ray dstr

# Compiler and linker
CC  := mpic++
LNK := mpic++

# Optimization level 
OPT := -O3 

# Include directories and flags
INC_H5 := -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I HighFive/include/
INC_H5MIN := -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial  -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I HighFive/include/

# Compiler flags
CPPFLAGS := -std=c++17 -lstdc++fs -Wall -g -fopenmp -I$(INCDIR)
debug: CPPFLAGS += -DDEBUG

SRC_TRACER := $(wildcard ${SRCDIR_TRACER}/*.cc)
SRC_INTERP := $(wildcard ${SRCDIR_INTERP}/*.cc)
SRC_DSTR   := $(wildcard ${SRCDIR_DSTR}/*.cc)

OBJ_TRACER := $(patsubst ${SRCDIR_TRACER}/%.cc, ${OBJDIR_TRACER}/%.o, $(SRC_TRACER))
OBJ_INTERP := $(patsubst ${SRCDIR_INTERP}/%.cc, ${OBJDIR_INTERP}/%.o, $(SRC_INTERP))
OBJ_DSTR   := $(patsubst ${SRCDIR_DSTR}/%.cc, ${OBJDIR_DSTR}/%.o, $(SRC_DSTR))


# Default Rule
default: directories dstr ray tracer
	@printf "\nEverything is Built\n"

# Rules
debug: default
	@printf "\nEverything is built in DEBUG mode\n"

ray: ${OBJ_INTERP}
	${LNK} -o $@ ${OBJ_INTERP} ${INC_H5} ${CPPFLAGS} 
	@printf "\nBuilt to interpolate ray\n\n" ;
dstr: ${OBJ_DSTR}
	${LNK} -o $@ ${OBJ_DSTR} ${INC_H5} ${CPPFLAGS} 
	@printf "\nBuilt to distribute particles\n\n" ;
tracer: ${OBJ_TRACER}
	${LNK} -o $@ ${OBJ_TRACER} ${INC_H5} ${CPPFLAGS}
	@printf "\nBuilt to trace particles\n" ;
		
# Phony Targets(Just recipes)
directories:
	@mkdir -p ${OFILESDIR} ${PLOTSDIR} ${OBJDIR_TRACER} ${OBJDIR_INTERP} ${OBJDIR_DSTR}
	@printf "Directories are made\n\n" ;
clean: 
	@rm  -rf ${OBJDIR} ${BIN}
	@echo Removed objects and executables ;
allclean:
	@rm -rf **.out ${BIN} ${OBJDIR} ${OFILESDIR} ${PLOTSDIR}
	@echo Removed all. Restored to default ;

.PHONY: directories clean allclean


# Compiling

# Tracer
${OBJDIR_TRACER}/%.o: ${SRCDIR_TRACER}/%.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5MIN}

# Distribution

${OBJDIR_DSTR}/functions.o: ${SRCDIR_TRACER}/functions.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5MIN}
${OBJDIR_DSTR}/struct_Particles.o: ${SRCDIR_TRACER}/struct_Particles.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5MIN}
${OBJDIR_DSTR}/distribution.o: ${SRCDIR_DSTR}/distribution.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5MIN}


# Interpolation
$(OBJDIR_INTERP)/%.o: $(SRCDIR_INTERP)/%.cc
	$(CC) $(CPPFLAGS) $(OPT) -c $< -o $@ ${INC_H5MIN}

