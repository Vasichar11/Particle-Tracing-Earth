BIN := tracer ray dstr
CC  := g++
LNK := g++
OPT := -O3 

#Ubuntu includes
INC_H5 := -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/vchariton/HighFive/include/
INC_H5LESS := -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial  -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/vchariton/HighFive/include/

#Centos includes
#INC_H5 := -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security  /usr/lib64/libhdf5.so.103 /usr/lib64/libhdf5_cpp.so.103 /usr/lib64/libhdf5_hl.so.100 /usr/lib64/libhdf5_hl_cpp.so.100 -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl, -I/home/bill-centos/HighFive/include/
#INC_H5LESS := -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security   -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl, -I/home/bill-centos/HighFive/include/



#Profiling with gprof
#add the -g -pg flag both for compilation and linking nad write to gmon.out. Then invoke gprof with the executable


CPPFLAGS := -std=c++17 -lstdc++fs -Wall -g -fopenmp 
debug: CPPFLAGS += -DDEBUG

SRCDIR_TRACER := src/tracer
SRCDIR_INTERP := src/interpolation
SRCDIR_DSTR   := src/distribution
OBJDIR_TRACER := obj/tracer
OBJDIR_INTERP := obj/interpolation
OBJDIR_DSTR   := obj/distribution
OBJDIR := obj
H5DIR  := h5files
SIMDIR := simulation_MM

SRC_TRACER := $(wildcard ${SRCDIR_TRACER}/*.cc)
SRC_INTERP := $(wildcard ${SRCDIR_INTERP}/*.cc)

OBJ_TRACER := $(patsubst ${SRCDIR_TRACER}/%.cc, ${OBJDIR_TRACER}/%.o, $(SRC_TRACER))
OBJ_INTERP := $(patsubst ${SRCDIR_INTERP}/%.cc, ${OBJDIR_INTERP}/%.o, $(SRC_INTERP))

SRC_DSTR := $(wildcard ${SRCDIR_DSTR}/*.cc)
OBJ_DSTR := $(patsubst ${SRCDIR_DSTR}/%.cc, ${OBJDIR_DSTR}/%.o, $(SRC_DSTR))



#Default Rule
default: directories dstr ray tracer
	@printf "\nEverything is builded\n"

debug: default
	@printf "\nEverything is built in DEBUG mode\n"

#Implicit Rules

ray: ${OBJ_INTERP}
	${LNK} -o $@ ${OBJ_INTERP} ${INC_H5} ${CPPFLAGS} 
	@printf "\nBuilded to interpolate ray\n\n" ;
dstr: ${OBJ_DSTR}
	${LNK} -o $@ ${OBJ_DSTR} ${INC_H5} ${CPPFLAGS} 
	@printf "\nBuilded to distribute particles\n\n" ;
tracer: ${OBJ_TRACER}
	${LNK} -o $@ ${OBJ_TRACER} ${INC_H5} ${CPPFLAGS}
	@printf "\nBuilded to trace particles\n" ;
	
#Phony Targets(Just recipes)
directories:
	@mkdir -p ${H5DIR} ${SIMDIR} ${OBJDIR_TRACER} ${OBJDIR_INTERP} ${OBJDIR_DSTR}
	@printf "Directories are made\n\n" ;
clean: 
	@rm  -rf ${OBJDIR} ${BIN}
	@echo Removed objects and executables ;
#Comment to avoid removing directories by mistake.
allclean:
	@rm -rf **.out ${BIN} ${OBJDIR} ${H5DIR} ${SIMDIR}
	@echo Removed all. Restored to default ;

.PHONY: clean allclean



#COMPILING

#TRACER
${OBJDIR_TRACER}/struct_Particles.o: ${SRCDIR_TRACER}/struct_Particles.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@   
${OBJDIR_TRACER}/struct_Species.o: ${SRCDIR_TRACER}/struct_Species.cc 	
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ 
${OBJDIR_TRACER}/struct_Telescope.o: ${SRCDIR_TRACER}/struct_Telescope.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@    
${OBJDIR_TRACER}/functions.o: ${SRCDIR_TRACER}/functions.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@
${OBJDIR_TRACER}/is_in_packet.o: ${SRCDIR_TRACER}/is_in_packet.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ 
${OBJDIR_TRACER}/time_rates.o: ${SRCDIR_TRACER}/time_rates.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@
${OBJDIR_TRACER}/RK_estimations.o: ${SRCDIR_TRACER}/RK_estimations.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@
${OBJDIR_TRACER}/read_hdf5.o: ${SRCDIR_TRACER}/read_hdf5.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5LESS}
${OBJDIR_TRACER}/no_wpi.o: ${SRCDIR_TRACER}/no_wpi.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@
${OBJDIR_TRACER}/bell_wpi.o: ${SRCDIR_TRACER}/bell_wpi.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@
${OBJDIR_TRACER}/li_wpi.o: ${SRCDIR_TRACER}/li_wpi.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5LESS}
${OBJDIR_TRACER}/main.o: ${SRCDIR_TRACER}/main.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5LESS}

#DISTRIBUTION, WITHOUT THE SYMBOLIC LINKS FOR SOURCE CODES struct_Particles.cc, functions.cc
${OBJDIR_DSTR}/distribution.o: ${SRCDIR_DSTR}/distribution.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5LESS}
${OBJDIR_DSTR}/struct_Particles.o: ${SRCDIR_TRACER}/struct_Particles.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5LESS}
${OBJDIR_DSTR}/functions.o: ${SRCDIR_TRACER}/functions.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5LESS}

#INTERPOLATION
${OBJDIR_INTERP}/read_csv.o: ${SRCDIR_INTERP}/read_csv.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@
${OBJDIR_INTERP}/interpolate.o: ${SRCDIR_INTERP}/interpolate.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@
${OBJDIR_INTERP}/Read_Ray_Write.o: ${SRCDIR_INTERP}/Read_Ray_Write.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5LESS}


