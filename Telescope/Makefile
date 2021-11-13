BIN := tracer ray_int
CC  := g++
LNK := g++ 
OPT := -O3
INC_H5 := -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/billy/HighFive/include/
INC_H5LESS := -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial  -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/billy/HighFive/include/
CPPFLAGS := -std=c++17 -lstdc++fs -Wall -g -fopenmp 

SRCDIR_TRACER := src/tracer
SRCDIR_INTERP := src/interpolation
OBJDIR_TRACER := obj/tracer
OBJDIR_INTERP := obj/interpolation
OBJDIR := obj
H5DIR  := h5files
SIMDIR := simulation_MM

SRC_TRACER := $(wildcard ${SRCDIR_TRACER}/*.cc)
SRC_INTERP := $(wildcard ${SRCDIR_INTERP}/*.cc)

OBJ_TRACER := $(patsubst ${SRCDIR_TRACER}/%.cc, ${OBJDIR_TRACER}/%.o, $(SRC_TRACER))
OBJ_INTERP := $(patsubst ${SRCDIR_INTERP}/%.cc, ${OBJDIR_INTERP}/%.o, $(SRC_INTERP))


#Default Rule
default: directories ray_int tracer
	@printf "\nEverything is builded\n"

#Implicit Rules
ray_int: ${OBJ_INTERP}
	${LNK} -o $@ ${OBJ_INTERP} ${INC_H5} ${CPPFLAGS} 
	@printf "Builded to interpolate ray\n\n" ;
tracer: ${OBJ_TRACER}
	${LNK} -o $@ ${OBJ_TRACER} ${INC_H5} ${CPPFLAGS}
	@printf "\nBuilded to trace particles\n" ;
	
#Phony Targets(Just recipes)
directories:
	@mkdir -p ${H5DIR} ${SIMDIR} ${OBJDIR_TRACER} ${OBJDIR_INTERP}
	@printf "Directories are made\n\n" ;
clean: 
	@rm  -rf ${OBJDIR} ${BIN}
	@echo Removed objects and executables ;
allclean:
	@rm -rf **.out ${BIN} ${OBJDIR} ${H5DIR} ${SIMDIR}
	@echo Removed all. Restored to default ;

.PHONY: clean allclean



#Compiling Objects 
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
${OBJDIR_TRACER}/read_hdf5.o: ${SRCDIR_TRACER}/read_hdf5.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5LESS}
${OBJDIR_TRACER}/bell_nowpi.o: ${SRCDIR_TRACER}/bell_nowpi.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@
${OBJDIR_TRACER}/bell_wpi.o: ${SRCDIR_TRACER}/bell_wpi.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@
${OBJDIR_TRACER}/li_wpi.o: ${SRCDIR_TRACER}/li_wpi.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5LESS}
${OBJDIR_TRACER}/main.o: ${SRCDIR_TRACER}/main.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5LESS}

${OBJDIR_INTERP}/read_csv.o: ${SRCDIR_INTERP}/read_csv.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@
${OBJDIR_INTERP}/interpolate.o: ${SRCDIR_INTERP}/interpolate.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@
${OBJDIR_INTERP}/Read_Ray_Write.o: ${SRCDIR_INTERP}/Read_Ray_Write.cc
	${CC} ${CPPFLAGS} ${OPT} -c $< -o $@ ${INC_H5LESS}



