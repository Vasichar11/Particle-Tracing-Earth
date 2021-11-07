BIN = tracer ray_int
CC  = g++ 
LNK = g++ 
OPT = -O3
CXXFLAGS = -std=c++17 -lstdc++fs -Wall -pg -fopenmp


INC_HIFIVE = -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/billy/HighFive/include/

INC_H5LESS = -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial  -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/billy/HighFive/include/

OBJ_TRACER =          main.o  bell_wpi.o   bell_nowpi.o  li_wpi.o   functions.o   time_rates.o   struct_Particles.o    struct_Telescope.o   struct_Species.o  is_in_packet.o    read_hdf5.o

OBJ_INT = Read_Ray_Write.o   read_csv.o   interpolate.o



#RULES
default:tracer
all:${BIN}
clean: 
	rm  *.o 
allclean:
	rm *.o ${BIN} *.out *.png


#For all code
struct_Particles.o: struct_Particles.cc
	${CC} ${CXXFLAGS} ${OPT} -c struct_Particles.cc   
struct_Species.o: struct_Species.cc 	#Only for Bell code
	${CC} ${CXXFLAGS} ${OPT} -c struct_Species.cc 
struct_Telescope.o: struct_Telescope.cc
	${CC} ${CXXFLAGS} ${OPT} -c struct_Telescope.cc    
functions.o: functions.cc
	${CC} ${CXXFLAGS} ${OPT} -c functions.cc
time_rates.o: time_rates.cc
	${CC} ${CXXFLAGS} ${OPT} -c time_rates.cc


#For Ray Interpolation
read_hdf5.o: read_hdf5.cc
	${CC} ${CXXFLAGS} ${OPT} -c read_hdf5.cc  ${INC_H5LESS}
read_csv.o: read_csv.cc
	${CC} ${CXXFLAGS} ${OPT} -c read_csv.cc
interpolate.o: interpolate.cc
	${CC} ${CXXFLAGS} ${OPT} -c interpolate.cc
Read_Ray_Write.o: Read_Ray_Write.cc
	${CC} ${CXXFLAGS} ${OPT} -c Read_Ray_Write.cc  ${INC_H5LESS}



bell_nowpi.o: bell_nowpi.cc
	${CC} ${CXXFLAGS} ${OPT} -c bell_nowpi.cc
bell_wpi.o: bell_wpi.cc
	${CC} ${CXXFLAGS} ${OPT} -c bell_wpi.cc
li_wpi.o: li_wpi.cc
	${CC} ${CXXFLAGS} ${OPT} -c li_wpi.cc ${INC_H5LESS}


main.o: main.cc
	${CC} ${CXXFLAGS} ${OPT} -c main.cc ${INC_H5LESS}










ray_int:  ${OBJ_INT}
	${LNK} -o ray_int ${OBJ_INT} ${INC_HIFIVE} ${CXXFLAGS}

tracer: ${OBJ_TRACER}
	${LNK} -o tracer ${OBJ_TRACER} ${INC_HIFIVE} ${CXXFLAGS} ${OMP}
