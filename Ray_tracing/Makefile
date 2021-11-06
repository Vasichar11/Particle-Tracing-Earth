BIN=tracer ray_int
CC=g++ 
LNK=g++  
OPT=-O3
CXXFLAGS=-std=c++17 -lstdc++fs -Wall -g -fopenmp


INC_HIFIVE= -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/billy/HighFive/include/


INC_H5FAST= -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/billy/HighFive/include/ -lhdf5


OBJ= Li_Ray_Telescope.o  struct_Particles.o  struct_Telescope.o  wpi_ray.o  read_hdf5.o  functions.o  time_rates.o  is_in_packet.o 



default:tracer

all:ray_int tracer

clean: 
	rm  *.o ${BIN}

allclean:
	rm *.o ${BIN} */*.h5



read_hdf5.o: read_hdf5.cc
	${CC} ${CXXFLAGS} ${OPT} -c read_hdf5.cc  ${INC_H5FAST}

struct_Particles.o: struct_Particles.cc
	${CC} ${CXXFLAGS} ${OPT} -c struct_Particles.cc   

struct_Telescope.o: struct_Telescope.cc
	${CC} ${CXXFLAGS} ${OPT} -c struct_Telescope.cc    ${INC_FAST}

functions.o: functions.cc
	${CC} ${CXXFLAGS} ${OPT} -c functions.cc

time_rates.o: time_rates.cc
	${CC} ${CXXFLAGS} ${OPT} -c time_rates.cc

wpi_ray.o: wpi_ray.cc
	${CC} ${CXXFLAGS} ${OPT} -c wpi_ray.cc ${INC_H5FAST}

Li_Ray_Telescope.o: Li_Ray_Telescope.cc
	${CC} ${CXXFLAGS} ${OPT} -c Li_Ray_Telescope.cc  ${INC_H5FAST}

read_csv.o: read_csv.cc
	${CC} ${CXXFLAGS} ${OPT} -c read_csv.cc

interpolate.o: interpolate.cc
	${CC} ${CXXFLAGS} ${OPT} -c interpolate.cc

Read_Ray_Write.o: Read_Ray_Write.cc
	${CC} ${CXXFLAGS} ${OPT} -c Read_Ray_Write.cc  ${INC_H5FAST}


ray_int:  Read_Ray_Write.o read_csv.o interpolate.o
	${LNK} -o ray_int Read_Ray_Write.o read_csv.o interpolate.o ${INC_H5FAST} ${CXXFLAGS}

tracer: ${OBJ}
	${LNK} -o tracer ${OBJ} ${INC_H5FAST} ${CXXFLAGS} ${OMP}

