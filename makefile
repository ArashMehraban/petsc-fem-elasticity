CFLAGS           =
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


all: main

main: main.o fe.o user.o exact.o chkopts
	-${CLINKER} -o main main.o fe.o user.o exact.o ${PETSC_KSP_LIB} ${PETSC_LIB} ${PETSC_SNES_LIB}
	${RM} main.o fe.o user.o exact.o




#myprogam: main.o foo.o
#    gcc -o myprogram main.o foo.o
#
#main.o: main.c foo.h
#    gcc -c main.c
#
#foo.o: foo.c foo.h
#    gcc -c foo.c
