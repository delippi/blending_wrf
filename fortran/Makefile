SHELL = /bin/sh
.SUFFIXES: .f .f90 .o

ifneq ("$(wildcard $(NETCDF)/lib/libnetcdff*)", "")
NETCDF_LIB=-L$(NETCDF)/lib -lnetcdff -lnetcdf -lcurl
else
NETCDF_LIB=-L$(NETCDF)/lib -lnetcdf
endif

UNAME    := $(shell uname)

F90_AIX      = xlf90
F77_AIX      = xlf
FFLAGS_AIX   = -qintsize=4 -qrealsize=8

F90_Linux    = pgf90
F77_Linux    = pgf77
FFLAGS_Linux = -r8 -i4

#F90    = $(F90_$(UNAME))
#F77    = $(F77_$(UNAME))
#FFLAGS = $(FFLAGS_$(UNAME))

#F90    = mpifrtpx
#F77    = mpifrtpx
F90    = ifort
F77    = ifort
FFLAGS = -r8 # #-fw -CcdRR8 -Kfast,parallel  -Qa,d,i,p,t,x -Koptmsg=2
#FFLAGS = -r8 -g -traceback -warn # #-fw -CcdRR8 -Kfast,parallel  -Qa,d,i,p,t,x -Koptmsg=2

LDFLAG = "" #-Kparallel

OBJ =  da_blending.o raymond.o

OBJ1=  da_diagnose.o balance.o

all:da_blending.exe da_diagnose.exe

da_blending.exe: $(OBJ)
	$(F90) $(LDFLAG) -o $@  $(OBJ) $(NETCDF_LIB)

da_diagnose.exe: $(OBJ1)
	$(F90) $(LDFLAG) -o $@  $(OBJ1) $(NETCDF_LIB)

.f90.o:
	$(F90) -c $(FFLAGS) -I$(NETCDF)/include $<

.f.o:
	$(F77) -c $(FFLAGS) $<

clean:
	rm *.o *.exe

