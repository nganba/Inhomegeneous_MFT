FC = ifort


#My Laptop
#LDFLAGS =  /usr/lib/lapack/liblapack.so.3gf  /usr/lib/libblas/libblas.so.3gf
#Fermion
#LDFLAGS =  -L/usr/local/mkl-8.1.1/lib/32 -lmkl -lmkl_lapack -lguide -pthread
#Harmony
#LDFLAGS =  -L/usr/local/mkl-8.1.1/lib/em64t -lmkl -lmkl_lapack -lguide -pthread
#glenn
#LDFLAGS = ${ACML}
#oakley
LDFLAGS = -L/usr/local/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64 -mkl -lmkl_sequential
#newton
#LDFLAGS = /home/nganba/lapack-3.3.1/lapack_LINUX.a  /home/nganba/lapack-3.3.1/blas_LINUX.a

main = ./subroutines

OBJ1 = $(main)/stdtypes.f90   \
$(main)/mtprng.f90   \
$(main)/hamiltonian.f90  \
$(main)/HF.f90 

OBJ2 = $(main)/stdtypes.f90   \
$(main)/mtprng.f90   \
$(main)/hamiltonian.f90  \
$(main)/magnetic_basis.f90  \
$(main)/generate_k_points_withTBC.f90 \
$(main)/magnetic_spectral_function.f90

OBJ3 = $(main)/stdtypes.f90   \
$(main)/mtprng.f90   \
$(main)/hamiltonian.f90  \
$(main)/getprobdensity.f90

OBJ4 = $(main)/stdtypes.f90   \
$(main)/mtprng.f90   \
$(main)/hamiltonian.f90  \
$(main)/sort.f90  \
$(main)/DOSautocorrelation.f90

#-------- Default target and other targets


HF: 
	$(FC) -o $@ $(OBJ1) $(LDFLAGS)

Akw:
	$(FC) -o $@ $(OBJ2) $(LDFLAGS)

probden:
	$(FC) -o $@ $(OBJ3) $(LDFLAGS)	

DOScorr:
	$(FC) -o $@ $(OBJ4) $(LDFLAGS)

clean:
	rm -f *.mod HF.o*

