FC   = gfortran
Flag = -fopenmp -ffree-form -ffree-line-length-none -O3 -Wall -std=f2003 -fcheck=bounds
files = mod_fftw.f03 mod_param.f03 mod_operator.f03 mod_io.f03 mod_solver.f03 main.f03

FFTW3_DIR=/opt/fftw3
# FFTW3_DIR=/usr/physics/fftw/fftw-3.3.3

# fftwC = -L${FFTW3_DIR}/lib -lfftw3 -lfftw3_omp -lm 
fftwC = -L${FFTW3_DIR}/lib -lfftw3 -lm 
fftwI = -I${FFTW3_DIR}/include

cfit: 	$(files) mod_fftw.mod mod_io.mod mod_param.mod mod_operator.mod mod_solver.mod
	$(FC) ${Flag} ${fftwC} ${fftwI} -o $@ $(files)

clean:
	rm -f cfit *.mod *.o

mod_io.mod: mod_param.f03 mod_operator.f03 mod_io.f03 
	$(FC) -c $^
mod_operator.mod: mod_param.f03 mod_operator.f03
	$(FC) -c $^
mod_solver.mod: mod_fftw.f03 mod_param.f03 mod_io.f03 mod_operator.f03 mod_solver.f03 
	$(FC) ${fftwI} ${fftwC} -c $^
mod_fftw.mod: mod_param.f03 mod_operator.f03 mod_io.f03 mod_fftw.f03 
	$(FC) ${fftwI} ${fftwC} -c $^ 
mod_param.mod: mod_param.f03
	$(FC) -c $^

%.o: %.f03
	$(FC) -c $<