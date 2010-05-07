#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!        /\        CYBO-Unstuctured Euler Solver  !!!
#!!!       /  \                                      !!!
#!!!      /    \      Chris M. Yu & Britton J. Olson !!!
#!!!     /______\     Department on Aero/Astro       !!!
#!!!    /\      /\    Stanford University            !!!
#!!!   /  \    /  \   Spring 2010                    !!!
#!!!  /    \  /    \                                 !!! 
#!!! /______\/______\                                !!!
#!!!                                                 !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

f90 = ifort
#f90 = gfortran
flags = -O3
exec = cybo

source = modules.o IO.o solver.o cybo.o

$(exec) : $(source)
	$(f90) $(source) $(flags) -o $(exec)
	mv $(exec) ../bin/.

%.o: %.f90
	$(f90) $(flags) -c $<

clean:
	rm -f *.o *.mod ../bin/$(exec)
