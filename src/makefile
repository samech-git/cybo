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
flags = 
exec = cybo

source = modules.o IO.o cybo.o

$(exec) : $(source)
	$(f90) $(source) $(flags) -o $(exec)

%.o: %.f90
	$(f90) $(flags) -c $<

clean:
	rm -f *.o *.mod
clobber: 
	rm -f *.o *.mod $(exec)