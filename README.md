# oned_rd
code to test 1D reaction diffusion

to compile, type: make

#EXECUTABLES
rd_pair_oned_rewgt
--calculates the association probability for a pair of particles in 1D vs time: p_assoc(t|x0), given an initial separation x0. Takes in x0 as a commandline input: e.g. if x0=1.1-->
./rd_pair_oned_rewght 1.1 

--Model parameters are hard-coded (ka, D, sigma)

rd_trap_rewgt_1D
--calculate the probability of being in a bound state between a trap molecule A surrounded by NB B molecules in 1D. theory is pBound= Keq*B0/(1+Keq*B0), where B0=NB/(L-sigma) (units of nm), and Keq=ka/kb (units of nm). outputs pBound vs time. Also outputs means of pBound for each repeat, performs multiple repeats. Takes NB as a commandline input: e.g. if NB=15
./rd_trap_rewgt_1D 15

rd_trap_Events_rewgt1D
--calculates bound probability, but stops when maxevents have converged, so at the simulation end, the pBound value should be closer to correct. Not as good for error estimates, and do not time average because outputs only when events occur.


