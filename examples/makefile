OBJS =  makefile ../src/Scintillans/BlockMatrix.hs ../src/Scintillans/Solver.hs ../src/Scintillans/Synchrotron.hs

all: BlockMatrixTest SolverTest IntAiTest WTest SynchrotronTestClassical PhotonEmissionTest SynchrotronClassicalTest

BlockMatrixTest: $(OBJS) BlockMatrixTest.hs
	stack build && stack ghc -- -Odph -rtsopts -eventlog -fno-liberate-case -funfolding-use-threshold1000 -funfolding-keeness-factor1000 -fllvm -optlo-O3 --make BlockMatrixTest.hs && ./BlockMatrixTest +RTS -l

SolverTest: $(OBJS) SolverTest.hs
	stack build && stack ghc -- --make SolverTest.hs && ./SolverTest

IntAiTest: $(OBJS) IntAiTest.hs
	stack build && stack ghc -- --make IntAiTest.hs && ./IntAiTest

WTest: $(OBJS) WTest.hs
	stack build && stack ghc -- --make WTest.hs && ./WTest

SynchrotronTestClassical: $(OBJ) SynchrotronTestClassical.hs
	stack build && stack ghc -- -Odph -rtsopts -eventlog -fno-liberate-case -funfolding-use-threshold1000 -funfolding-keeness-factor1000 -fllvm -optlo-O3 --make SynchrotronTestClassical.hs && ./SynchrotronTestClassical +RTS -l

PhotonEmissionTest: $(OBJ) PhotonEmissionTest.hs
	stack build && stack ghc -- -O2 --make PhotonEmissionTest.hs && ./PhotonEmissionTest

SynchrotronClassicalTest: $(OBJ) SynchrotronClassicalTest.hs
	stack build && stack ghc -- -O2 --make SynchrotronClassicalTest.hs && ./SynchrotronClassicalTest

clean:
	rm -I BlockMatrixTest SolverTest IntAiTest WTest SynchrotronTestClassical PhotonEmissionTest SynchrotronClassicalTest
