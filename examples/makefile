OBJS = ../src/Scintillans/BlockMatrix.hs ../src/Scintillans/Solver.hs makefile

all: BlockMatrixTest SolverTest IntAiTest WTest

BlockMatrixTest: $(OBJS) BlockMatrixTest.hs
	stack build && stack ghc -- -Odph -rtsopts -eventlog -fno-liberate-case -funfolding-use-threshold1000 -funfolding-keeness-factor1000 -fllvm -optlo-O3 --make BlockMatrixTest.hs && ./BlockMatrixTest +RTS -l

SolverTest: $(OBJS) SolverTest.hs
	stack build && stack ghc -- --make SolverTest.hs && ./SolverTest

IntAiTest: $(OBJS) IntAiTest.hs
	stack build && stack ghc -- --make IntAiTest.hs && ./IntAiTest

WTest: $(OBJS) WTest.hs
	stack build && stack ghc -- --make WTest.hs && ./WTest

ProbabilityTestClassical: $(OBJ) ProbabilityTestClassical.hs
	stack build && stack ghc -- --make ProbabilityTestClassical.hs && ./ProbabilityTestClassical
