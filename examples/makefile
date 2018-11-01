OBJS = ../src/Scintillans/BlockMatrix.hs ../src/Scintillans/Solver.hs makefile

all: BlockMatrixTest SolverTest

BlockMatrixTest: $(OBJS) BlockMatrixTest.hs
	stack build && stack ghc -- -Odph -rtsopts -eventlog -fno-liberate-case -funfolding-use-threshold1000 -funfolding-keeness-factor1000 -fllvm -optlo-O3 --make BlockMatrixTest.hs && ./BlockMatrixTest +RTS -l

SolverTest: $(OBJS) SolverTest.hs
	stack build && stack ghc -- --make SolverTest.hs && ./SolverTest
