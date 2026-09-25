# Linux / WSL:          make
# Windows (MinGW-w64):  mingw32-make CC=gcc EXE=.exe
# Without OpenMP:       make OMP=
# Frozen uniform reference (regression parity, SPEC_nonuniform §20.1 0-1a):  make ref
CC ?= gcc
OMP ?= -fopenmp
EXE ?=
CFLAGS ?= -O3 -march=native -std=c99 -D_POSIX_C_SOURCE=200809L -Wall -Wextra $(OMP)
LDLIBS = -lm

fdtd3d_oblique$(EXE): fdtd3d_oblique.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

ref: fdtd3d_oblique_ref$(EXE)

fdtd3d_oblique_ref$(EXE): ref/fdtd3d_oblique_uniform_ref.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f fdtd3d_oblique fdtd3d_oblique.exe fdtd3d_oblique_ref fdtd3d_oblique_ref.exe

.PHONY: clean ref
