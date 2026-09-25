# Linux / WSL:          make
# Windows (MinGW-w64):  mingw32-make CC=gcc EXE=.exe
# Without OpenMP:       make OMP=
CC ?= gcc
OMP ?= -fopenmp
EXE ?=
CFLAGS ?= -O3 -march=native -std=c99 -D_POSIX_C_SOURCE=200809L -Wall -Wextra $(OMP)
LDLIBS = -lm

fdtd3d_oblique$(EXE): fdtd3d_oblique.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f fdtd3d_oblique fdtd3d_oblique.exe

.PHONY: clean
