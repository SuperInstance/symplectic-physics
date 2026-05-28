# Compiler and flags
FC = gfortran
FFLAGS = -O2 -Wall -std=f2008 -fcheck=all
LDFLAGS = -llapack -lblas

# Directories
SRCDIR = src
TESTDIR = tests
BUILDDIR = build

# Source files
SRCS = $(SRCDIR)/symplectic.f90 \
       $(SRCDIR)/hamiltonian.f90 \
       $(SRCDIR)/integrators.f90 \
       $(SRCDIR)/conservation.f90 \
       $(SRCDIR)/nbody.f90

# Object files
OBJS = $(BUILDDIR)/symplectic.o \
       $(BUILDDIR)/hamiltonian.o \
       $(BUILDDIR)/integrators.o \
       $(BUILDDIR)/conservation.o \
       $(BUILDDIR)/nbody.o

# Module files get built in BUILDDIR
MODFLAGS = -J$(BUILDDIR)

.PHONY: all test clean lib

all: lib test

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Build module objects (order matters for dependencies)
$(BUILDDIR)/symplectic.o: $(SRCDIR)/symplectic.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

$(BUILDDIR)/hamiltonian.o: $(SRCDIR)/hamiltonian.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

$(BUILDDIR)/integrators.o: $(SRCDIR)/integrators.f90 $(BUILDDIR)/hamiltonian.o | $(BUILDDIR)
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

$(BUILDDIR)/conservation.o: $(SRCDIR)/conservation.f90 $(BUILDDIR)/hamiltonian.o | $(BUILDDIR)
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

$(BUILDDIR)/nbody.o: $(SRCDIR)/nbody.f90 $(BUILDDIR)/hamiltonian.o | $(BUILDDIR)
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

# Static library
lib: $(OBJS)
	ar rcs $(BUILDDIR)/libsymplectic.a $(OBJS)

# Test executable
test: $(BUILDDIR)/test_symplectic
	./$(BUILDDIR)/test_symplectic

$(BUILDDIR)/test_symplectic: $(TESTDIR)/test_symplectic.f90 $(OBJS)
	$(FC) $(FFLAGS) $(MODFLAGS) -I$(BUILDDIR) $< $(OBJS) -o $@ $(LDFLAGS)

clean:
	rm -rf $(BUILDDIR)
