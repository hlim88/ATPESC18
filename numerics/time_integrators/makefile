# Define the MFEM directories here
MFEM_DIR=/projects/ATPESC2018/FASTMath/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/mfem-develop-shrxlv7z4whbh2476qpdhuwtf6frzatv
PETSC_DIR=/projects/ATPESC2018/FASTMath/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/petsc-develop-xufshydc44bard4v7qinrl2nfetl7fcm

CONFIG_MK = $(MFEM_DIR)/share/mfem/config.mk
-include $(CONFIG_MK)
.SUFFIXES: .o .cpp .mk
.PHONY: all clean
.cpp.o:
	$(MFEM_CXX) $(MFEM_FLAGS) -c $<

all: elasticity-snes advection-ode

elasticity-snes: elasticity-snes.o
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS) $(LDFLAGS)

advection-ode: advection-ode.o
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS) $(LDFLAGS)

clean:
	rm -f *.o *~ elasticity-snes advection-ode adv-init.0* adv-final.0* adv-mesh.0* elastic_energy.0* velocity.0* deformed.0*
	rm -rf *.dSYM *.TVD.*breakpoints elas_* adv_*
