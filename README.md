# Grid 

This branch holds the modifications to Grid made in preparation for the paper
"Symplectic lattice gauge theories in the Grid framework: domain wall fermions and continuum extrapolations",
to appear.

Notable changes:

 - Code used for HMC in `HMC/MobiusFundnf2.cc`
 - Fermion instantiations added in `Grid/qcd/action/fermion/instantiations`,
   and corresponding typedefs to enable Hadrons to use them.
 - Code used for propagator test in `tests/Test_propagator.cc`
 - Code used for preliminary (unsuccessful) eigenvalue tests in `tests/lanczos`
