NAMESPACE_BEGIN(Grid);

template <class Gimpl> class WuppertalSmearing : public Gimpl
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField GaugeLorentz;

  template<typename T>
  static void WuppertalSmear(const std::vector<LatticeColourMatrix>& U,
                            T& chi,
                            const Real& step, int Iterations, int orthog)
  {
    GridBase *grid = chi.Grid();
    T psi(grid);

    Real coeff = step;

    int dims = Nd;
    if( orthog < Nd ) dims=Nd-1;
    double norm = 1/(1 + 2.0*dims*coeff);
    std::cout<< GridLogMessage << "Normalization factor : "<< norm <<std::endl;

    for(int n = 0; n < Iterations; ++n) {
      psi = chi;
      for(int mu=0;mu<Nd;mu++) {
        if ( mu != orthog ) {
          psi = psi + coeff*(Gimpl::CovShiftForward(U[mu],mu,chi));
          psi = psi + coeff*(Gimpl::CovShiftBackward(U[mu],mu,chi));
        }
      }
      chi = norm*psi;
    }
  }
};

NAMESPACE_END(Grid);

