//transform.h

namespace NNPDF_APFELgrid
{
    // LHA-style flavour basis
    enum {TBAR,BBAR,CBAR,SBAR,UBAR,DBAR,GLUON,D,U,S,C,B,T,PHT};
    // NNPDF_APFELgrid-style EVLN basis
    enum evlnBasis {  EVLN_GAM, EVLN_SNG, EVLN_GLU, EVLN_VAL, EVLN_V3,
                      EVLN_V8, EVLN_V15, EVLN_V24, EVLN_V35,  EVLN_T3,
                      EVLN_T8, EVLN_T15, EVLN_T24, EVLN_T35 };

  /**
   * Rotate flavour basis PDFs to evolution basis
   * \param LHA the les houches 13 pdfs
   * \return evln the lha in the evln basis
   */
  template<class inType, class outType>
  void LHA2EVLN(const inType *LHA, outType *EVLN)
  {
    const outType uplus = LHA[U] + LHA[UBAR];
    const outType uminus = LHA[U] - LHA[UBAR];

    const outType dplus = LHA[D] + LHA[DBAR];
    const outType dminus = LHA[D] - LHA[DBAR];

    const outType cplus = LHA[C] + LHA[CBAR];
    const outType cminus = LHA[C] - LHA[CBAR];

    const outType splus = LHA[S] + LHA[SBAR];
    const outType sminus = LHA[S] - LHA[SBAR];

    const outType tplus = LHA[T] + LHA[TBAR];
    const outType tminus = LHA[T] - LHA[TBAR];

    const outType bplus = LHA[B] + LHA[BBAR];
    const outType bminus = LHA[B] - LHA[BBAR];

    EVLN[0]= LHA[PHT]; // photon
    EVLN[1]=(uplus + dplus + cplus + splus + tplus + bplus); //Singlet
    EVLN[2]=(LHA[GLUON]); // Gluon

    EVLN[3]=( uminus + dminus + sminus + cminus + bminus + tminus ); //V
    EVLN[4]=( uminus - dminus ); // V3
    EVLN[5]=( uminus + dminus - 2*sminus); // V8
    EVLN[6]=( uminus + dminus + sminus - 3*cminus); //V15
    EVLN[7]=( uminus + dminus + sminus + cminus - 4*bminus ); //V24
    EVLN[8]=( uminus + dminus + sminus + cminus + bminus - 5*tminus); // V35

    EVLN[9]=(  uplus - dplus ); // T3
    EVLN[10]=( uplus + dplus - 2*splus ); // T8
    EVLN[11]=( uplus + dplus + splus - 3*cplus ); //T15
    EVLN[12]=( uplus + dplus + splus + cplus - 4*bplus ); //T24
    EVLN[13]=( uplus + dplus + splus + cplus + bplus - 5*tplus ); // T35
  }

}
