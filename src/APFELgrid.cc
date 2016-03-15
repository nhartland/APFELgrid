// The MIT License (MIT)

// Copyright (c) 2016 Nathan Hartland

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "APFELgrid/APFELgrid.h"
#include "APFELgrid/fkgenerator.h"

#include "APFEL/APFEL.h"

#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_igrid.h"

#include <math.h>

namespace APFELgrid{

  // ************************ Public utility functions **************************
  // These are available through the APFELgrid header, for use when initialising
  // APFEL kinematics such that they are suitable for the relevant APPLgrid.


  // Returns the minimum value in the APPLgrid x-grid
  // If nonzero is true, returns the smallest x-value regardless of weight
  // if false, returns the smallest x-value associated with a non-zero weight
  double get_appl_Xmin(appl::grid const& g, const bool& nonzero)
  {
    double xmin = 1.0;
    for(int i=0; i<g.nloops()+1; i++)
      for (int j=0; j<g.Nobs(); j++)
      {
        appl::igrid *igrid = const_cast<appl::igrid*>(g.weightgrid(i, j));
        for (int ix1=0; ix1<igrid->Ny1(); ix1++)
          for (int ix2=0; ix2<igrid->Ny2(); ix2++)
            for (int t=0; t<igrid->Ntau(); t++) 
              for (int ip=0; ip<g.subProcesses(i); ip++)
                {
                  // Associated weight
                  const double point_weight = (*(const SparseMatrix3d*) igrid->weightgrid(ip))(t,ix1,ix2);
                  if ( point_weight != 0 || !nonzero) 
                  {
                    xmin = std::min(xmin, igrid->fx(igrid->gety1(ix1)));
                    xmin = std::min(xmin, igrid->fx(igrid->gety2(ix2)));
                  }
                }
      }
    
    return xmin;
  }

  // Returns the minimum (Q2min) and maximum (Q2max) reach in scale of the applgrid g
  void get_appl_Q2lims(appl::grid const& g, double& Q2min, double& Q2max)
  {
    Q2max = g.weightgrid(0, 0)->getQ2max();
    Q2min = g.weightgrid(0, 0)->getQ2min();

    for(int i=0; i<g.nloops()+1; i++)
      for (int j=0; j<g.Nobs(); j++)
      {
        appl::igrid const *igrid = g.weightgrid(i, j);
        Q2max = std::max(Q2max, igrid->getQ2max());
        Q2min = std::max(Q2min, igrid->getQ2min());
      }
  }

  // ************************ Internal utility functions **************************
  // These functions are used as helpers for the main FK table combination routine.

  // Splits an APPLgrid PDF string into it's component parts in order that the
  // appropriate genPDF may be loaded.
  std::vector<std::string> splitpdf ( std::string const& str )
  {
    std::vector<std::string> outvec;
    std::stringstream ss(str); std::string s;
    while (getline(ss, s, ':')) 
      outvec.push_back(s);
    return outvec;
  }

  // Returns the perturbative order to be used in the combination.
  // This is primarily set in APFEL. If the requested perturbative
  // order in APFEL is higher than the maximum available in the APPLgrid,
  // a warning is shown, and the combination is performed at max(APPLgrid),
  // with the full accuracy specified in APFEL. To acheive the full requested
  // accuracy K-factors must be applied.
  double get_ptord( appl::grid const& g)
  {
    const int apfel_pto = APFEL::GetPerturbativeOrder();
    const int appl_pto  = g.nloops();
    if (appl_pto < apfel_pto)
    {
      std::cout << "Warning: requested perturbative order not available in APPLgrid." <<std::endl;
      std::cout << "Using perturbative order " << apfel_pto << " in APFEL and " << appl_pto <<" in APPLgrid"<<std::endl;
    }
    return std::min(appl_pto, apfel_pto) + 1; // +1 to use as maximum in APPLgrid loop
  }

  // Translates 'loop' order to appl::grid index
  // This is specifically in order to translate aMC@NLO four-part grids
  // into LO and NLO components.
  // aMC@NLO convolution uses Born = grid 3
  //                          NLO  = grid 0
  int get_grid_idx( appl::grid const& g, int const& pto )
  {
    if (g.calculation() == appl::grid::AMCATNLO)
      return (pto==0) ? 3:0;
    return pto;
  }

  // Returns the APPLgrid PDF object associated with the ith subgrid of g
  appl::appl_pdf* get_appl_pdf( appl::grid const& g, int const& i )
  {
    const std::string pdfnames = g.getGenpdf();
    std::vector<std::string> pdfvec = splitpdf( pdfnames );
    const size_t isubproc = pdfvec.size() == 1 ? 0:i;
    return appl::appl_pdf::getpdf( pdfvec[isubproc] );
  }

  // Returns the minimum and maxiumum x-grid points for a specified subgrid slice.
  // igrid is the requested subgrid, nsubproc the number of subprocesses held within igrid.
  // tau specified the bin in scale to be investigated, and alpha specifies the bin in x1.
  // nxlow and nxhigh return the minumum and maxiumum bins in x2 respectively.
  void get_igrid_limits(appl::igrid const* igrid, int const& nsubproc, int const& tau, int const& alpha, int& nxlow, int& nxhigh)
  {
    // Need to remove const due to APPLgrid non-const methods
    appl::igrid* igrid_nc = const_cast<appl::igrid*>(igrid);
    nxlow=igrid->Ny2(); nxhigh=0;
    tsparse1d<double> *ts1;
    for (int tsp=0; tsp<nsubproc; tsp++)
      if ((*(const SparseMatrix3d*) igrid_nc->weightgrid(tsp))[tau] != NULL)
        if (( ts1 = (*(const SparseMatrix3d*) igrid_nc->weightgrid(tsp))[tau]->trimmed(alpha) ) != NULL)
        {                
          nxlow  = std::min(ts1->lo(), nxlow);
          nxhigh = std::max(ts1->hi(), nxhigh);
        }
  }

  // Computes the normalisation required for translating APPLgrid weights to FK ones
  // e.g including factors of alpha_S, bin width etc.
  // g is the APPLgrid being combined with evolution factors.
  // pto specified the perturbative order being combined, as the value of alpha_S in the current bin,
  // and x1/x2 specify the numerical values of the PDF x-values for the first and second PDF respectively.
  double compute_wgt_norm(appl::grid const& g, int const& d, double const& pto, double const& as, double const& x1, double const& x2)
  {
    // PDF x and bin width normalisation
    double norm = 1.0/(x1*x2*g.deltaobs(d));

    // Normalisation by number of runs
    appl::grid& g_nc = const_cast<appl::grid&>(g);
    if ( !g.getNormalised() && g_nc.run() )
      norm*=1.0/double(g_nc.run());

    // Factor of alpha_S
    const double LO = g.leadingOrder();
    if (g.calculation() == appl::grid::AMCATNLO)
      norm*=pow(as*(4.0*M_PI), LO+pto );
    else
      norm*=pow(as/(2.0*M_PI), LO+pto );

    return norm;
  }

  // Allocates the arrays used for evolution factors
  double*** alloc_evfactor()
  {
    const int nxin = APFEL::nIntervals();
    std::cout << "Alloc " << nxin <<std::endl;
    double*** f = new double**[nxin];
    for (int i=0; i<nxin; i++)
    {
      f[i] = new double*[14];  // These are in EVLN basis (photon!)
      for (size_t j=0; j<14; j++)
        f[i][j] = new double[13]; // These are in APPLGRID basis (no photon!)
    }
    return f;
  }

  // Frees evolution factor array f
  void free_evfactor( double*** f )
  {
    const int nxin = APFEL::nIntervals();
    for (int i=0; i<nxin; i++)
    {    
      for (size_t j=0; j<14; j++)
        delete[] f[i][j];
      delete[] f[i];
    }
    delete[] f;
  }

  // Computes the requested evolution factors, for initial scale Q0, target scale Q1 and output x-value xo.
  void compute_evfactors( double const& Q0, double const& Q1, double const& xo, double*** fA )
  {
    // APFEL parameters
    const int nxin = APFEL::nIntervals();
    const double Q0diff = fabs(Q0-APFEL::GetMuF0());
    const double Q1diff = fabs(Q1-APFEL::GetMuF());

    // Recalculate if not cached
    if ( (Q0diff > 1E-10) || (Q1diff > 1E-10) )
        APFEL::EvolveAPFEL(Q0,Q1);

    for (int xi = 0; xi < nxin; xi++)
      for (size_t fi = 0; fi < 14; fi++)
        for(int i=0; i<13; i++)
          fA[xi][fi][i] = APFEL::ExternalEvolutionOperator(std::string("Ev2Ph"),i-6,fi,xo,xi);
  }


  // ************************ FK Table computation **************************
  // These functions provide the tools to initialise and generate FK tables.

  // Generates a new FKGenerator class, given a base appl::grid, initial scale and setname.
  // Physics and interpolation parameters are obtained directly from APFEL
  NNPDF::FKGenerator* generate_FK( appl::grid const& g, double const& Q0, std::string const& setname)
  {
    // Generate FKTable header
    NNPDF::FKHeader FKhead;
    FKhead.AddTag(NNPDF::FKHeader::BLOB, "GridDesc", g.getDocumentation());
    FKhead.AddTag(NNPDF::FKHeader::GRIDINFO, "SETNAME", setname);
    FKhead.AddTag(NNPDF::FKHeader::GRIDINFO, "NDATA", g.Nobs());
    FKhead.AddTag(NNPDF::FKHeader::GRIDINFO, "HADRONIC", true);
    FKhead.AddTag(NNPDF::FKHeader::VERSIONS, "APFEL", APFEL::GetVersion());
    FKhead.AddTag(NNPDF::FKHeader::THEORYINFO, "Q0", Q0 );

    // x-grid header
    const int nx = APFEL::nIntervals();
    FKhead.AddTag(NNPDF::FKHeader::GRIDINFO, "NX", nx );
    std::stringstream xGheader;
    for (int i=0; i<nx; i++)
      xGheader << std::setprecision(16) << std::scientific << APFEL::xGrid(i) <<std::endl;
    FKhead.AddTag(NNPDF::FKHeader::BLOB, "xGrid", xGheader.str());

    // Full flavourmap
    std::stringstream fMapHeader;
    for (int i=0; i<14; i++)
    {
      for (int i=0; i<14; i++)
        fMapHeader << "1 ";
      fMapHeader<<std::endl;
    }
    FKhead.AddTag(NNPDF::FKHeader::BLOB, "FlavourMap", fMapHeader.str());

    std::stringstream IO; FKhead.Print(IO);
    return new NNPDF::FKGenerator( IO );
  }

  // Performs the combination of an APPLgrid g with evolution factors provided
  // by APFEL, resulting in a new FK table. Required arguments are the name of the produced table 'name', 
  // the appl::grid g, the initial scale for the FK tables Q0, and a boolean specifying whether or not the
  // APPLgrid has PDF weights enabled (this is at the moment impossible to tell from the APPLgrid API) 
  NNPDF::FKTable* computeFK( std::string const& name, appl::grid const& g, double const& Q0, bool const& pdfwgt )
  {
    // Set APFEL scale limits
    double Qmin, Qmax;
    get_appl_Q2lims(g, Qmin, Qmax);
    Qmin = std::sqrt(Qmin); Qmax = std::sqrt(Qmax);
    APFEL::SetQLimits( std::min(Q0, Qmin), Qmax );
    
    // Initialise APFEL
    APFEL::LockGrids(true);
    APFEL::SetFastEvolution(false);
    APFEL::EnableEvolutionOperator(true);
    APFEL::InitializeAPFEL();
    APFEL::EvolveAPFEL(Q0, Q0);

    // Setup FK table and evolution factor arrays
    NNPDF::FKGenerator* FK = generate_FK(g, Q0, name);
    double*** fA = alloc_evfactor();
    double*** fB = alloc_evfactor();
  
    for (int d=0; d<g.Nobs(); d++)
    {
      std::cout << d<<"/"<<g.Nobs() << " points completed "<<std::endl;
      for (size_t pto=0; pto < get_ptord(g); pto++)
      {
        const int gidx = get_grid_idx(g, pto);          // APPLgrid grid index
        appl::appl_pdf *genpdf = get_appl_pdf(g, gidx); // APPLgrid pdf generator

        // define subprocess weight array W, and parton density array H
        const size_t nsubproc = g.subProcesses(gidx);
        double *W = new double[nsubproc];
        double *H = new double[nsubproc];
        
        // Fetch grid pointer
        appl::igrid const *igrid = g.weightgrid(gidx, d);
        for (int t=0; t<igrid->Ntau(); t++) // Loop over scale bins
        {
          const double Q   = sqrt( igrid->fQ2( igrid->gettau(t)) );
          const double as  = APFEL::AlphaQCD(Q);
          
          for (int a=0; a<igrid->Ny1(); a++  ) // Loop over x1 bins
          {
            // Get trimmed limits
            int nxlow, nxhigh;
            get_igrid_limits(igrid, nsubproc, t, a, nxlow, nxhigh);

            // Compute evolution factors for first PDF, only if there are nonzero weights
            const double x1 = igrid->fx(igrid->gety1(a));
            if (nxlow <= nxhigh)
              compute_evfactors(Q0, Q, x1, fA);
            
            for (int b=nxlow; b<=nxhigh; b++) // Loop over x2 bins
            {
              // fetch weight values
              bool nonzero=false;
              for (size_t ip=0; ip<nsubproc; ip++)
                if (( W[ip] = (*(const SparseMatrix3d*) const_cast<appl::igrid*>(igrid)->weightgrid(ip))(t,a,b) )!=0)
                  nonzero=true;
              
              // If nonzero, perform combination
              if (nonzero)
              {
                // Calculate normalisation factor
                const double x2 = igrid->fx(igrid->gety2(b));
                const double pdfnrm =  pdfwgt ? igrid->weightfun(x1)*igrid->weightfun(x2) : 1.0;
                const double norm = pdfnrm*compute_wgt_norm(g, d, pto, as, x1, x2);
                
                // Compute evolution factors for second PDF
                compute_evfactors(Q0, Q, x2, fB);

                const size_t nxin = APFEL::nIntervals();
                for (size_t i=0; i<nxin; i++) // Loop over input pdf x1
                for (size_t j=0; j<nxin; j++) // Loop over input pdf x2
                for (size_t k=0; k<14; k++) // loop over flavour 1
                for (size_t l=0; l<14; l++) // loop over flavour 2
                  {
                    // Rotate to subprocess basis and fill
                    genpdf->evaluate(fA[i][k],fB[j][l],H);
                    for (size_t ip=0; ip<nsubproc; ip++)
                      if (W[ip] != 0 and H[ip] != 0)
                        FK->Fill( d, i, j, k, l, norm*W[ip]*H[ip] );
                  }
              }
            }
          }
        }
        // Free subprocess arrays
        delete[] W;
        delete[] H;
      }
    }

    // Cleanup evolution factors
    free_evfactor(fA);
    free_evfactor(fB);
  
    return FK;
  }


}