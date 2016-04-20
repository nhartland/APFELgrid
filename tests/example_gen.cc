// APFELgrid
// =========
// FastKernel table generation
// ---------------------------
// Here we shall demonstrate with a working example how an **APPLgrid** may be
// combined with evolution factors provided by **APFEL** to produce a **FastKernel**
// (**FK**) table

// To do this we require **APFEL**
#include "APFEL/APFEL.h"
// The **APFELgrid** plugin
#include "APFELgrid/APFELgrid.h"
#include "APFELgrid/fastkernel.h"
// and **APPLgrid**
#include "appl_grid/appl_grid.h"

// We will also demonstrate how said **FK** table may be written to and read from file.
// To do so we will need the following headers.
#include <iostream>
#include <fstream>

// Beginning the main loop
int main(int argc, char* argv[]) {

  // We start by initialising the general evolution parameters of **APFEL**.
  // Here we are using the exact solution of NLO DGLAP with a strong coupling
  // of 0.118 evaluated at the Z mass. Our evolution and alpha_S running is limited
  // to five flavours
  APFEL::SetTheory("QCD");
  APFEL::SetPDFEvolution("exactalpha");
  APFEL::SetAlphaEvolution("exact");

  APFEL::SetAlphaQCDRef(0.118, 91.2);
  APFEL::SetPerturbativeOrder(1);
  APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetMaxFlavourPDFs(5);

  // Next we read the **APPLgrid** through the usual procedure
  const std::string gridfile = "./tests/atlas-Z0-rapidity.root";
  const appl::grid  g(gridfile);

  // Once the **APPLgrid** is read, we can specify the initial-scale grid of x-points
  // to be used for the interpolation of the evolution factors. To do this
  // we make use of the **APFELgrid** helper function *get_appl_Xmin* which returns the lowest
  // x-value required by the **APPLgrid** in question.
  APFEL::SetNumberOfGrids(2);
  APFEL::SetGridParameters(1,15,5, APFELgrid::get_appl_Xmin(g, true) );
  APFEL::SetGridParameters(2,15,5,1e-1);

  // With the x-grid specified, we can now generate the **FK** table. Firstly we specify
  // an initial scale from which the evolution is to be performed. Here we are using an initial
  // scale of *Q0 = 1 GeV*.
  const double Q0 = 1.0;

  // The **FK** table is then generated with a single call to *computeFK*. Its arguments are
  // + *Q0*           - the requested initial scale
  // + *"ATLASZRAP"*  - the requested name of the **FK** table
  // + *g*            - the **APPLgrid** to be combined
  // + *gridfile*     - the path to the **APPLgrid** to be combined

  // *computeFK* will then return an *NNPDF::FKTable\<double\>*. The template parameter specifies
  // the precision of the internal representation of the FK table, along with the precision of its
  // output after convolution. For generation purposes this is fixed at double precision, however
  // the grid may be later read at single precision in order to enjoy the benefits of faster convolution.
  NNPDF::FKTable<double>* FK = APFELgrid::computeFK(Q0, "ATLASZRAP", g, gridfile);

  // Now the user may add various tags to the FK table as they wish, to help identify e.g
  // the theory parameters used in its generation later. For example:
  FK->AddTag(NNPDF::FKHeader::THEORYINFO, "TAG_KEY", "TAG_VALUE");
  FK->AddTag(NNPDF::FKHeader::THEORYINFO, "PertubativeOrder_asString", "NLO");
  FK->AddTag(NNPDF::FKHeader::THEORYINFO, "PertubativeOrder_asInt", 1);

  // Finally we write the new FK table to file, exit and close the main loop
  std::ofstream outfile; outfile.open("./tests/atlas-Z0-rapidity.fk");
  FK->Print(outfile); outfile.close();
  
  exit(0);
}



