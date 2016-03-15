// example_gen.cc
// This example demonstrates how an APPLgrid can be combined
// with evolution factors provided by APFEL to produce a FastKernel (FK)
// table.

#include <iostream>
#include <fstream>
#include <string>

#include "APFEL/APFEL.h"
#include "APFELgrid/APFELgrid.h"
#include "APFELgrid/fastkernel.h"

#include "appl_grid/appl_grid.h"

static const double Q0 = 1.0; // 1 GeV initial scale

int main(int argc, char* argv[]) {

  // Initialise APFEL
  APFEL::SetTheory(string("QCD"));
  APFEL::SetAlphaQCDRef(0.118, 91.2);
  APFEL::SetPerturbativeOrder(1);

  APFEL::SetPDFEvolution("exactalpha");
  APFEL::SetAlphaEvolution("exact");

  // Needs to match with the APPLgrid
  APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetMaxFlavourPDFs(5);

  // Read APPLgrid
  appl::grid g("./tests/atlas-Z0-rapidity.root");

  // Set x-grid
  APFEL::SetNumberOfGrids(2);
  APFEL::SetGridParameters(1,15,5, APFELgrid::get_appl_Xmin(g, true) );
  APFEL::SetGridParameters(2,15,5,1e-1);

  // Generate FK table
  NNPDF::FKTable* FK = APFELgrid::computeFK("ATLASZRAP", g, Q0, false);

  // Write FK table to file
  std::ofstream outfile; outfile.open("./tests/atlas-Z0-rapidity.fk");
  FK->Print(outfile); outfile.close();
  
  exit(0);
}



