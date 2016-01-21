// test0.cc
// Test APFEL initialisation

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <cstdio>

#include "APFEL/APFEL.h"

static const double Q0 = 1.0; // 1 GeV initial scale

int main(int argc, char* argv[]) {

  // Theory, perturbative order of evolution
  APFEL::SetTheory(string("QCD"));
  APFEL::SetAlphaQCDRef(0.118, 91.2);

  APFEL::SetPDFEvolution("exactalpha");
  APFEL::SetAlphaEvolution("exact");

  // Needs to match with the APPLgrid
  APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetMaxFlavourPDFs(5);

  APFEL::SetQLimits( Q0, 1000 );

  // Initialise
  APFEL::LockGrids(true);
  APFEL::EnableEvolutionOperator(true);
  APFEL::InitializeAPFEL();

    // Initialise FK table x-grid
  std::cout << APFEL::nIntervals() <<std::endl;

   exit(0);
}

