// example_conv.cc
// This example demonstrates how an FK table may be convoluted with a given PDF

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <cstdio>

#include "LHAPDF/LHAPDF.h"

#include "appl_grid/appl_grid.h"
#include "NNPDF/fastkernel.h"
#include "NNPDF/transform.h"

extern "C" void evolvepdf_(const double& , const double& , double* );
static double* lha_pdf;
void fkpdf (const double& x, const double& Q, const size_t& n, NNPDF::real* pdf)
{
  evolvepdf_(x,Q,lha_pdf);
  NNPDF::LHA2EVLN<double>(lha_pdf, pdf);
}

int main(int argc, char* argv[]) {
	// Initialise temp. array
	lha_pdf = new double[13];

	// Read FK table
	std::ifstream infile; infile.open("test.dat");
	NNPDF::FKTable FK(infile);

	// Initialise PDF set
	LHAPDF::initPDFSet("CT10", LHAPDF::LHGRID, 0);
	NNPDF::real* results = new NNPDF::real[FK.GetNData()];
	FK.Convolute(fkpdf, 1, results);

	for (size_t i=0; i < FK.GetNData(); i++)
		std::cout << results[i] <<std::endl;

	delete[] results;
	delete[] lha_pdf;
	exit(0);
}



