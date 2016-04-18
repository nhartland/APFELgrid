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
#include "APFELgrid/fastkernel.h"
#include "APFELgrid/transform.h"

typedef float ctype;

extern "C" void evolvepdf_(const double& , const double& , double* );
static double* lha_pdf;
void fkpdf (const double& x, const double& Q, const size_t& n, ctype* pdf)
{
  evolvepdf_(x,Q,lha_pdf);
  NNPDF::LHA2EVLN<double, ctype>(lha_pdf, pdf);
}

int main(int argc, char* argv[]) {
	// Initialise temp. array
	lha_pdf = new double[13];

	// Read FK table
	std::ifstream infile; infile.open("./tests/atlas-Z0-rapidity.fk");
	NNPDF::FKTable<ctype> FK(infile);

	// Initialise PDF set
	LHAPDF::initPDFSet("NNPDF30_nlo_as_0118", LHAPDF::LHGRID, 0);
	ctype* results = new ctype[FK.GetNData()];

	const int nIte = 1E5;
	timeval startTime; gettimeofday(&startTime, NULL);
	for (int i=0; i<nIte; i++) FK.Convolute(fkpdf, 1, results);
	timeval endTime; gettimeofday(&endTime, NULL);

	for (int i=0; i < FK.GetNData(); i++)
		std::cout << results[i] <<std::endl;

	double seconds  = endTime.tv_sec  - startTime.tv_sec;
    double useconds = endTime.tv_usec - startTime.tv_usec;
    std::cout << "Timer: " << (1E6f*seconds + useconds)/(nIte*FK.GetNData()) <<" us per point"<<std::endl;

	delete[] results;
	delete[] lha_pdf;
	exit(0);
}



