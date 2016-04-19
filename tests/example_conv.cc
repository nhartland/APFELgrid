// APFELgrid
// =========
// FastKernel table convolution
// ----------------------------
// This example demonstrates how an **FK** table may be read from file and convoluted with a given PDF.

// For this demonstration, we need some standard headers
#include <iostream>
#include <cstdlib>

// Along with **LHAPDF** to provide initial scale PDFs
#include "LHAPDF/LHAPDF.h"

// And the required **APFELgrid** headers
#include "APFELgrid/fastkernel.h"
#include "APFELgrid/transform.h"

// We start with some boilerplate, a typedef named *ctype* so that we may switch between
// double and single precision convolutions simply, and a handle for the **LHAPDF** (v5-style)
// interface for **PDF** evolution (as per standard **APPLgrid** procedure).
typedef float ctype;
extern "C" void evolvepdf_(const double& , const double& , double* );

// Now we setup an equivalent bit of boilerplate for the **FK** convolution.
// Unlike **APPLgrids**, **FK** tables require PDFs in the DGLAP or Evolution (EVLN) basis
// (see the reference manual for details). **APFELgrid** provides a helper utility to
// perform the rotation to the EVLN basis from the **LHAPDF** basis. 
// No implicit rotation is performed, as this would require prior knowledge of the
// internal representation of PDFs used in the fit.

// The function *fkpdf* therefore interfaces the **LHAPDF** call *evolvepdf_* with
// the **FK** convolution by means of this rotation. The double *lha_pdf* holds the
// intermediate values in the **LHAPDF** basis. Here the member argument *n* is unused.
static double* lha_pdf;
void fkpdf (const double& x, const double& Q, const size_t& n, ctype* pdf)
{
  evolvepdf_(x,Q,lha_pdf);
  NNPDF::LHA2EVLN<double, ctype>(lha_pdf, pdf);
}


// With the boilerplate completed, we start the main loop by initialising the *lha_pdf* array
int main(int argc, char* argv[]) {
	lha_pdf = new double[13];

	// The **FK** table is then read from file, and a PDF set is initialised
	std::ifstream infile; infile.open("./tests/atlas-Z0-rapidity.fk");
	NNPDF::FKTable<ctype> FK(infile);
	LHAPDF::initPDFSet("NNPDF30_nlo_as_0118", LHAPDF::LHGRID, 0);

	// We now allocate an array (of type *ctype*) for results and perform the convolution.
	// The arguments to *FK::Convolute* here are
	// + *fkpdf* - a pointer to a function providing PDFs for individual values of *x*, *Q* and member *n* as defined above
	// + *1* - Here we only require one PDF member to be convoluted, but this may be increased if required.
	// + *results* - The results array
	ctype* results = new ctype[FK.GetNData()];
	FK.Convolute(fkpdf, 1, results);

	// For a single member convolution, the results of the product can be simply displayed as so
    for (int i=0; i < FK.GetNData(); i++)
		std::cout << results[i] <<std::endl;

	// Finally we clean up and end the program.
	delete[] results;
	delete[] lha_pdf;
	exit(0);
}



