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

#include <string>
#include "fastkernel.h"

// Forward decls
namespace appl { class grid; }

namespace APFELgrid
{
  // ************************ Public utility functions **************************
  // These functions are for use when initialising APFEL kinematics such that
  // they are suitable for the relevant APPLgrid.

  // Returns the minimum value in the APPLgrid x-grid
  // If nonzero is true, returns the smallest x-value regardless of weight
  // if false, returns the smallest x-value associated with a non-zero weight
  double get_appl_Xmin(appl::grid const& g, const bool& nonzero);

  // Returns the minimum (Q2min) and maximum (Q2max) reach in scale of the applgrid g
  void get_appl_Q2lims(appl::grid const& g, double& Q2min, double& Q2max);

  // ************************ FK table computation **************************
  // Performs the combination of an APPLgrid g with evolution factors provided
  // by APFEL, resulting in a new FK table. Required arguments are the initial
  // scale for the FK tables (Q0), the name of the produced table (name), the appl::grid (g),
  // the path to the appl::grid (gridfile) and an optional appl::grid directory (directory).
  NNPDF_APFELgrid::FKTable<double>* computeFK( double const& Q0, std::string const& name, appl::grid const& g, std::string const& gridfile, std::string directory = "grid" );
}

namespace NNPDF_APFELgrid
{
    /**
    * \class FKGenerator
    * \brief Class for filling FastKernel tables
    */
    class FKGenerator : public FKTable<double>
    {
      private:
        FKGenerator();                          //!< Disable default constructor
        FKGenerator( const FKGenerator& );      //!< Disable default copy-constructor
        FKGenerator& operator=(const FKGenerator&); //!< Disable copy-assignment

      public:
        FKGenerator( std::istream& ); // Constructor
        virtual ~FKGenerator(); //!< Destructor

         // Fill the FKTable
        void Fill(  int const& d,     // Datapoint index
                    int const& ix1,   // First x-index
                    int const& ix2,   // Second x-index
                    size_t const& ifl1,  // First flavour index
                    size_t const& ifl2,  // Second flavour index
                    double const& isig  // FK Value
                  );

        // DIS version of Fill
        void Fill(  int const& d,     // Datapoint index
                    int const& ix,    // x-index
                    size_t const& ifl,   // flavour index
                    double const& isig  // FK Value
                  );
    };
}
