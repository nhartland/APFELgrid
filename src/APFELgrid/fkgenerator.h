// The MIT License (MIT)

// Copyright (c) 2016 Stefano Carrazza, Luigi Del Debbio, Nathan Hartland

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

#pragma once

#include <string>
#include <vector>
#include "fastkernel.h"

namespace NNPDF
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