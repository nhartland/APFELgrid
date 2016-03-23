// The MIT License (MIT)

// Copyright (c) 2016 Nathan Hartland, Stefano Carrazza, Luigi Del Debbio

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

#include <stdio.h>
#include <iostream>
#include <cstdlib>

#include "APFELgrid/fkgenerator.h"
#include "APFELgrid/exceptions.h"

namespace NNPDF
{

  FKGenerator::FKGenerator( std::istream& is ):
  FKTable(is) 
  {

  };

  FKGenerator::~FKGenerator()
  {

  }

  // Fill the FKTable - beware there is no error checking performed!
  void FKGenerator::Fill(   int const& d,     // Datapoint index
                            int const& ix1,   // First x-index
                            int const& ix2,   // Second x-index
                            size_t const& ifl1,  // First flavour index
                            size_t const& ifl2,  // Second flavour index
                            double const& fk  // FK Value
                          )
  {

    if (fk==0)
      return;

    if (d >= fNData)
      throw RangeError("FKGenerator::Fill","datapoint " + ToString(d) + "out of bounds.");

    if (ix1 >= fNx)
      throw RangeError("FKGenerator::Fill","xpoint " + ToString(ix1) + " out of bounds.");

    if (ix2 >= fNx)
      throw RangeError("FKGenerator::Fill","xpoint " + ToString(ix2) + " out of bounds.");

    if (ifl1 >= 14)
      throw RangeError("FKGenerator::Fill","flavour " + ToString(ifl1) + " out of bounds.");

    if (ifl2 >= 14)
      throw RangeError("FKGenerator::Fill","flavour " + ToString(ifl2) + " out of bounds.");

    // pointer to FKTable segment
    const size_t jFL = 14*ifl1 + ifl2;
    const size_t iSig = d*fDSz+jFL*fTx+ix1*fNx+ix2 ;

    // Assign FK Table
    fSigma[iSig] += fk;

    return;
  };

  // DIS version of Fill
  void FKGenerator::Fill( int const& d,     // Datapoint index
                          int const& ix,    // x-index
                          size_t const& ifl,   // flavour index
                          double const& fk  // FK Value
                        )
  {

    if (fk==0)
      return;

    if (d >= fNData)
      throw RangeError("FKGenerator::Fill","datapoint " + ToString(d) + " out of bounds.");

    if (ix >= fNx)
      throw RangeError("FKGenerator::Fill","xpoint " + ToString(ix) + " out of bounds.");

    if (ifl >= 14)
      throw RangeError("FKGenerator::Fill","flavour " + ToString(ifl) + " out of bounds.");

    // pointer to FKTable segment
    const size_t iSig = d*fDSz+ifl*fNx+ix;

    // Assign FK Table
    fSigma[iSig] += fk;

    return;
  };

}
