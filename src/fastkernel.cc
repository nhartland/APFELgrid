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
#include <fstream>
#include <iostream>
#include <ios>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iterator>

#include "APFELgrid/fastkernel.h"
#include "APFELgrid/exceptions.h"

namespace NNPDF
{

  // Default verbosity level
  bool FKTable::Verbose = true;

  // ***************** Linear convolution *********************


#ifdef SSE_CONV
  #include <pmmintrin.h>

  static inline void convolute(const real* __restrict__ x, const real* __restrict__ y, real& retval, int const& n)
  {
    __m128 acc = _mm_setzero_ps();

    const __m128* a = (const __m128*) x;
    const __m128* b = (const __m128*) y;

    for (int i=0; i<n/4; i++)
      acc = _mm_add_ps(acc, _mm_mul_ps(*(a+i),*(b+i)));

    const __m128 shuffle1 = _mm_shuffle_ps(acc, acc, _MM_SHUFFLE(1,0,3,2));
    const __m128 tmp1     = _mm_add_ps(acc, shuffle1);
    const __m128 shuffle2 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(2,3,0,1));
    const __m128 tmp2     = _mm_add_ps(tmp1, shuffle2);

    _mm_store_ss(&retval,tmp2);
    _mm_empty();

    return;
  }

#else

  // Basic convolution
  static inline void convolute(const real* __restrict__ pdf,const real* __restrict__ sig,real& retval,int const& n)
  {
    for (int i = 0; i < n; i++)
      retval += pdf[i]*sig[i];
  }

#endif

  // ******************* Helpers ******************************

  /**
   * Split std::string into real std::vector
   */
  static void rsplit(std::vector<real>& results, std::string const& input)
  {
    std::stringstream strstr(input);
    std::istream_iterator<real> it(strstr);
    std::istream_iterator<real> end;
    
    results.assign(it, end);
    return;
  }

  // *********************************************************

  // Section delineators for FK headers
  static const int FK_DELIN_SEC = std::char_traits<char>::to_int_type('_');
  static const int FK_DELIN_BLB = std::char_traits<char>::to_int_type('{');
  static const int FK_DELIN_KEY = std::char_traits<char>::to_int_type('*');
  // *********************************************************





  /**
  * @brief Prototype Constructor for FKParser class - to be used only for constructing new FK tables
  * @param filename The FK table filename
  */
  FKHeader::FKHeader()
  {

  }

  /**
  * @brief Constructor for FK header parsing class
  * @param str The input stream
  */
  FKHeader::FKHeader(std::istream& str)
  {
     Read(str);
  }

    /**
  * @brief Constructor for FK header parsing class
  * @param filename The FK table filename
  */
  FKHeader::FKHeader(std::string const& filename)
  {
    std::ifstream instream;
    instream.open(filename);

    if (!instream.good())
        throw FileError("FKHeader::FKHeader","cannot open FK grid file: " + filename);

    Read(instream);
  }

  FKHeader::FKHeader(FKHeader const& ref):
  fVersions(ref.fVersions),
  fGridInfo(ref.fGridInfo),
  fTheoryInfo(ref.fTheoryInfo)
  {
 
  }

  FKHeader::~FKHeader()
  {

  }

  void FKHeader::Read(std::istream& is) 
  {
    if (!is.good())
        throw FileError("FKHeader::FKHeader","cannot open FK grid file ");

    int peekval = (is >> std::ws).peek();
    while  ( peekval == FK_DELIN_SEC ||
             peekval == FK_DELIN_BLB )
    {
      std::string sectionTitle; getline( is, sectionTitle );
      // Remove delineating character and trailing underscores
      sectionTitle = sectionTitle.substr(1, sectionTitle.size());
      sectionTitle.erase(std::remove(sectionTitle.begin(), sectionTitle.end(), '_'), sectionTitle.end());

      // Terminating case
      if (sectionTitle.compare("FastKernel") == 0)
        return;

      // Blob read
      if (peekval == FK_DELIN_BLB)
        ReadBlob(is, sectionTitle);
      else
        ReadKeyValue(is, GetSection(sectionTitle));

      // Re-peek
      peekval = (is >> std::ws).peek();
    }
  }

  void FKHeader::Print(std::ostream& out) const
  {
    out << SectionHeader("GridDesc", BLOB)<<std::endl;
    PrintBlob(out, "GridDesc");

    out << SectionHeader("VersionInfo", VERSIONS)<<std::endl;
    PrintKeyValue(out, VERSIONS);

    // Print nonstandard blobs
    keyMap::const_iterator iBlob = fBlobString.begin();
    for (;iBlob != fBlobString.end(); iBlob++)
      if ((*iBlob).first.compare("GridDesc") != 0)
      if ((*iBlob).first.compare("FlavourMap") != 0)
      if ((*iBlob).first.compare("xGrid") != 0)
        {
          out << SectionHeader((*iBlob).first.c_str(), BLOB)<<std::endl;
          PrintBlob(out, (*iBlob).first);
        }

    out << SectionHeader("GridInfo", GRIDINFO)<<std::endl;
    PrintKeyValue(out, GRIDINFO);

    out << SectionHeader("FlavourMap", BLOB)<<std::endl;
    PrintBlob(out, "FlavourMap");

    out << SectionHeader("TheoryInfo", THEORYINFO)<<std::endl;
    PrintKeyValue(out, THEORYINFO);
  
    out << SectionHeader("xGrid", BLOB)<<std::endl;
    PrintBlob(out, "xGrid");

    out << SectionHeader("FastKernel", BLOB)<<std::endl;
  }

  void FKHeader::AddTag( section sec, std::string const& key, std::string const& value)
  { 
    keyMap *tMap = GetMap(sec);
    keyMap::const_iterator iMap = (*tMap).find(key);
      
    if (iMap != (*tMap).end())
      throw FileError("FKHeader::AddTag","key clash: " + key);

    // trim trailing characters
    const size_t trimpos = std::min(value.size(), value.find_last_not_of(" \t\n\r")+1);
    (*tMap).insert(std::pair<std::string,std::string>(key,value.substr(0, trimpos)));
  }

  bool FKHeader::HasTag( section sec, std::string const& key) const
  {
    const keyMap *tMap = GetMap(sec);
    keyMap::const_iterator iMap = (*tMap).find(key);

    if (iMap != (*tMap).end())
      return true;
    return false;
  }

  std::string FKHeader::GetTag( section sec, std::string const& key) const
  { 
      const keyMap *tMap = GetMap(sec);
      keyMap::const_iterator iMap = (*tMap).find(key);
      
      if (iMap != (*tMap).end())
          return (*iMap).second;
      else
          throw FileError("FKHeader::GetTag","key " + key + " not found in header!");

      return std::string();
  }

  FKHeader::section FKHeader::GetSection(std::string const& title) const
  {
    if (title.compare("VersionInfo") == 0) return VERSIONS;
    if (title.compare("GridInfo") == 0) return GRIDINFO;
    if (title.compare("TheoryInfo") == 0) return THEORYINFO;
    throw FileError("FKHeader::GetSection","Unrecognised section title: " + title);
    return BLOB;
  }

  const FKHeader::keyMap* FKHeader::GetMap( section const& sec ) const
  {
    switch (sec)
    {
      case VERSIONS:       return &fVersions;
      case GRIDINFO:      return &fGridInfo;
      case THEORYINFO:    return &fTheoryInfo;
      case BLOB:          return &fBlobString;
    }

    return NULL;
  }

  void FKHeader::RemTag( section sec, std::string const& key )
  {
    keyMap *tMap = GetMap(sec);
    keyMap::iterator iTag = (*tMap).find(key);
    if (iTag != (*tMap).end())
      (*tMap).erase(iTag);
    else
      throw FileError("FKHeader::RemTag", "key " + key + " not found in header!");
  }


  void FKHeader::PrintKeyValue( std::ostream& os, section sec ) const
  {
    const keyMap *tMap = GetMap(sec);
    for (keyMap::const_iterator iPair = (*tMap).begin();
         iPair != (*tMap).end(); iPair++)
      os << "*" << (*iPair).first <<": "<<(*iPair).second<<std::endl;
  }

  void FKHeader::ReadKeyValue( std::istream& is, section sec )
  {
    // While the next char is a key-value pair indicator
    while ((is >> std::ws).peek() == FK_DELIN_KEY )
    {
      // Read Key-Value pair
      std::string keypair; getline( is, keypair );
      const std::string key = keypair.substr(1, keypair.find(':')-1);
      std::string val = keypair.substr(keypair.find(':')+1, keypair.size());

      // Remove leading spaces from value
      const size_t startpos = val.find_first_not_of(' ');
      val = val.substr(startpos, val.size());

      // Add the tag
      AddTag(sec, key, val);
    }

    return;

  }

  void FKHeader::PrintBlob(std::ostream& os, std::string blobname) const
  {
    keyMap::const_iterator iBlob = fBlobString.find(blobname);
    if (iBlob != fBlobString.end())
    { os << (*iBlob).second << std::endl; }
    else
      throw InitError("FKHeader::PrintBlob","Blob " + blobname + " not initialised");
  }

  void FKHeader::ReadBlob(std::istream& is, std::string blobname)
  {
     // While the next char is not a delineator
    std::stringstream blobstream;
    while ((is >> std::ws).peek() != FK_DELIN_SEC &&
           (is >> std::ws).peek() != FK_DELIN_BLB )
    {
      std::string blobline; getline( is, blobline );
      blobstream << blobline <<std::endl;
    }

    AddTag(BLOB, blobname, blobstream.str());
  }

  std::string FKHeader::SectionHeader(const char* title, section sec) const
  {
    std::string secdeln = sec == BLOB ? std::string("{"):std::string("_");
    std::string sechead = secdeln + title;
    sechead += std::string(60-sechead.size(),'_');
    return sechead;
  }


  // *********************************************************

  /**
   * @brief Constructor for FK Table
   * @param filename The FK table filename
   * @param cFactors A vector of filenames for potential C-factors
   */
  FKTable::FKTable( std::string const& filename, 
                    std::vector<std::string> const& cFactors):
  FKHeader(filename),
  fDataName(    GetTag        (GRIDINFO,       "SETNAME")),
  fDescription( GetTag        (BLOB,   "GridDesc")),
  fNData(       GetTag<int>   (GRIDINFO,   "NDATA")),
  fQ20(     pow(GetTag<double>(THEORYINFO, "Q0"),2)),
  fHadronic(    GetTag<bool>  (GRIDINFO,   "HADRONIC")),
  fNonZero(parseNonZero()),  // All flavours
  fFlmap(fHadronic ? new int[2*fNonZero]:new int[fNonZero]),
  fNx(          GetTag<int>   (GRIDINFO,   "NX")),
  fTx(fHadronic ? fNx*fNx:fNx),
  fRmr(fTx*fNonZero % convoluteAlign),
  fPad((fRmr == 0) ? 0:convoluteAlign - fRmr ),
  fDSz( fTx*fNonZero + fPad ),
  fXgrid(new double[fNx]),
  fSigma( new real[fDSz*fNData]),
  fHasCFactors(cFactors.size()),
  fcFactors(new double[fNData])
  {
    std::ifstream is(filename);
    NNPDF::FKHeader headSkip(is);
    InitialiseFromStream(is, cFactors);
  };

  /**
   * @brief Constructor for FK Table
   * @param is An input stream for reading from
   * @param cFactors A vector of filenames for potential C-factors
   */
  FKTable::FKTable( std::istream& is, std::vector<std::string> const& cFactors ):
  FKHeader(is),
  fDataName(    GetTag        (GRIDINFO,       "SETNAME")),
  fDescription( GetTag        (BLOB,   "GridDesc")),
  fNData(       GetTag<int>   (GRIDINFO,   "NDATA")),
  fQ20(     pow(GetTag<double>(THEORYINFO, "Q0"),2)),
  fHadronic(    GetTag<bool>  (GRIDINFO,   "HADRONIC")),
  fNonZero(parseNonZero()),  // All flavours
  fFlmap(fHadronic ? new int[2*fNonZero]:new int[fNonZero]),
  fNx(          GetTag<int>   (GRIDINFO,   "NX")),
  fTx(fHadronic ? fNx*fNx:fNx),
  fRmr(fTx*fNonZero % convoluteAlign),
  fPad((fRmr == 0) ? 0:convoluteAlign - fRmr ),
  fDSz( fTx*fNonZero + fPad ),
  fXgrid(new double[fNx]),
  fSigma( new real[fDSz*fNData]),
  fHasCFactors(cFactors.size()),
  fcFactors(new double[fNData])
  {
    InitialiseFromStream(is, cFactors);
  };

  /**
   * @brief Method for initialisation from stream
   * @param is the input stream after reading the FK header
   */
  void FKTable::InitialiseFromStream( std::istream& is, std::vector<std::string> const& cFactors )
  {
   if (Verbose)
    {
     std::cout <<SectionHeader(fDataName.c_str(),GRIDINFO)<<std::endl;
     std::cout << fDescription << std::endl;
    }

    // Read FlavourMap from header
    const int nFL = 14;
    std::stringstream fmBlob(GetTag(BLOB,"FlavourMap"));

    if (!fHadronic) // DIS flavourmap
    {
      bool maxMap[nFL]; // Maximal size flavourmap
      for (int i=0; i<nFL; i++)
        fmBlob >> maxMap[i];
     
      // Build FLMap
      int index = 0;
      for (int i=0; i<nFL; i++)
        if (maxMap[i]) { fFlmap[index] = i; index++; }
    }
    else // Hadronic flavourmap
    {
      bool maxMap[nFL][nFL]; // Maximal size flavourmap
      for (int i=0; i<nFL; i++)
        for (int j=0; j<nFL; j++)
          fmBlob >> maxMap[i][j];

      // Build FLMap
      int index = 0;
      for (int i=0; i<nFL; i++)
        for (int j=0; j<nFL; j++)
          if (maxMap[i][j])
          {
            fFlmap[2*index] = i;
            fFlmap[2*index+1] = j;
            index++;
          }
    }

    // Read x grid - need more detailed checks
    std::stringstream xBlob(GetTag(BLOB,"xGrid"));
    for (int i=0; i<fNx; i++)
      xBlob >> fXgrid[i];

    // Sanity checks
    if ( fNData <= 0 )
      throw RangeError("FKTable::FKTable","Number of datapoints is set to: " + std::to_string(fNData) );
    
    if ( fNx <= 0 )
      throw RangeError("FKTable::FKTable","Number of x-points is set to: " + std::to_string(fNx) );
    
    if ( fNonZero <= 0 )
      throw RangeError("FKTable::FKTable","Number of nonzero flavours is set to: " + std::to_string(fNonZero) );

    if (Verbose)
      std::cout << fNData << " Data Points "
      << fNx << " X points "
      << fNonZero << " active flavours"<<std::endl;


    // Zero sigma array -> also zeros pad quantities
    for (int i=0; i<fDSz*fNData; i++)
      fSigma[i]=0;

    // Read Cfactors
    for (int i = 0; i < fNData; i++) fcFactors[i] = 1.0;
      for (size_t i=0; i<cFactors.size(); i++)
       ReadCFactors(cFactors[i]);

    // Read FastKernel Table
    std::string line;
    std::vector<real> datasplit;

    if (fHadronic) {
      while (getline(is,line))
      {
        rsplit(datasplit,line);
        const int d = datasplit[0]; // datapoint
        const int a = datasplit[1]; // x1 index
        const int b = datasplit[2]; // x2 index
              
        for (int j=0; j<fNonZero; j++)
          {
            const int targetFlIndex = 14*fFlmap[2*j] + fFlmap[2*j+1]+3;
            fSigma[ d*fDSz+j*fTx+a*fNx+b ] = fcFactors[d]*datasplit[targetFlIndex];
          }
      }
    } else { // DIS
      while (getline(is,line))
      {
        rsplit(datasplit,line);
        
        const int d = datasplit[0];
        const int a = datasplit[1];
        
        for (int j=0; j<fNonZero; j++)
          fSigma[ d*fDSz+j*fNx+a ] = fcFactors[d]*datasplit[fFlmap[j]+2];

      }
    }  
  }



  /**
   * @brief FKTable copy-constructor.
   * @param set  The FK table to be copied
   */
  FKTable::FKTable(FKTable const& set):
  FKHeader(set),
  fDataName(set.fDataName),
  fNData(set.fNData),
  fQ20(set.fQ20),
  fHadronic(set.fHadronic),
  fNonZero(set.fNonZero),
  fFlmap(fHadronic ? (new int[2*fNonZero]):(new int[fNonZero])),
  fNx(set.fNx),
  fTx(set.fTx),
  fRmr(set.fRmr),
  fPad(set.fPad),
  fDSz(set.fDSz),
  fXgrid(new double[fNx]),
  fSigma(new real[fDSz*fNData]),
  fHasCFactors(set.fHasCFactors),
  fcFactors(new double[fNData])
  {
     // Copy X grid
    for (int i = 0; i < fNx; i++)
      fXgrid[i] = set.fXgrid[i];
      
    // Copy flavour map
    if (IsHadronic()) // Hadronic
    {
      for (int fl = 0; fl < fNonZero; fl++)
      {
        fFlmap[2*fl] = set.fFlmap[2*fl];
        fFlmap[2*fl+1] = set.fFlmap[2*fl+1];
      }
    }
    else  // DIS
    {
      for (int fl = 0; fl < fNonZero; fl++)
        fFlmap[fl] = set.fFlmap[fl];
    }
    
    // Copy reduced FK table
    for (int i = 0; i < fNData; i++)
      {
        for (int j = 0; j < fDSz; j++)
          fSigma[i*fDSz + j] = set.fSigma[i*fDSz + j];
        fcFactors[i] = set.GetCFactors()[i];
      }
  }

  /**
   * @brief Masked FK Table copy constructor.
   * @param set  The FK table to be filtered
   * @param mask The mask from filter
   */
  FKTable::FKTable(FKTable const& set, std::vector<int> const& mask):
  FKHeader(set),
  fDataName(set.fDataName),
  fNData(mask.size()),
  fQ20(set.fQ20),
  fHadronic(set.fHadronic),
  fNonZero(set.fNonZero),
  fFlmap(fHadronic ? (new int[2*fNonZero]):(new int[fNonZero])),
  fNx(set.fNx),
  fTx(set.fTx),
  fRmr(set.fRmr),
  fPad(set.fPad),
  fDSz(set.fDSz),
  fXgrid(new double[fNx]),
  fSigma(new real[fDSz*fNData]),
  fHasCFactors(set.fHasCFactors),
  fcFactors(new double[fNData])
  {
     if (fNData == 0)
       throw RangeError("FKTable::FKTable", "datapoints cut to 0!");

    // Copy X grid
    for (int i = 0; i < fNx; i++)
      fXgrid[i] = set.fXgrid[i];
      
    // Copy flavour map
    if (IsHadronic()) // Hadronic
    {
      for (int fl = 0; fl < fNonZero; fl++)
      {
        fFlmap[2*fl] = set.fFlmap[2*fl];
        fFlmap[2*fl+1] = set.fFlmap[2*fl+1];
      }
    } else  // DIS
    {
      for (int fl = 0; fl < fNonZero; fl++)
        fFlmap[fl] = set.fFlmap[fl];
    }
    
    // Copy reduced FK table
    for (int i = 0; i < fNData; i++)
      {
        for (int j = 0; j < fDSz; j++)
          fSigma[i*fDSz + j] = set.fSigma[mask[i]*fDSz + j];
        fcFactors[i] = set.GetCFactors()[mask[i]];
      }
  }

    /**
   * @brief FKTable destructor
   */
  FKTable::~FKTable()
  {
    delete[] fSigma;
    delete[] fFlmap;
    delete[] fXgrid;
    delete[] fcFactors;
  }


  /**
   * @brief FKTable print to ostream
   */
  void FKTable::Print(std::ostream& os)
  {
    if (Verbose)
      std::cout << "****** Exporting FKTable: "<<fDataName << " ******"<< std::endl;

    // Verify current flavourmap
    std::string nflmap;
    std::string cflmap = GetTag(BLOB, "FlavourMap");
    const bool optmap = OptimalFlavourmap(nflmap);
    if (!optmap)
    {
      RemTag(BLOB, "FlavourMap");
      AddTag(BLOB, "FlavourMap", nflmap);
    }

    // Print header
    FKHeader::Print(os);

    // Check that stream is still good
    if (!os.good())
      throw FileError("FKTable::Print","no good outstream!");

    if (fHasCFactors != 0)
    {
      std::cout << "FKTable::Print Warning: EXPORTING AN FKTABLE COMBINED WITH C-FACTORS" << std::endl;
      std::cout << "                        PLEASE ENSURE THAT THIS IS INTENTIONAL!" << std::endl;
    }

    // Write FastKernel Table
    if (fHadronic) 
    {
      for (int d=0; d<fNData; d++)
        for(int a=0; a<fNx; a++ )
          for (int b=0; b<fNx; b++)
          {
            bool isNonZero = false;
            std::stringstream outputline;
            for (int i=0; i<14; i++)
              for (int j=0; j<14; j++)
              {
                const int iSigma = GetISig(d,a,b,i,j);

                // Set precision
                if (iSigma == -1)
                  outputline << std::fixed <<std::setprecision(0);
                else if (fSigma[iSigma] == 0.0)
                  outputline << std::fixed <<std::setprecision(0);
                else
                  outputline << std::scientific << std::setprecision(16);

                outputline << (iSigma == -1 ? 0:fSigma[iSigma]) << "\t";

                if (iSigma != -1) if (fSigma[iSigma] != 0) isNonZero = true;
              }

            if (isNonZero) 
              os << d << "\t" << a <<"\t" << b <<"\t" <<  outputline.str() <<std::endl;
          }
    } else { // DIS
      for (int d=0; d<fNData; d++)
        for(int a=0; a<fNx; a++ )
          {
            bool isNonZero = false;
            std::stringstream outputline;
            for (int i=0; i<14; i++)
              {
                const int iSigma = GetISig(d,a,i);

                // Set precision
                if (iSigma == -1)
                  outputline << std::fixed <<std::setprecision(0);
                else if (fSigma[iSigma] == 0.0)
                  outputline << std::fixed <<std::setprecision(0);
                else
                  outputline << std::scientific << std::setprecision(16);

                outputline << (iSigma == -1 ? 0:fSigma[iSigma]) << "\t";
                if (iSigma != -1) if (fSigma[iSigma] != 0) isNonZero = true;
              }

            if (isNonZero)
              os << d << "\t" << a <<"\t" << outputline.str() <<std::endl;
          }
    }      

    // Restore current map
    if (!optmap)
    {
      RemTag(BLOB, "FlavourMap");
      AddTag(BLOB, "FlavourMap", cflmap);
    }
  
    return;    
  }

  
  void FKTable::ReadCFactors(std::string const& cfilename)
  {
    std::fstream g;
    g.open(cfilename.c_str(), std::ios::in);
    if (g.fail())
      throw FileError("FKTable::FKTable","cannot open cfactor file: " + cfilename);

    // Read through header
    std::string line;
    int nDelin = 0;
    while (nDelin < 2)
    {
      const int peekval = (g >> std::ws).peek();
      if (peekval == FK_DELIN_KEY)
        nDelin++;
      getline(g,line);
    }

    // Read C-factors
    getline(g,line); real tmp;
    for (int i = 0; i < fNData; i++)
    {
      g >> tmp;
      fcFactors[i] *= tmp;
    }

    g.close();
  }

  // Perform convolution
  void FKTable::Convolute(extern_pdf inpdf, size_t const& Npdf, real* out)
  {
    // Fetch PDF array
    real *pdf = 0;
    int err = posix_memalign(reinterpret_cast<void **>(&pdf), 16, sizeof(real)*fDSz*Npdf);
    if (err != 0)
      throw RangeError("FKTable::Convolute","ThPredictions posix_memalign " + std::to_string(err));
    CachePDF(inpdf, Npdf, pdf);

    // Calculate observables
    //#pragma omp parallel for
    for (int i = 0; i < fNData; i++)
      for (size_t n = 0; n < Npdf; n++)
      {
        out[i*Npdf + n] = 0;
        convolute(pdf+fDSz*n,fSigma+fDSz*i,out[i*Npdf + n],fDSz);
      }

    // Delete pdfs
    free(reinterpret_cast<void *>(pdf));
    return;
  }

    // Perform convolution
  void FKTable::CachePDF(extern_pdf inpdf, size_t const& NPDF, real* pdf)
  {
    // prepare PDF representation
    int index = 0;
    const int NFL = 14;
    real* EVLN = new real[fNx*NFL];
    
    for (size_t n = 0; n < NPDF; n++)
    {
      for (int i = 0; i < fNx; i++)
        inpdf(fXgrid[i], sqrt(fQ20), n, &EVLN[i*NFL]);
      
      for (int fl=0; fl<fNonZero; fl++)
      {
        const int fl1 = fFlmap[2*fl];
        const int fl2 = fFlmap[2*fl+1];
        
        for (int i = 0; i < fNx; i++)
          for (int j = 0; j < fNx; j++)
            pdf[index++] = EVLN[i*NFL+fl1]*EVLN[j*NFL+fl2];
      }
      
      for (int i=0; i<fPad; i++) // Add padding PDF
        pdf[index++] = 0;
      
    }
    
    delete[] EVLN;
    return;
  }


  // GetFKValue returns the appropriate point of the FK table
  int FKTable::GetISig( int const& d,     // Datapoint index
                        int const& a,   // First x-index
                        int const& b,   // Second x-index
                        int const& ifl1,  // First flavour index
                        int const& ifl2   // Second flavour index
                      ) const
  {

    if (!fHadronic)
      throw EvaluationError("FKTable::GetISig","Hadronic call for DIS table!");

    // Identify which nonZero flavour ifl1 and ifl2 correspond to
    int j = -1;
    for (int i=0; i < fNonZero; i++)
      if ( ifl1 == fFlmap[2*i] && ifl2 == fFlmap[2*i+1])
      {
        j = i;
        break;
      }

    // Not in FLmap
    if (j == -1)
      return -1;

    // Not in dataset
    if (d >= fNData)
      return -1;

    // Not in x-grid
    if (a >= fNx || b >= fNx)
      return -1;

    // Return pointer to FKTable segment
    return d*fDSz+j*fTx+a*fNx+b ;
  }

  // DIS version of GetFKValue
  int FKTable::GetISig(  int const& d,     // Datapoint index
                         int const& a,    // x-index
                         int const& ifl    // flavour index
                       ) const
  {

    if (fHadronic)
      throw EvaluationError("FKTable::GetISig","DIS call for Hadronic table!");

    // Identify which nonZero flavour ifl corresponds to
    int j = -1;
    for (int i=0; i < fNonZero; i++)
      if ( ifl == fFlmap[i] )
      {
        j = i;
        break;
      }

    // Not in FLmap
    if (j == -1)
      return -1;

    // Not in dataset
    if (d >= fNData)
      return -1;

    // Not in x-grid
    if (a >= fNx)
      return -1;

    return d*fDSz+j*fNx+a;
  }

  bool FKTable::OptimalFlavourmap(std::string& flmap) const  //!< Determine and return the optimal flavour map
  {
      bool nonZero[fNonZero];
      for (int j=0; j < fNonZero; j++)
      {
        nonZero[j] = false;
        for (int d=0; d < fNData; d++)
        {
          if (nonZero[j]) break;

          for (int a=0; a<fTx; a++)
            if (fSigma[d*fDSz+j*fTx+a] != 0)
            {
              nonZero[j] = true;
              break;
            }
        }
      }

      std::stringstream outmap;
      if (fHadronic)
      {
        for (int i=0; i<14; i++)
        {
          for (int j=0; j<14; j++)
          {
            bool found = false;
            for (int k=0; k<fNonZero; k++)
              if (nonZero[k] == true)
                if (fFlmap[2*k] == i && fFlmap[2*k+1] == j)
                  found = true;
            if (found) outmap << 1 <<" ";
            else outmap << 0 <<" ";
          }
          outmap << std::endl;
        }
      }
      else
      {
        for (int i=0; i<14; i++)
        {
          bool found = false;
          for (int k=0; k<fNonZero; k++)
            if (nonZero[k] == true)
              if (fFlmap[k] == i)
                found = true;
          if (found) outmap << 1 <<" ";
          else outmap << 0 <<" ";
        }
        
        outmap << std::endl;
      }

    flmap = outmap.str();

    bool anyZeros = false;
    for (int i=0; i<fNonZero; i++)
      if (nonZero[i] == false)
        anyZeros = true;

    return !anyZeros;
  }

  int FKTable::parseNonZero()
  {
    // Fetch the flavourmap blob
    std::stringstream fmBlob(GetTag(BLOB,"FlavourMap"));
    int nNonZero = 0;

    for (int i=0; i<14; i++)
    {

      if (!fmBlob.good())
        throw RangeError("FKTable::parseNonZero","FlavourMap formatting error!");

      if (!fHadronic)
      {
          bool iNonZero = false; fmBlob >> iNonZero;
          nNonZero += iNonZero;
      }
      else
      {
        for (int j=0; j<14; j++)
          if (fmBlob.good())
          {
            bool iNonZero = false; fmBlob >> iNonZero;
            nNonZero += iNonZero;
          }
          else
            throw RangeError("FKTable::parseNonZero","FlavourMap formatting error!");
      }
    }

    return nNonZero;
  }

}
