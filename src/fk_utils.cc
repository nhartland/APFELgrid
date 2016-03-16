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

#include "APFELgrid/fk_utils.h"
#include "APFELgrid/exceptions.h"

namespace NNPDF
{

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

}
