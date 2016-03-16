// The MIT License (MIT)

// Copyright (c) Stefano Carrazza, Luigi Del Debbio, Nathan Hartland

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
#include <istream>
#include <fstream>
#include <sstream>

#include <map>

#include "common.h"

namespace NNPDF
{
    /**
    * Replacement for std::tostring
    */
    template <typename T>
    std::string ToString(T val)
    {
        std::stringstream stream;
        stream << val;
        return stream.str();
    }

    /**
    * Split std::string into std::vector
    */
    template<typename T>
    static void split(std::vector<T>& results, std::string const& input)
    {
        std::stringstream strstr(input);
        std::istream_iterator<T> it(strstr);
        std::istream_iterator<T> end;
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
    * \class FKHeader 
    * \brief Class for parsing the header of FK tables
    */
    class FKHeader
    {
    public:
        FKHeader();   //!< Prototype constructor
        FKHeader(std::istream&);                  //!< Construct from istream
        FKHeader(std::string const& filename);   //!< Construct from file
        FKHeader(FKHeader const&);               //!< Copy-construct
        ~FKHeader();                             //!< Destructor

        void Read(std::istream&);  //!< Read FKTable header from ostream
        void Print(std::ostream&) const; //!< Print FKTable header to ostream

        typedef std::map<std::string, std::string> keyMap;
        typedef enum {VERSIONS, GRIDINFO, THEORYINFO, BLOB} section;

        // ********************************* Tag Manip *********************************

        template<typename T>
        void AddTag( section sec, std::string const& key, T const& value)
        { AddTag(sec, key,ToString(value)); }
        void AddTag( section sec, std::string const& key, const char* value)
        { AddTag(sec, key,ToString(value));}
        void AddTag( section sec, std::string const& key, std::string const& value);

        // ********************************* Tag Getters *********************************

        bool HasTag( section sec, std::string const& key ) const; 
        std::string GetTag( section sec, std::string const& key) const;
        template<typename T> T GetTag( section sec, std::string const& key) const
        { return static_cast<T>(atof(GetTag(sec,key).c_str())); }

    protected:
        // Printing helper functions
        std::string SectionHeader(const char* title, section) const;

        void PrintKeyValue(std::ostream&, section) const;
        void ReadKeyValue( std::istream& is, section sec );

        void PrintBlob(std::ostream& os, std::string blobname) const;
        void ReadBlob(std::istream& is, std::string blobname);

        // Map fetchers helper functions
        section GetSection(std::string const&) const; //!< Fetch the appropriate section for a title
        const keyMap* GetMap(section const& sec) const; //!< Fetch the appropriate map for a section
        keyMap* GetMap(section const& sec)              //!< Fetch the appropriate map for a section
        { return const_cast<keyMap*>(const_cast<const FKHeader*>(this)->GetMap(sec)); };

        void RemTag( section sec, std::string const& key ); //!< Remove existing tags

        // ********************************* Attributes *********************************

        // Key-value pairs
        keyMap fVersions;
        keyMap fGridInfo;
        keyMap fTheoryInfo;
        keyMap fBlobString;
    };
}