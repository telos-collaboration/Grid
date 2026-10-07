/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/parallelIO/IldgIOtypes.h

Copyright (C) 2015

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */

/*!
 * \file
 * \brief Defines objects for storing Grid/Scidac/ILDG metadata.
 */

/// \cond DO_NOT_DOCUMENT
#ifndef GRID_ILDGTYPES_IO_H
#define GRID_ILDGTYPES_IO_H

#ifdef HAVE_LIME
extern "C" { // for linkage
#include "lime.h"
}

NAMESPACE_BEGIN(Grid);
/// \endcond

/*! 
 * \remark Data representation of records that enter ILDG and Scidac formats
 */
#define GRID_FORMAT               "grid-format"
#define ILDG_FORMAT               "ildg-format"
#define ILDG_BINARY_DATA          "ildg-binary-data"
#define ILDG_DATA_LFN             "ildg-data-lfn"
#define SCIDAC_CHECKSUM           "scidac-checksum"
#define SCIDAC_PRIVATE_FILE_XML   "scidac-private-file-xml"
#define SCIDAC_FILE_XML           "scidac-file-xml"
#define SCIDAC_PRIVATE_RECORD_XML "scidac-private-record-xml"
#define SCIDAC_RECORD_XML         "scidac-record-xml"
#define SCIDAC_BINARY_DATA        "scidac-binary-data"
// Unused SCIDAC records names; could move to support this functionality
#define SCIDAC_SITELIST           "scidac-sitelist"

////////////////////////////////////////////////////////////
const int GRID_IO_SINGLEFILE = 0; // hardcode lift from QIO compat
const int GRID_IO_MULTIFILE  = 1; // hardcode lift from QIO compat
const int GRID_IO_FIELD      = 0; // hardcode lift from QIO compat
const int GRID_IO_GLOBAL     = 1; // hardcode lift from QIO compat
////////////////////////////////////////////////////////////

/*! \remarks QIO uses mandatory "private" records fixed format.
 *  private is in principle "opaque" however it can't be changed 
 *  now because that would break existing file compatability, 
 *  so should be correct to assume the undocumented but defacto file structure.
 *  \brief required for Scidac file compatibility.
 */
struct emptyUserRecord : Serializable { 
  GRID_SERIALIZABLE_CLASS_MEMBERS(emptyUserRecord,int,dummy);
  emptyUserRecord() { dummy=0; };
};

/*!
 * Example metadata header: Scidac private file xml
\code{.xml}
<?xml version="1.0" encoding="UTF-8"?><scidacFile><version>1.1</version><spacetime>4</spacetime><dims>16 16 16 32 </dims><volfmt>0</volfmt></scidacFile>
\endcode
 */
/*! 
 *  \brief Serializable class to hold Scidac metadata.
 *  \param version Scidac lattice format version number.
 *  \param spacetime
 *  \param dims
 *  \param volfmt
 */
struct scidacFile : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(scidacFile,
                                  double, version,
                                  int, spacetime,
				  std::string, dims, // must convert to int
                                  int, volfmt);

  std::vector<int> getDimensions(void) { 
    std::stringstream stream(dims);
    std::vector<int> dimensions;
    int n;
    while(stream >> n){
      dimensions.push_back(n);
    }
    return dimensions;
  }

  void setDimensions(Coordinate dimensions) { 
    char delimiter = ' ';
    std::stringstream stream;
    for(int i=0;i<dimensions.size();i++){ 
      stream << dimensions[i];
      if ( i != dimensions.size()-1) { 
	stream << delimiter <<std::endl;
      }
    }
    dims = stream.str();
  }

  // Constructor provides Grid
  scidacFile() =default; // default constructor
  scidacFile(GridBase * grid){
    version      = 1.0;
    spacetime    = grid->_ndimension;
    setDimensions(grid->FullDimensions()); 
    volfmt       = GRID_IO_SINGLEFILE;
  }

};

/*!
 * Example metadata header: scidac-private-record-xml
\code{.xml}
<scidacRecord>
<version>1.1</version><date>Tue Jul 26 21:14:44 2011 UTC</date><recordtype>0</recordtype>
<datatype>QDP_D3_ColorMatrix</datatype><precision>D</precision><colors>3</colors><spins>4</spins>
<typesize>144</typesize><datacount>4</datacount>
</scidacRecord>
\endcode
 */
/*! 
 *  \brief Serializable class to hold Scidac metadata.
 *  \param version Scidac lattice format version number.
 *  \param date
 *  \param recordtype 
 *  \param datatype 
 *  \param precision 64 or 32.
 *  \param colors
 *  \param spins
 *  \param typesize 
 *  \param datacount 
 */
struct scidacRecord : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(scidacRecord,
                                  double, version,
                                  std::string, date,
				  int, recordtype,
				  std::string, datatype,
				  std::string, precision,
				  int, colors,
				  int, spins,
				  int, typesize,
				  int, datacount);

  scidacRecord()
  : version(1.0), recordtype(0), colors(0), spins(0), typesize(0), datacount(0)
  {}
};

/*! 
 *  \param version ILDG Binary Format Version.
 *  \param field   String representing gauge group and Nc.
 *  \param precision 64 or 32.
 *  \param rows Number of rows saved per gauge field. 
 *  \brief Serializable class to hold ILDG metadata.
 */
struct ildgFormat : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(ildgFormat,
				  double, version,
				  std::string, field,
				  int, precision,
				  int, rows,
				  int, lx,
				  int, ly,
				  int, lz,
				  int, lt);
  ildgFormat() { version=1.0; };
};

/*! 
 * \brief Serializable class to hold USQCD metadata.
 * \param version USQCD lattice format version
 * \param plaq average plaquette value
 * \param linktr average link trace value
 * \param info 
 */
struct usqcdInfo : Serializable { 
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(usqcdInfo,
				  double, version,
				  double, plaq,
				  double, linktr,
				  std::string, info);
  usqcdInfo() { 
    version=1.0; 
  };
};
/*!
 *  \brief Serializable Class for holding Scidac checksums.
 */
struct scidacChecksum : Serializable { 
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(scidacChecksum,
				  double, version,
				  std::string, suma,
				  std::string, sumb);
  scidacChecksum() { 
    version=1.0; 
  };
};
/// \cond DO_NOT_DOCUMENT
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Type:           scidac-file-xml         <title>MILC ILDG archival gauge configuration</title>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Type:           
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////
// Scidac private file xml 
// <?xml version="1.0" encoding="UTF-8"?><scidacFile><version>1.1</version><spacetime>4</spacetime><dims>16 16 16 32 </dims><volfmt>0</volfmt></scidacFile> 
////////////////////////                                                                                                                                                                              

#if 0
////////////////////////////////////////////////////////////////////////////////////////
// From http://www.physics.utah.edu/~detar/scidac/qio_2p3.pdf
////////////////////////////////////////////////////////////////////////////////////////
struct usqcdPropFile : Serializable { 
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(usqcdPropFile,
				  double, version,
				  std::string, type,
				  std::string, info);
  usqcdPropFile() { 
    version=1.0; 
  };
};
struct usqcdSourceInfo : Serializable { 
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(usqcdSourceInfo,
				  double, version,
				  std::string, info);
  usqcdSourceInfo() { 
    version=1.0; 
  };
};
struct usqcdPropInfo : Serializable { 
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(usqcdPropInfo,
				  double, version,
				  int, spin,
				  int, color,
				  std::string, info);
  usqcdPropInfo() { 
    version=1.0; 
  };
};
#endif

NAMESPACE_END(Grid);
/// \endcond

#endif
#endif
