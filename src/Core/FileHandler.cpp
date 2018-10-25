/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* FileHandler Program File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : FileHandler.cpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/FileHandler.hpp"
    
FileHandler::FileHandler( )
    : fileName( NULL )
    , doRead( 0 )
    , doWrite( 0 )
    , overWrite( 0 )
{

    /* Default Constructor */

} /* End of FileHandler::FileHandler */


FileHandler::FileHandler( const char* fileName_, bool doWrite_, bool doRead_, bool overWrite_ )
{
    
    /* Constructor */

    fileName = fileName_;
    doRead = doRead_;
    doWrite = doWrite_;
    overWrite = overWrite_;
    
    isOpen = 0;

    NcError err( NcError::silent_nonfatal );

    UpdateMode();

} /* End of FileHandler::FileHandler */


FileHandler::~FileHandler( )
{

    /* Destructor */

} /* End of FileHandler::~FileHandler */


void FileHandler::UpdateMode( )
{

    switch ( doWrite ) {
        case 0:
            switch ( doRead ) {
                case 0:
                    /* Neither read nor write */
                    std::cout << " netCDF fileName: " << fileName << "\n";
                    std::cout << " -> doRead and doWrite set to 0" << "\n";
                case 1:
                    /* Read only */
                    mode = NcFile::ReadOnly;
            }
        case 1:
            switch ( overWrite ) {
                case 0:
                    /* Write without overwriting */
                    mode = NcFile::New;
                case 1:
                    /* Write with overwriting */
                    mode = NcFile::Replace;
            }
    }

} /* End of FileHandler::UpdateMode */


FileHandler::FileHandler( const FileHandler &f )
{

    /* Copy file name*/ 
    fileName = f.fileName;

    /* Copy logicals */
    doRead = f.doRead;
    doWrite = f.doWrite;
    overWrite = f.overWrite;
  
    /* Copy file */
    isOpen = f.isOpen;

} /* End of FileHandler::FileHandler */


FileHandler& FileHandler::operator=( const FileHandler &f )
{

    if ( &f == this )
        return *this;
    
    /* Assign file name*/ 
    fileName = f.fileName;

    /* Assign logicals */
    doRead = f.doRead;
    doWrite = f.doWrite;
    overWrite = f.overWrite;

    /* Assign file */
    isOpen = f.isOpen;
    return *this;

} /* End of FileHandler::operator= */


NcFile FileHandler::openFile( )
{

    NcFile dataFile( fileName, mode );

    isOpen = 1;
    
    if ( !( dataFile.is_valid() ) ) {
        std::cout << " In FileHandler::FileHandler: Couldn't open: " << fileName << "\n";
        std::cout << " -> doRead   : " << doRead << "\n";
        std::cout << " -> doWrite  : " << doWrite << "\n";
        std::cout << " -> overWrite: " << overWrite << "\n";
        std::cout << " -> mode: " << mode << "\n";
        isOpen = 0;
    }
    return dataFile;

} /* End of FileHandler::openFile */


void FileHandler::closeFile( NcFile &dataFile )
{

    isOpen = 0;

    if ( !(dataFile.close()) ) {
        std::cout << "In FileHandler::closeFile: Couldn't close " << fileName << "\n";
        if ( dataFile.is_valid() ) {
            std::cout << "In FileHandler::closeFile: is_valid returns: " << dataFile.is_valid() << "\n";
            isOpen = 1;
        }
    }

} /* End of FileHandler::closeFile */


int FileHandler::getNumDim( NcFile &dataFile ) const
{

    return dataFile.num_dims( );

} /* End of FileHandler::getNumDim */


int FileHandler::getNumVar( NcFile &dataFile ) const
{

    return dataFile.num_vars( );

} /* End of FileHandler::getNumVar */


int FileHandler::getNumAtt( NcFile &dataFile ) const
{

    return dataFile.num_atts( );

} /* End of FileHandler::getNumAtt */

NcDim* FileHandler::addDim( NcFile &dataFile, const char* dimName, long size ) const
{

    /* If size = 0, dimension is unlimited */

    NcDim *varDim;
    /* Define the dimensions. netCDF will hand back an NcDim object. */
    if ( !(varDim = dataFile.add_dim( dimName, size )) ) {
        std::cout << "In FileHandler::addDim: defining dimension failed for " << dimName << "\n";
    }

    return varDim;

} /* End of FileHandler::addDim */

template <class T>
int FileHandler::addConst( NcFile &dataFile, T *inputVar, const char* varName, long size, const char* type, const char* unit, const char* varFullName ) const
{

    /* Convert input type to netCDF type */
    NcType varType;
    if ( strcmp(type, "double") == 0 ) {
        varType = ncDouble;
    } else if ( strcmp(type, "float") == 0 ) {
        varType = ncFloat;
    } else if ( strcmp(type, "long") == 0 ) {
        varType = ncLong;
    } else if ( strcmp(type, "int") == 0 ) {
        varType = ncInt;
    } else if ( strcmp(type, "short") == 0 ) {
        varType = ncShort;
    } else if ( strcmp(type, "char") == 0 ) {
        varType = ncChar;
    } else if ( strcmp(type, "byte") == 0 ) {
        varType = ncByte;
    } else {
        varType = ncNoType;
        std::cout << "In FileHandler::addConst: varType takes an undefined value. varType: " << varType << " in " << fileName << "\n";
        return NC_ERROR;
    }

    /* Create netCDF variables which hold the actual specified variable */
    NcVar *var;
    if ( !(var = dataFile.add_var( varName, varType )) ) {
        std::cout << "In FileHandler::addConst: defining variable failed for " << varName << " in " << fileName << "\n";
        return NC_ERROR;
    }

    /* Define unit attributes for variables. This attaches a text attribute to each of the coordinate variables, containing the unit */
    if ( !var->add_att("unit", unit) ) {
        std::cout << "In FileHandler::addConst: unit definition failed for " << varName << " ( unit: [" << unit << "]) in " << fileName << "\n";
        return NC_ERROR;
    }
    
    /* Define long name attributes for variables. This attaches a text attribute to each of the coordinate variables, containing the long name */
    if ( !var->add_att("long_name", varFullName) ) {
        std::cout << "In FileHandler::addConst: full name definition failed for " << varName << " ( full name: [" << varFullName << "]) in " << fileName << "\n";
        return NC_ERROR;
    }

    /* Write the variable data. The arrays of data are the same size as the netCDF variables we have defined, and below we write it in one step */
    if ( !var->put(inputVar, size) ) {
        std::cout << "In FileHandler::addConst: writing variable failed for " << varName << " in " << fileName << "\n";
        return NC_ERROR;
    }

    std::cout << " -> Variable " << varName << " has been written to " << fileName << "!" << "\n";
    return NC_SUCCESS;

} /* End of FileHandler::addConst */

template <class T>
int FileHandler::addVar( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim, const char* type, const char* unit, const char* varFullName ) const
{

    /* Convert input type to netCDF type */
    NcType varType;
    if ( strcmp(type, "double") == 0 ) {
        varType = ncDouble;
    } else if ( strcmp(type, "float") == 0 ) {
        varType = ncFloat;
    } else if ( strcmp(type, "long") == 0 ) {
        varType = ncLong;
    } else if ( ( strcmp(type, "int") == 0 ) || ( strcmp(type, "unsigned int") == 0 ) ) {
        varType = ncInt;
    } else if ( strcmp(type, "short") == 0 ) {
        varType = ncShort;
    } else if ( strcmp(type, "char") == 0 ) {
        varType = ncChar;
    } else if ( strcmp(type, "byte") == 0 ) {
        varType = ncByte;
    } else {
        varType = ncNoType;
        std::cout << "In FileHandler::addVar: varType takes an undefined value. varType: " << varType << " in " << fileName << "\n";
        return NC_ERROR;
    }

    /* Create netCDF variables which hold the actual specified variable */
    NcVar *var;
    if ( !(var = dataFile.add_var( varName, varType, varDim )) ) {
        std::cout << "In FileHandler::addVar: defining variable failed for " << varName << " in " << fileName << "\n";
        return NC_ERROR;
    }

    /* Define unit attributes for variables. This attaches a text attribute to each of the coordinate variables, containing the unit */
    if ( !var->add_att("unit", unit) ) {
        std::cout << "In FileHandler::addVar: unit definition failed for " << varName << " ( unit: [" << unit << "]) in " << fileName << "\n";
        return NC_ERROR;
    }

    /* Define long name attributes for variables. This attaches a text attribute to each of the coordinate variables, containing the long name */
    if ( !var->add_att("long_name", varFullName) ) {
        std::cout << "In FileHandler::addConst: full name definition failed for " << varName << " ( full name: [" << varFullName << "]) in " << fileName << "\n";
        return NC_ERROR;
    }
    
    /* Write the variable data. The arrays of data are the same size as the netCDF variables we have defined, and below we write it in one step */
    if ( !var->put( inputVar, varDim->size() ) ) {
        std::cout << "In FileHandler::addVar: writing variable failed for " << varName << " in " << fileName << "\n";
        return NC_ERROR;
    }

    std::cout << " -> Variable " << varName << " has been written to " << fileName << "!" << "\n";
    return NC_SUCCESS;

} /* End of FileHandler::addVar */

template <class T>
int FileHandler::addVar2D( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const char* type, const char* unit, const char* varFullName ) const
{

    /* Convert input type to netCDF type */
    NcType varType;
    if ( strcmp(type, "double") == 0 ) {
        varType = ncDouble;
    } else if ( strcmp(type, "float") == 0 ) {
        varType = ncFloat;
    } else if ( strcmp(type, "long") == 0 ) {
        varType = ncLong;
    } else if ( strcmp(type, "int") == 0 ) {
        varType = ncInt;
    } else if ( strcmp(type, "short") == 0 ) {
        varType = ncShort;
    } else if ( strcmp(type, "char") == 0 ) {
        varType = ncChar;
    } else if ( strcmp(type, "byte") == 0 ) {
        varType = ncByte;
    } else {
        varType = ncNoType;
        std::cout << "In FileHandler::addVar2D: varType takes an undefined value. varType: " << varType << " in " << fileName << "\n";
        return NC_ERROR;
    }

    /* Create netCDF variables which hold the actual specified variable */
    NcVar *var;
    if ( !(var = dataFile.add_var( varName, varType, varDim1, varDim2 )) ) {
        std::cout << "In FileHandler::addVar2D: defining variable failed for " << varName << " in " << fileName << "\n";
        return NC_ERROR;
    }

    /* Define unit attributes for variables. This attaches a text attribute to each of the coordinate variables, containing the unit */
    if ( !var->add_att("unit", unit) ) {
        std::cout << "In FileHandler::addVar2D: unit definition failed for " << varName << " ( unit: [" << unit << "]) in " << fileName << "\n";
        return NC_ERROR;
    }

    /* Define long name attributes for variables. This attaches a text attribute to each of the coordinate variables, containing the long name */
    if ( !var->add_att("long_name", varFullName) ) {
        std::cout << "In FileHandler::addConst: full name definition failed for " << varName << " ( full name: [" << varFullName << "]) in " << fileName << "\n";
        return NC_ERROR;
    }
    
    /* Write the variable data. The arrays of data are the same size as the netCDF variables we have defined, and below we write it in one step */
    if ( !var->put( inputVar, varDim1->size(), varDim2->size() ) ) {
        std::cout << "In FileHandler::addVar2D: writing variable failed for " << varName << " in " << fileName << "\n";
        return NC_ERROR;
    }

    std::cout << " -> Variable " << varName << " has been written to " << fileName << "!" << "\n";
    return NC_SUCCESS;

} /* End of FileHandler::addVar2D */

int FileHandler::addAtt( NcFile &dataFile, const char* attName, const char* attValue ) const
{

    if ( !dataFile.add_att( attName, attValue ) ) {
        std::cout << "In FileHandler::addAtt: adding attribute " << attName << " with value " << attValue << "failed! \n";
        return NC_ERROR;
    }

    return NC_SUCCESS;

} /* End of FileHandler::addAtt */

NcVar* FileHandler::getVar( NcFile &dataFile, const char* varName ) const
{

    NcVar* var;

    /* Give back a pointer to the requested NcVar */
    if ( !(var = dataFile.get_var( varName )) ) {
        std::cout << "In FileHandler::getVar: getting variable " << varName << " failed for " << fileName << "\n";
    }
        
    return var;

} /* End of FileHandler::getVar */


char* FileHandler::getAtt( NcVar* var ) const
{

    /* Each of the netCDF variables has a "unit" attribute */
    NcAtt *att;
    char *unit;

    if ( !(att = var->get_att("unit")) ) {
        std::cout << "In FileHandler::getAtt: getting atribute 'unit' failed in " << fileName << "\n";
    }
    unit = att->as_string(0);
    delete att;

    return unit;

} /* End of FileHandler::getAtt */


const char* FileHandler::getFileName( ) const
{

    return fileName;

} /* End of FileHandler::getFileName */


bool FileHandler::isFileOpen( ) const
{

    return isOpen;

} /* End of FileHandler::isOpen */

template int FileHandler::addVar<double>( NcFile &dataFile, double *inputVar, const char* varName, const NcDim *varDim, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addVar<int>( NcFile &dataFile, int *inputVar, const char* varName, const NcDim *varDim, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addVar<float>( NcFile &dataFile, float *inputVar, const char* varName, const NcDim *varDim, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addConst<double>( NcFile &dataFile, double *inputVar, const char* varName, long size, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addConst<int>( NcFile &dataFile, int *inputVar, const char* varName, long size, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addConst<float>( NcFile &dataFile, float *inputVar, const char* varName, long size, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addVar2D<double>( NcFile &dataFile, double *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addVar2D<int>( NcFile &dataFile, int *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addVar2D<float>( NcFile &dataFile, float *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addVar<const double>( NcFile &dataFile, const double *inputVar, const char* varName, const NcDim *varDim, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addVar<const int>( NcFile &dataFile, const int *inputVar, const char* varName, const NcDim *varDim, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addVar<const float>( NcFile &dataFile, const float *inputVar, const char* varName, const NcDim *varDim, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addConst<const double>( NcFile &dataFile, const double *inputVar, const char* varName, long size, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addConst<const int>( NcFile &dataFile, const int *inputVar, const char* varName, long size, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addConst<const float>( NcFile &dataFile, const float *inputVar, const char* varName, long size, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addVar2D<const double>( NcFile &dataFile, const double *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addVar2D<const int>( NcFile &dataFile, const int *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const char* type, const char* unit, const char* varFullName ) const;
template int FileHandler::addVar2D<const float>( NcFile &dataFile, const float *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const char* type, const char* unit, const char* varFullName ) const;

/* End of FileHandler.cpp */
