/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* FileHandler Header File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : FileHandler.hpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef FILEHANDLER_H_INCLUDED
#define FILEHANDLER_H_INCLUDED

#include <iostream>
#include <cstring>
#include <netcdfcpp.h>

template <typename T>
class FileHandler
{

    public:

        FileHandler( );
        FileHandler( const char* fileName, bool doWrite, bool doRead = 0, bool overWrite_ = 0 );
        ~FileHandler( );
        void UpdateMode();
        FileHandler( const FileHandler &f );
        FileHandler& operator=( const FileHandler &f );
        NcFile openFile( );
        void closeFile( NcFile &dataFile );
        int getNumDim( NcFile &dataFile ) const;
        int getNumVar( NcFile &dataFile ) const;
        int getNumAtt( NcFile &dataFile ) const;
        int addVar( NcFile &dataFile, T *inputVar, const char* varName, long size, const char* type, const char* unit = "" );
        NcVar* getVar( NcFile &dataFile, NcToken token ) const; /* NcToken -> const char* */
        const char* getFileName( ) const;
        bool isFileOpen( ) const;

        /* File name */
        const char* fileName;

        /* Logicals */
        bool doRead;
        bool doWrite;
        bool overWrite;

        /* File mode */
        NcFile::FileMode mode;

        bool isOpen;

    private:

};


static const int NC_ERROR = 0;
static const int NC_SUCCESS = 1;

template <typename T>
FileHandler<T>::FileHandler( )
    : fileName( NULL )
    , doRead( 0 )
    , doWrite( 0 )
    , overWrite( 0 )
{

    /* Default Constructor */

} /* End of FileHandler<T>::FileHandler */

template <typename T>
FileHandler<T>::FileHandler( const char* fileName_, bool doWrite_, bool doRead_, bool overWrite_ )
{
    
    /* Constructor */

    fileName = fileName_;
    doRead = doRead_;
    doWrite = doWrite_;
    overWrite = overWrite_;
    
    isOpen = 0;

    NcError err( NcError::silent_nonfatal );

} /* End of FileHandler<T>::FileHandler */

template <typename T>
FileHandler<T>::~FileHandler( )
{

    /* Destructor */

} /* End of FileHandler<T>::~FileHandler */

template <typename T>
void FileHandler<T>::UpdateMode( )
{

    switch ( doWrite ) {
        case 0:
            switch ( doRead ) {
                case 0:
                    /* Neither read nor write */
                    std::cout << " netCDF fileName: " << fileName << std::endl;
                    std::cout << " -> doRead and doWrite set to 0" << std::endl;
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

} /* End of FileHandler<T>::UpdateMode */

template <typename T>
FileHandler<T>::FileHandler( const FileHandler &f )
{

    /* Copy file name*/ 
    fileName = f.fileName;

    /* Copy logicals */
    doRead = f.doRead;
    doWrite = f.doWrite;
    overWrite = f.overWrite;
  
    /* Copy file */
    isOpen = f.isOpen;

} /* End of FileHandler<T>::FileHandler */

template <typename T>
FileHandler<T>& FileHandler<T>::operator=( const FileHandler &f )
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

} /* End of FileHandler<T>::operator= */

template <typename T>
NcFile FileHandler<T>::openFile( )
{

    NcFile dataFile( fileName, mode );

    isOpen = 1;
    
    if ( !( dataFile.is_valid() ) ) {
        std::cout << " In FileHandler<T>::FileHandler: Couldn't open: " << fileName << std::endl;
        std::cout << " -> doRead   : " << doRead << std::endl;
        std::cout << " -> doWrite  : " << doWrite << std::endl;
        std::cout << " -> overWrite: " << overWrite << std::endl;
        isOpen = 0;
    }
    return dataFile;

} /* End of FileHandler<T>::openFile */

template <typename T>
void FileHandler<T>::closeFile( NcFile &dataFile )
{

//    dataFile.close( );
    isOpen = 0;

    if ( !(dataFile.close()) ) {
        std::cout << "In FileHandler<T>::closeFile: Couldn't close " << fileName << std::endl;
        std::cout << "In FileHandler<T>::closeFile: is_valid returns: " << dataFile.is_valid() << std::endl;
        isOpen = 1;
    }

} /* End of FileHandler<T>::closeFile */

template <typename T>
int FileHandler<T>::getNumDim( NcFile &dataFile ) const
{

    return dataFile.num_dims( );

} /* End of FileHandler<T>::getNumDim */

template <typename T>
int FileHandler<T>::getNumVar( NcFile &dataFile ) const
{

    return dataFile.num_vars( );

} /* End of FileHandler<T>::getNumVar */

template <typename T>
int FileHandler<T>::getNumAtt( NcFile &dataFile ) const
{

    return dataFile.num_atts( );

} /* End of FileHandler<T>::getNumAtt */

template <typename T>
int FileHandler<T>::addVar( NcFile &dataFile, T *inputVar, const char* varName, long size, const char* type, const char* unit )
{

    /* Define the dimensions. netCDF will hand back an NcDim object. */
    NcDim *varDim;
    if ( !(varDim = dataFile.add_dim( varName, size )) ) {
        std::cout << "In FileHandler<T>::addVar: defining dimension failed for " << varName << " in " << fileName << std::endl;
        return NC_ERROR;
    }

    /* Convert input type to netCDF type */
    NcType varType;
    if ( strcmp(type, "double") ) {
        varType = ncDouble;
    } else if ( strcmp(type, "float") ) {
        varType = ncFloat;
    } else if ( strcmp(type, "long") ) {
        varType = ncLong;
    } else if ( strcmp(type, "int") ) {
        varType = ncInt;
    } else if ( strcmp(type, "short") ) {
        varType = ncShort;
    } else if ( strcmp(type, "char") ) {
        varType = ncChar;
    } else if ( strcmp(type, "byte") ) {
        varType = ncByte;
    } else {
        varType = ncNoType;
        std::cout << "In FileHandler<T>::addVar: varType takes an undefined value. varType: " << varType << " in " << fileName << std::endl;
        return NC_ERROR;
    }

    /* Create netCDF variables which hold the actual specified variable */
    NcVar *var;
    if ( !(var = dataFile.add_var( varName, varType, varDim )) ) {
        std::cout << "In FileHandler<T>::addVar: defining variable failed for " << varName << " in " << fileName << std::endl;
        return NC_ERROR;
    }

    /* Define units attributes for variables. This attaches a text attribute to each of the coordinate variables, containing the units */
    if ( !var->add_att("units", unit) ) {
        std::cout << "In FileHandler<T>::addVar: unit definition failed for " << varName << " ( unit: [" << unit << "]) in " << fileName << std::endl;  
        return NC_ERROR;
    }

    /* Write the variable data. The arrays of data are the same size as the netCDF variables we have defined, and below we write it in one step */
    if ( !var->put(inputVar, size) ) {
        std::cout << "In FileHandler<T>::addVar: writing variable failed for " << varName << " in " << fileName << std::endl;  
        return NC_ERROR;
    }

    std::cout << " Variable " << varName << " has been written to " << fileName << "!" << std::endl;
    return NC_SUCCESS;

}

template <typename T>
NcVar* FileHandler<T>::getVar( NcFile &dataFile, NcToken token ) const
{

    return dataFile.get_var( token );

} /* End of FileHandler<T>::getVar */

template <typename T>
const char* FileHandler<T>::getFileName( ) const
{

    return fileName;

} /* End of FileHandler<T>::getFileName */

template <typename T>
bool FileHandler<T>::isFileOpen( ) const
{

    return isOpen;

} /* End of FileHandler<T>::isOpen */



#endif /* FILEHANDLER_H_INCLUDED */

