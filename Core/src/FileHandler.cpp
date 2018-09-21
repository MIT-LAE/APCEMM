/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* FileHandler Header File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : FileHandler.cpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "FileHandler.hpp"

static const int NC_ERR = 2;

FileHandler::FileHandler( )
{

    /* Default Constructor */

} /* End of FileHandler::FileHandler */

FileHandler::FileHandler( const char* fileName_, bool doRead_, bool doWrite_, bool overWrite_ )
{
    
    /* Constructor */

    fileName = fileName_;
    doRead = doRead_;
    doWrite = doWrite_;
    overWrite = overWrite_;
    
    isOpen = 0;

    NcError err( NcError::silent_nonfatal );

    NcFile::FileMode mode;

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


    NcFile dataFile = NcFile( fileName, mode );

    if ( dataFile.is_valid() ) {
        std::cout << " In FileHandler::FileHandler: Couldn't open: " << fileName << std::endl;
        std::cout << " -> doRead   : " << doRead << std::endl;
        std::cout << " -> doWrite  : " << doWrite << std::endl;
        std::cout << " -> overWrite: " << overWrite << std::endl;
    }

} /* End of FileHandler::FileHandler */

FileHandler::~FileHandler( )
{

    /* Destructor */

} /* End of FileHandler::~FileHandler */

FileHandler::FileHandler( const FileHandler &f )
{

    /* Copy file name*/ 
    fileName = f.fileName;

    /* Copy logicals */
    doRead = f.doRead;
    doWrite = f.doWrite;
    overWrite = f.overWrite;

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
    return *this;

} /* End of FileHandler::operator= */

const char* FileHandler::getFileName( ) const
{

    return fileName;

} /* End of FileHandler::getFileName */


/* End of FileHandler.cpp */
