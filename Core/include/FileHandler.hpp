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
#include <netcdfcpp.h>


class FileHandler
{

    public:

        FileHandler( );
        FileHandler( const char* fileName = NULL, bool doRead = 0, bool doWrite = 0, bool overWrite_ = 0 );
        ~FileHandler( );
        FileHandler( const FileHandler &f );
        FileHandler& operator=( const FileHandler &f );
        const char* getFileName( ) const;

    protected:

        /* File name */
        const char* fileName;

        /* Logicals */
        bool doRead;
        bool doWrite;
        bool overWrite;

        bool isOpen;

    private:

};

#endif /* FILEHANDLER_H_INCLUDED */

