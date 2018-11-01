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

static const int NC_ERROR = 2;
static const int NC_SUCCESS = 1;

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
        NcDim* addDim( NcFile &dataFile, const char* dimName, long size = 0 ) const;
        template <class T>
        int addConst( NcFile &dataFile, T *inputVar, const char* varName, long size, const char* type, const char* unit, const char* varFullName ) const;
        template <class T>
        inline int addConst( NcFile &dataFile, T *inputVar, const char* varName, long size, const char* type, const char* unit ) const {
            return addConst( dataFile, inputVar, varName, size, type, unit, varName );
        }
        template <class T>
        inline int addConst( NcFile &dataFile, T *inputVar, const char* varName, long size, const char* type) const {
            return addConst( dataFile, inputVar, varName, size, type, "-", varName );
        }
        template <class T>
        int addVar( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim, const char* type, const char* unit, const char* varFullName ) const;
        template <class T>
        inline int addVar( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim, const char* type, const char* unit ) const {
            return addVar( dataFile, inputVar, varName, varDim, type, unit, varName );
        }
        template <class T>
        inline int addVar( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim, const char* type ) const {
            return addVar( dataFile, inputVar, varName, varDim, type, "-", varName );
        }
        template <class T>
        int addVar2D( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const char* type, const char* unit, const char* varFullName ) const;
        template <class T>
        inline int addVar2D( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const char* type, const char* unit ) const {
            return addVar2D( dataFile, inputVar, varName, varDim1, varDim2, type, unit, varName );
        }
        template <class T>
        inline int addVar2D( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const char* type ) const {
            return addVar2D( dataFile, inputVar, varName, varDim1, varDim2, type, "-", varName );
        }
        template <class T>
        int addVar3D( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const NcDim *varDim3, const char* type, const char* unit, const char* varFullName ) const;
        template <class T>
        inline int addVar3D( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const NcDim *varDim3, const char* type, const char* unit ) const {
            return addVar3D( dataFile, inputVar, varName, varDim1, varDim2, varDim3, type, unit, varName );
        }
        template <class T>
        inline int addVar3D( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const NcDim *varDim3, const char* type ) const {
            return addVar3D( dataFile, inputVar, varName, varDim1, varDim2, varDim3, type, "-", varName );
        }
        template <class T>
        int addVar4D( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const NcDim *varDim3, const NcDim *varDim4, const char* type, const char* unit, const char* varFullName ) const;
        template <class T>
        inline int addVar4D( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const NcDim *varDim3, const NcDim *varDim4, const char* type, const char* unit ) const {
            return addVar4D( dataFile, inputVar, varName, varDim1, varDim2, varDim3, varDim4, type, unit, varName );
        }
        template <class T>
        inline int addVar4D( NcFile &dataFile, T *inputVar, const char* varName, const NcDim *varDim1, const NcDim *varDim2, const NcDim *varDim3, const NcDim *varDim4, const char* type ) const {
            return addVar4D( dataFile, inputVar, varName, varDim1, varDim2, varDim3, varDim4, type, "-", varName );
        }
        int addAtt( NcFile &dataFile, const char* attName, const char* attValue ) const;
        NcVar* getVar( NcFile &dataFile, const char* varName ) const; /* NcToken -> const char* */
        char* getAtt( NcVar* var ) const;
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

#endif /* FILEHANDLER_H_INCLUDED */
