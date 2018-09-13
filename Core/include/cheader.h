#ifndef CHEADER_H_INCLUDED
#define CHEADER_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

int KPP_Main( double varArray[], double fixArray[], double currentT, double dt, \
              double airDens, double temperature, double pressure, \
              double sinLAT, double cosLAT, double sinDEC, double cosDEC, \
              double rtols, double atols );

#ifdef __cplusplus
}
#endif


#endif /* CHEADER_H_INCLUDED */

