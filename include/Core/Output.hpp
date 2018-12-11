/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Output Header File                                               */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 8/12/2018                                 */
/* File                 : Output.hpp                                */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef OUTPUT_H_INCLUDED
#define OUTPUT_H_INCLUDED

#define DO_SAVE_NOx                     1 /* Save NOx      */
#define DO_SAVE_NOy                     1 /* Save NOy      */
#define DO_SAVE_CO2                     1 /* Save CO2      */
#define DO_SAVE_PPN                     0 /* Save PPN      */
#define DO_SAVE_BrNO2                   0 /* Save BrNO2    */
#define DO_SAVE_IEPOX                   0 /* Save IEPOX    */
#define DO_SAVE_PMNN                    0 /* Save PMNN     */
#define DO_SAVE_N2O                     0 /* Save N2O      */
#define DO_SAVE_N                       0 /* Save N        */
#define DO_SAVE_PAN                     1 /* Save PAN      */
#define DO_SAVE_ALK4                    0 /* Save ALK4     */
#define DO_SAVE_MAP                     0 /* Save MAP      */
#define DO_SAVE_MPN                     1 /* Save MPN      */
#define DO_SAVE_Cl2O2                   0 /* Save Cl2O2    */
#define DO_SAVE_ETP                     0 /* Save ETP      */
#define DO_SAVE_HNO2                    1 /* Save HNO2     */
#define DO_SAVE_C3H8                    0 /* Save C3H8     */
#define DO_SAVE_RA3P                    0 /* Save RA3P     */
#define DO_SAVE_RB3P                    0 /* Save RB3P     */
#define DO_SAVE_OClO                    0 /* Save OClO     */
#define DO_SAVE_ClNO2                   0 /* Save ClNO2    */
#define DO_SAVE_ISOP                    0 /* Save ISOP     */
#define DO_SAVE_HNO4                    1 /* Save HNO4     */
#define DO_SAVE_MAOP                    0 /* Save MAOP     */
#define DO_SAVE_MP                      0 /* Save MP       */
#define DO_SAVE_ClOO                    0 /* Save ClOO     */
#define DO_SAVE_RP                      0 /* Save RP       */
#define DO_SAVE_BrCl                    0 /* Save BrCl     */
#define DO_SAVE_PP                      0 /* Save PP       */
#define DO_SAVE_PRPN                    0 /* Save PRPN     */
#define DO_SAVE_SO4                     0 /* Save SO4      */
#define DO_SAVE_SO4_L                   0 /* Save SO4_L    */
#define DO_SAVE_Br2                     0 /* Save Br2      */
#define DO_SAVE_ETHLN                   0 /* Save ETHLN    */
#define DO_SAVE_MVKN                    0 /* Save MVKN     */
#define DO_SAVE_R4P                     0 /* Save R4P      */
#define DO_SAVE_C2H6                    0 /* Save C2H6     */
#define DO_SAVE_RIP                     0 /* Save RIP      */
#define DO_SAVE_VRP                     0 /* Save VRP      */
#define DO_SAVE_ATOOH                   0 /* Save ATOOH    */
#define DO_SAVE_IAP                     0 /* Save IAP      */
#define DO_SAVE_DHMOB                   0 /* Save DHMOB    */
#define DO_SAVE_MOBA                    0 /* Save MOBA     */
#define DO_SAVE_MRP                     0 /* Save MRP      */
#define DO_SAVE_N2O5                    1 /* Save N2O5     */
#define DO_SAVE_ISNOHOO                 0 /* Save ISNOHOO  */
#define DO_SAVE_ISNP                    0 /* Save ISNP     */
#define DO_SAVE_ISOPNB                  0 /* Save ISOPNB   */
#define DO_SAVE_IEPOXOO                 0 /* Save IEPOXOO  */
#define DO_SAVE_MACRNO2                 0 /* Save MACRNO2  */
#define DO_SAVE_ROH                     0 /* Save ROH      */
#define DO_SAVE_MOBAOO                  0 /* Save MOBAOO   */
#define DO_SAVE_DIBOO                   0 /* Save DIBOO    */
#define DO_SAVE_PMN                     0 /* Save PMN      */
#define DO_SAVE_ISNOOB                  0 /* Save ISNOOB   */
#define DO_SAVE_INPN                    0 /* Save INPN     */
#define DO_SAVE_H                       0 /* Save H        */
#define DO_SAVE_BrNO3                   0 /* Save BrNO3    */
#define DO_SAVE_PRPE                    0 /* Save PRPE     */
#define DO_SAVE_MVKOO                   0 /* Save MVKOO    */
#define DO_SAVE_Cl2                     0 /* Save Cl2      */
#define DO_SAVE_ISOPND                  0 /* Save ISOPND   */
#define DO_SAVE_HOBr                    0 /* Save HOBr     */
#define DO_SAVE_HOBr_L                  0 /* Save HOBr_L   */
#define DO_SAVE_A3O2                    0 /* Save A3O2     */
#define DO_SAVE_PROPNN                  0 /* Save PROPNN   */
#define DO_SAVE_GLYX                    0 /* Save GLYX     */
#define DO_SAVE_MAOPO2                  0 /* Save MAOPO2   */
#define DO_SAVE_CH4                     1 /* Save CH4      */
#define DO_SAVE_GAOO                    0 /* Save GAOO     */
#define DO_SAVE_B3O2                    0 /* Save B3O2     */
#define DO_SAVE_ACET                    0 /* Save ACET     */
#define DO_SAVE_MACRN                   0 /* Save MACRN    */
#define DO_SAVE_CH2OO                   0 /* Save CH2OO    */
#define DO_SAVE_MGLYOO                  0 /* Save MGLYOO   */
#define DO_SAVE_VRO2                    0 /* Save VRO2     */
#define DO_SAVE_MGLOO                   0 /* Save MGLOO    */
#define DO_SAVE_MACROO                  0 /* Save MACROO   */
#define DO_SAVE_PO2                     0 /* Save PO2      */
#define DO_SAVE_CH3CHOO                 0 /* Save CH3CHOO  */
#define DO_SAVE_MAN2                    0 /* Save MAN2     */
#define DO_SAVE_ISNOOA                  0 /* Save ISNOOA   */
#define DO_SAVE_H2O2                    1 /* Save H2O2     */
#define DO_SAVE_PRN1                    0 /* Save PRN1     */
#define DO_SAVE_ETO2                    0 /* Save ETO2     */
#define DO_SAVE_KO2                     0 /* Save KO2      */
#define DO_SAVE_RCO3                    0 /* Save RCO3     */
#define DO_SAVE_HC5OO                   0 /* Save HC5OO    */
#define DO_SAVE_GLYC                    0 /* Save GLYC     */
#define DO_SAVE_ClNO3                   0 /* Save ClNO3    */
#define DO_SAVE_RIO2                    0 /* Save RIO2     */
#define DO_SAVE_R4N1                    0 /* Save R4N1     */
#define DO_SAVE_HOCl                    0 /* Save HOCl     */
#define DO_SAVE_HOCl_L                  0 /* Save HOCl_L   */
#define DO_SAVE_ATO2                    0 /* Save ATO2     */
#define DO_SAVE_HNO3                    1 /* Save HNO3     */
#define DO_SAVE_HNO3_S                  0 /* Save HNO3_S   */
#define DO_SAVE_HNO3_L                  0 /* Save HNO3_L   */
#define DO_SAVE_ISN1                    0 /* Save ISN1     */
#define DO_SAVE_MAO3                    0 /* Save MAO3     */
#define DO_SAVE_MRO2                    0 /* Save MRO2     */
#define DO_SAVE_INO2                    0 /* Save INO2     */
#define DO_SAVE_HAC                     0 /* Save HAC      */
#define DO_SAVE_HC5                     0 /* Save HC5      */
#define DO_SAVE_MGLY                    0 /* Save MGLY     */
#define DO_SAVE_ISOPNBO2                0 /* Save ISOPNBO2 */
#define DO_SAVE_ISOPNDO2                0 /* Save ISOPNDO2 */
#define DO_SAVE_R4O2                    0 /* Save R4O2     */
#define DO_SAVE_R4N2                    0 /* Save R4N2     */
#define DO_SAVE_BrO                     0 /* Save BrO      */
#define DO_SAVE_RCHO                    0 /* Save RCHO     */
#define DO_SAVE_MEK                     0 /* Save MEK      */
#define DO_SAVE_ClO                     0 /* Save ClO      */
#define DO_SAVE_MACR                    0 /* Save MACR     */
#define DO_SAVE_SO2                     1 /* Save SO2      */
#define DO_SAVE_MVK                     0 /* Save MVK      */
#define DO_SAVE_ALD2                    0 /* Save ALD2     */
#define DO_SAVE_MCO3                    0 /* Save MCO3     */
#define DO_SAVE_CH2O                    0 /* Save CH2O     */
#define DO_SAVE_H2O                     1 /* Save H2O      */
#define DO_SAVE_H2O_S                   0 /* Save H2O_S    */
#define DO_SAVE_H2O_L                   0 /* Save H2O_L    */
#define DO_SAVE_Br                      0 /* Save Br       */
#define DO_SAVE_NO                      1 /* Save NO       */
#define DO_SAVE_NO3                     1 /* Save NO3      */
#define DO_SAVE_Cl                      0 /* Save Cl       */
#define DO_SAVE_O                       0 /* Save O        */
#define DO_SAVE_O1D                     0 /* Save O1D      */
#define DO_SAVE_O3                      1 /* Save O3       */
#define DO_SAVE_HO2                     1 /* Save HO2      */
#define DO_SAVE_NO2                     1 /* Save NO2      */
#define DO_SAVE_OH                      1 /* Save OH       */
#define DO_SAVE_HBr                     0 /* Save HBr      */
#define DO_SAVE_HBr_L                   0 /* Save HBr_L    */
#define DO_SAVE_HCl                     0 /* Save HCl      */
#define DO_SAVE_HCl_L                   0 /* Save HCl_L    */
#define DO_SAVE_CO                      1 /* Save CO       */
#define DO_SAVE_MO2                     0 /* Save MO2      */

#endif /* OUTPUT_H_INCLUDED */
