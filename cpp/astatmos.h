/*------------------------------------------------------------------------------
*
*                              ASTATMOS.FOR
*
*  This file contains the routines to generate the binary file for
*  a period of time of the atmosphere. It has the capability to input
*  the 3-day, 45-day, and yearly prediction values so a continuous file
*  is generated. The routines also include methods to find the F10.7 and
*  ap values at a given time. Finally, the integrated driver for the
*  Raytheon version of TRACE is included. This is a single interafec for
*  the Jacchia, and MSIS models.
*
*                          Companion code for
*             Fundamentals of Astrodyanmics and Applications
*                                   2002
*                            by David Vallado
*
*     (H)               email valladodl@worldnet.att.net
*     (W) 303-344-6037, email davallado@west.raytheon.com
*
*     *****************************************************************
*
*  Current :
*            26 Apr 02  David Vallado
*                         Fixes to predicted values
*  Changes :
*            12 Mar 02  David Vallado
*                         Fixes to calling parameters
*            20 Jan 02  David Vallado
*                         Original baseline
*
*     *****************************************************************
*
*  Links with astutil, astmath, asttime,
*             jach64, jach70, msiscom, msis86, msis90, msis00, DTM
*
*     The main program must also open the Atmosrec.rec file so that it is
*     available during the program.
*
*
*      SUBROUTINE UpDateAtmosRec
*
*      SUBROUTINE WriteOutAtmos
*
*      SUBROUTINE FindAtmosParam( MJD, MFME, Interp,
*     &                           F107, F107Ctr81,Ap, AvgAp, Kp, SumKp )
*
*      SUBROUTINE InterfaceAtmos( JDE, MFME, AtmosModel,r,rSun,hkm,
*     &                           erkm,pi,Interp,  rhoden )
*
*
*      AtmosRec = RECORD
*                   MJD            : REAL
*                   Year, Mon, Day : INTEGER
*                   Kp(8)          : INTEGER
*                   SumKp          : INTEGER
*                   Ap(8)          : INTEGER
*                   AvgAp          : REAL
*                   F107           : REAL
*                   F107Type       : CHAR
*                   F107Ctr81      : REAL
*                   F107Ctr81Type  : CHAR
*                 END
*
* ----------------------------------------------------------------------------*/
#include "astmath.h"

typedef struct atmosrec
{
  Real MJD;
  SINT Year, Mon, Day;
  SINT Kp[8];
  SINT SumKp;
  SINT Ap[8];
  Real AvgAp;
  Real F107;
  char F107Type;
  Real F107Ctr81;
  char F107Ctr81Type;
} AtmosRectype;
typedef AtmosRectype AtmosRec;

void FindAtmosParam(Real, Real, char, Real, Real, SINT*, Real, SINT*, SINT&);
void InterfaceAtmos
    (
      Real&, Real&, char*, Vector&, Vector&, Real&, Real&, char&, Real&, Real&
    );
void UpDateAtmosRec(void);
void WriteOutAtmos(void);
