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
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include "astatmos.h"
#include "asttime.h"
#include "astutil.h"

void FindAtmosParam
    (
      Real& MJD, Real& MFME, char& Interp, Real& F107,
      Real& F107Ctr81, SINT *Ap, Real& AvgAp, SINT *Kp, SINT& SumKp
    )
{
}

void InterfaceAtmos
    (
      Real& JDE, Real& MFME, char *AtmosModel, Vector& R, Vector& rSun,
      Real& hkm, Real& erkm, char& Interp, Real& RhoDen, Real& RhoDenp
    )
{
}

/*----------------------------------------------------------------------------
*
*                           SUBROUTINE UpdateAtmosRec
*
*  This subroutine updates the atmospheric parameters and creates the binary
*    datafile for use in OD. The subroutine also determines the centered F10.7
*    value so the program doesn't need to find this at runtime.
*    Note, the 81-day average is centered on the date of interest. If an
*    application needs the last 81-days, as with USSPACECOM, etc., one can
*    simply change the findatmos values to reference the record 40 days before.
*    These values will be the same. In this routine, also note that the 81-day
*    average is preset to the individual daily value during the first pass, it's
*    then updated in the final loop.
*
*  Author        : David Vallado                      303-344-6037   20 Jan 2002
*
*  Inputs          Description                          Range / Units
*    Year        - Year                                 1900 .. 2100
*    None        -
*    None        -
*    None        -
*    None        -
*    None        -
*    None        -
*    None        -
*    None        -
*
*  Outputs       :
*    None        -
*    None        -
*    None        -
*    None        -
*    None        -
*    None        -
*    None        -
*    None        -
*
*  Locals        :
*    None        -
*
*  Coupling      :
*    None        -
*
*  References    :
*    Vallado       2001,
*
*-----------------------------------------------------------------------------*/
void UpDateAtmosRec(void)
{
  Real  F107, JD, MJD,F107A[45], Cp, JDc, F107Ctr81, First81,
        F107Sum, JDc1, JDT[45], AvgAp, TempF107Ctr81;
  SINT  Year, Yr, Mon, Day, Kp[8], Bartels, NBartels, SumKp, Ap[8],Apa[45],
        C9, ISSN, FQual, RecNum, D[45], M[45], Y[45], GetIntMon, TemMon,
        RecNumf, ErrCode, TotRecNum, LMonth[12];
  UINT  i, j, k;
  char  FType, FCtrType, TempFCtrType;
  char  Dat[49], SunSpot[28], Month[3], MonStr[45], Date[8];
  char *AtmFileNameNOAA = "Atmosrec.rec";
  FILE *AtmFileNOAA;
  char *AtmFileNameNew = "atmosall.txt";
  FILE *AtmFileNew;
  char *AtmFileNameDSD = "dsd.txt";
  FILE *AtmFileDSD;
  char *AtmFileNameDGD = "dgd.txt";
  FILE *AtmFileDGD;
  char *AtmFileName45df = "45df.txt";
  FILE *AtmFile45df;
  AtmosRectype WriteRec;
  char  dummy[121];

  for (i = 0; i < 12; i++)
    LMonth[i] = 31;
  LMonth[ 2] = 28;
  LMonth[ 4] = 30;
  LMonth[ 6] = 30;
  LMonth[ 9] = 30;
  LMonth[11] = 30;

  if ((AtmFileNOAA = fopen(AtmFileNameNOAA, "w")) == NULL)
  {
    printf("Unable to open output record file %s\n", AtmFileNameNOAA);
    exit(0);
  }

  /* 
  ---- Read in actual data from ngdc.noaa.gov (within 1 month of current) -----
  */
  if ((AtmFileNew = fopen(AtmFileNameNew, "r")) == NULL)
  {
    printf("Unable to open update record set %s\n", AtmFileNameNew);
    exit(0);
  }

  RecNum = 0;
  while (!feof(AtmFileNew))
  {
    fscanf(AtmFileNew, "%d %d %d %d %d", &Yr, &Mon, &Day, &Bartels, &NBartels);
    for (i = 0; i < 8; i++)
      fscanf(AtmFileNew, "%d", Kp[i]);
    fscanf(AtmFileNew, "%d", &SumKp);
    for (i = 0; i < 8; i++)
      fscanf(AtmFileNew, "%d", Ap[i]);
    fscanf(AtmFileNew, "%f %f %d %d %f %d",
                       &AvgAp, &Cp, &C9, &ISSN, &F107, &FQual);
    
    F107Ctr81 = F107;

    if (Yr < 50)
      Year = Yr + 2000;
    else
      Year = Yr + 1900;

    JulianDay(Year, Mon, Day, 0, 0, 0.0, JD);
    MJD = JD - 2400000.5;

    RecNum++;

    // Write the record
    WriteRec.MJD = MJD;
    WriteRec.Year = Year;
    WriteRec.Mon = Mon;
    WriteRec.Day = Day;
    for (i = 0; i < 8; i++)
      WriteRec.Kp[i] = Kp[i];
    WriteRec.SumKp = SumKp;
    for (i = 0; i < 8; i++)
      WriteRec.Ap[i] = Ap[i];
    WriteRec.AvgAp = AvgAp;
    WriteRec.F107 = F107;
    WriteRec.F107Type = 'A';
    WriteRec.F107Ctr81 = F107Ctr81;
    WriteRec.F107Ctr81Type = 'P';
    write(fileno(AtmFileNOAA), (char *)&WriteRec, sizeof(AtmosRectype));
  }

  printf("Finished setup (curr    ) .rec to %14.4f%5d%3d%3d%5d\n",
         JD, Year, Mon, Day, RecNum);
  fclose(AtmFileNew);

  /*
  --------- Read in up to last 30 days of actual data from sel.noaa.gov -------
  This should be set up to search for the first date that's after
  the previous file read. Then, read from the DSD until the end,
  and simply read from the DGD until the DSD is done. The dgd will
  sometimes have some -1 fields at the end.
  */
  if ((AtmFileDSD = fopen(AtmFileNameDSD, "r")) == NULL)
  {
    printf("Unable to open DSD file %s\n", AtmFileNameDSD);
    exit(0);
  }
  if ((AtmFileDGD = fopen(AtmFileNameDGD, "r")) == NULL)
  {
    printf("Unable to open DGD file %s\n", AtmFileNameDGD);
    exit(0);
  }
  for (i = 0; i < 12; i++)
    fgets(dummy, 120, AtmFileDGD);

  // Set temporary value
  JDc = JD - 1.0;
  // Find the current date in the predicted files. Do separately
  // because they may not start on the same date
  i = 0;
  while ((JDc1 < JD) && (i < 29) && (!feof(AtmFileDGD)))
  {
    fscanf(AtmFileDSD, "%d %d %d %f\n", &Year, &Mon, &Day, &F107);
    JulianDay(Year, Mon, Day, 0, 0, 0.0, JDc);
    i++;
  }
  i = 0;
  while ((JDc1 < JD) && (i < 29) && (!feof(AtmFileDGD)))
  {
    fscanf(AtmFileDGD, "4d%3d%3d%49c%3d", &Year, &Mon, &Day, Dat, &Apa[0]);
    JulianDay(Year, Mon, Day, 0, 0, 0.0, JDc);
    i++;
  }

  // there may be a problem here if dsd and dgd don't start on the same date
  while ((!feof(AtmFileDSD)) && (!feof(AtmFileDGD)))
  {
    fscanf(AtmFileDSD, "%d %d %d %f", &Year, &Mon, &Day, &F107);
    fscanf(AtmFileDGD, "%4d%3d%3d%49c%3d", &Year, &Mon, &Day, Dat, &Apa[0]);

    for (UINT j = 0; j < 8; j++)
    {
      Ap[j] = Apa[0];
      Kp[j] = (SINT)(log(Ap[j]) * 1.7f - 1.6);
    }
    AvgAp  = Apa[0];
    SumKp  = Kp[7] * 8;
    RecNum++;
    printf("%d %d %d %d %7f\n", RecNum, Year, Mon, Day, Ap);
    JulianDay( Year, Mon, Day, 0, 0, 0.0, JD);
    MJD = JD - 2400000.5;
    F107Ctr81 = F107;

    // Write the record
    WriteRec.MJD = MJD;
    WriteRec.Year = Year;
    WriteRec.Mon = Mon;
    WriteRec.Day = Day;
    for (i = 0; i < 8; i++)
      WriteRec.Kp[i] = Kp[i];
    WriteRec.SumKp = SumKp;
    for (i = 0; i < 8; i++)
      WriteRec.Ap[i] = Ap[i];
    WriteRec.AvgAp = AvgAp;
    WriteRec.F107 = F107;
    WriteRec.F107Type = 'A';
    WriteRec.F107Ctr81 = F107Ctr81;
    WriteRec.F107Ctr81Type = 'P';
    write(fileno(AtmFileNOAA), (char *)&WriteRec, sizeof(AtmosRectype));
  }

  printf("Finished setup (last 30d) .rec to %14.4f%5d%3d%3d%5d\n",
         JD, Year, Mon, Day, RecNum);
  fclose(AtmFileDSD);
  fclose(AtmFileDGD);

  // ------- Read in next 3 days of predicted data from daypre.txt -------
  // ------- Read in up to next 45 days of predicted data from sel.noaa.gov ---
  if ((AtmFile45df = fopen(AtmFileName45df, "r")) == NULL)
  {
    printf("Unable to open DSD file %s\n", AtmFileName45df);
    exit(0);
  }
  for (i = 0; i < 10; i++)
    fgets(dummy, 120, AtmFile45df);

  // Read all 45 days of data in
  k = 1;
  for (j = 0; j < 9; j++)
  {
  }

}

void WriteOutAtmos(void)
{
}
