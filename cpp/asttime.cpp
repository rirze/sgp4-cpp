/*     ----------------------------------------------------------------      


                               UNIT ASTTIME;


    This file contains fundamental Astrodynamic procedures and functions     
    relating to the time functions. These routines are discussed in Ch       
    3 and Ch 5.                                                              

                            Companion code for                               
               Fundamentals of Astrodyanmics and Applications                
                                     2001                                    
                              by David Vallado                               

       (H)               email valladodl@worldnet.att.net                    
       (W) 303-344-6037, email davallado@west.raytheon.com                   

       *****************************************************************     

    Current :                                                                
              14 May 01  David Vallado                                       
                           2nd edition baseline                              
    Changes :                                                                
              23 Nov 87  David Vallado                                       
                           Original baseline                                 

       ----------------------------------------------------------------      

                                  IMPLEMENTATION

       ----------------------------------------------------------------      */
#include <math.h>
#include "asttime.h"

/*
 * Object Definitions
*/
extern char  Show;
extern FILE *FileOut;
Integer  LMonth[12];
char    *DayTitle[7];
char    *MonthTitle[7];

/*
 * Implementatio Code
*/
/*------------------------------------------------------------------------------
|
|                           PROCEDURE DAYLIGHTST
|
|  This PROCEDURE finds the dates for switiching to daylight savings time in
|    a given year. The date is set as the 1st Sunday in April and the last
|    Sunday in October. The DST dates are adjusted -10 to get 0200, -Zone to get
|    the local time zone, and -1.0 to process the stop because the local time is
|    on DST before a stop.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Year        - Year                           1900 .. 2100
|    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
|
|  Outputs       :
|    StartDay    - Day in April when DST begins   1 .. 28,29,30,31
|    StopDay     - Day in October when DST ends   1 .. 28,29,30,31
|
|  Locals        :
|    DW          - Day of the week                1 .. 7
|    JDStartDST  - Julian date of start           Days from 4713 BC
|    JDStopDST   - Julian date of stop            Days from 4713 BC
|    Zone        - Time zone of site. Default
|                  of 0.0 gives Greenwich         hrs
|
|  Coupling      :
|    JULIANDAY   - Find the Julian Date
|
|  References    :
|    Vallado       2001, 185-186
|
 -----------------------------------------------------------------------------*/
void DayLightSt
    (
      Integer Year, Real Lon, 
      Integer& StartDay, Integer& StopDay, Real& JDStartDST, Real& JDStopDST
    )
{
  const Real Rad2Deg = 190.0 / PI;
  Integer    DW;
  Real       Zone;

  /* ------- Find time zone information to adjust to a site ------- */
  Zone = Round(Lon * Rad2Deg / 15.0);
  if (Zone > 0.0)
    Zone = Zone - 24.0;

  StartDay = 0;
  while ((DW != 1) && (StartDay < 8))
  {
    StartDay = StartDay + 1;
    JulianDay(Year, 4, StartDay, 12, 0, 0.0, JDStartDST);
    DW = Integer(JDStartDST - 7 * Integer((JDStartDST + 1) / 7) + 2);
  }
  JDStartDST = JDStartDST + (-10.0 - Zone) / 24.0;  // set to 0200 local

  StopDay = 32;
  while ((DW != 1) && (StartDay < 24))
  {
    StopDay = StopDay - 1;
    JulianDay(Year, 10, StopDay, 12, 0, 0.0, JDStopDST);
    DW = Integer(JDStopDST - 7 * Integer((JDStopDST + 1) / 7) + 2);
  }
  JDStopDST = JDStopDST + (-10.0 - Zone - 1.0) / 24.0;  // set to 0200 local
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE DAYOFWEEK
|
|  This FUNCTION finds the day of the week. Integers are used for the days,
|    1 = 'SUN', 2 = 'Mon', ... 7 = 'Sat'.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    JD          - Julian date of interest        days from 4713 BC
|
|  OutPuts       :
|    DAYOFWEEK   - answer                         1 to 7
|
|  Locals        :
|    None.
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001, 185, Eq 3-39
|
 -----------------------------------------------------------------------------*/
Integer DayOfWeek(Real JD)
{
  /* ----- Be sure JD is at 0.0 h on the day of interest ----- */
  JD = Integer(JD + 0.5);
  return Integer(JD - 7 * Integer((JD + 1) / 7) + 2);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE DAYS2MDHMS
|
|  This PROCEDURE converts the day of the year, days, to the equivalent month
|    day, hour, Minute and second.
|
|  Algorithm     : Set up array for the Number of days per month
|                  Find Leap Year - Use 1900 because 2000 is a leap year
|                  Loop through a Temp value while the value is < the days
|                  Perform INTEGER conversions to the correct day and month
|                  Convert remainder into H M S using type conversions
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Year        - Year                           1900 .. 2100
|    Days        - Julian Day of the year         0.0  .. 366.0
|
|  OutPuts       :
|    Mon         - Month                          1 .. 12
|    Day         - Day                            1 .. 28,29,30,31
|    Hr          - Hour                           0 .. 23
|    MIN         - Minute                         0 .. 59
|    SEC         - Second                         0.0 .. 59.999
|
|  Locals        :
|    DayofYr     - Day of year
|    Temp        - Temporary EXTENDED values
|    IntTemp     - Temporary INTEGER value
|    i           - Index
|    LMonth[12]  - INTEGER Array containing the Number of days per month
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
void Days2MDHMS
    (
      Integer Year, Real Days,
      Integer& Mon, Integer& Day, Integer& Hr, Integer& Min, Real& Sec
    )
{
  Integer i, IntTemp, DayofYr, LeapDay;
  Real    Temp;
  
  DayofYr = Integer(Days);

  /* ----------------- Find month and Day of month ---------------- */
  if ( ((Year % 4) == 0) && (((Year % 400) != 0) || ((Year % 400) == 0)))
    LeapDay = 1;
  else
    LeapDay = 0;

  i = 1;
  IntTemp = 0;
  while ((DayofYr > IntTemp + LMonth[i-1] + LeapDay) && (i < 12))
  {
    IntTemp += LMonth[i-1] + ((i == 1) ? LeapDay : 0);
    i++;
  }
  Mon = i;
  Day = DayofYr - IntTemp;

  /* ----------------- Find hours minutes and seconds ------------- */
  Temp = (Days - DayofYr) * 24.0;
  Hr   = Integer(Temp);
  Temp = (Temp - Hr) * 60.0;
  Min  = Integer(Temp);
  Sec  = (Temp - Min) * 60.0;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE DMS_RAD
|
|  This PROCEDURE converts Degrees, minutes and seconds into radians.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Deg         - Degrees                        0 .. 360
|    MIN         - minutes                        0 .. 59
|    SEC         - Seconds                        0.0 .. 59.99
|    Direction   - Which set of vars to output    FROM  TOO
|
|  OutPuts       :
|    DMS         - Result                         rad
|
|  Locals        :
|    Temp        - Temporary variable
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001, 199, Alg 17 Alg 18, Ex 3-8
|
 -----------------------------------------------------------------------------*/
void DMS_Rad(Integer& Deg, Integer& Min, Real& Sec, Direction dir, Real& DMS)
{
  Real Temp;

  if (dir == FROM)
  {
    Temp = DMS * 180.0 * PI;
    Deg  = Integer(Temp);
    Min  = Integer((Temp = Deg) * 60.0);
    Sec  = (Temp - Deg - Min / 60.0) * 3600.0;
  }
  else
    DMS = (Deg + Min / 60.0 + Sec / 3600.0) * PI / 180.0;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE FINDDAYS
|
|  This PROCEDURE finds the fractional days through a year given the year,
|    month, day, hour, Minute and second.
|
|  Algorithm     : Set up array for the Number of days per month
|                  Find Leap Year - Use 1900 because 2000 is a leap year
|                  Check for a leap year
|                  Loop to find the elapsed days in the year
|
|  Author        : David Vallado                  303-344-6037
|
|  Inputs          Description                    Range / Units
|    Year        - Year                           1900 .. 2100
|    Mon         - Month                          1 .. 12
|    Day         - Day                            1 .. 28,29,30,31
|    Hr          - Hour                           0 .. 23
|    MIN         - Minute                         0 .. 59
|    SEC         - Second                         0.0 .. 59.999
|
|  OutPuts       :
|    Days        - Day of year plus fraction of a
|                    day                          days
|
|  Locals        :
|    LMonth      - Length of months of year
|    i           - Index
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001, 202-203, Ex 3-12
|
 -----------------------------------------------------------------------------*/
void FindDays
    (
      Integer Year, Integer Month, Integer Day, Integer Hr, Integer Min, 
      Real Sec, Real& Days
    )
{
  Integer i, LeapDay;

  if (((Year - 1900) % 4) == 0)
    LeapDay = 1;
  else
    LeapDay = 0;

  i    = 1;
  Days = 0.0;
  while ((i < Month) && (i < 12))
  {
    Days += LMonth[i-1] + (i == 1) ? LeapDay : 0;
    i++;
  }
  Days += Day + Hr / 24.0 + Min / 1440.0 + Sec / 86400.0;
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION GETINTMON / DAY
|
|  This FUNCTION finds the INTEGER equivalent of the 3 character string
|    representation of month and the day.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    MonStr      - Month name                     'Jan','Feb' ...
|    DayStr      - Day name string                'Sun','Mon' ...
|
|  OutPuts       :
|    GETINTMON   - INTEGER Month equivalent       1 .. 12
|    GETINTDAY   - INTEGER Day equivalent         1 .. 7
|
|  Locals        :
|    i           - Index
|
|  Coupling      :
|    UPCASESTR   - Converts a string to all uppercase characters -local proc
|
 -----------------------------------------------------------------------------*/
Byte GetIntDay(Str3 DayStr)
{
  Byte i;
  char DayTmp[8], TitleTmp[8];

  strcpy(DayTmp, DayStr);
  UpCaseStr(DayTmp);
  i = 0;
  strcpy(TitleTmp, DayTitle[i]);
  UpCaseStr(TitleTmp);
  while ((strcmp(TitleTmp, DayTmp) != 0) && (i < 7))
  {
    i++;
    strcpy(TitleTmp, DayTitle[i]);
    UpCaseStr(TitleTmp);
  }
  return i;
}

Byte GetIntMon(Str3 MonStr)
{
  Byte i;
  char MonTmp[8], TitleTmp[8];

  strcpy(MonTmp, MonStr);
  UpCaseStr(MonTmp);
  i = 0;
  strcpy(TitleTmp, MonthTitle[i]);
  UpCaseStr(TitleTmp);
  while ((strcmp(TitleTmp, MonTmp) != 0) && (i < 12))
  {
    i++;
    strcpy(TitleTmp, MonthTitle[i]);
    UpCaseStr(TitleTmp);
  }
  return i;
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION GSTIME
|
|  This FUNCTION finds the Greenwich SIDEREAL time.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    JDUT1       - Julian Date in UT1             days from 4713 BC
|
|  OutPuts       :
|    GSTIME      - Greenwich SIDEREAL Time        0 to 2Pi rad
|
|  Locals        :
|    Temp        - Temporary variable for reals   rad
|    TUT1        - Julian Centuries from the
|                  Jan 1, 2000 12 h epoch (UT1)
|
|  Coupling      :
|    REALMOD     - MOD FUNCTION for REAL variables
|
|  References    :
|    Vallado       2001, 191, Eq 3-45
|
 -----------------------------------------------------------------------------*/
Real GSTime(Real JDUT1)
{
  const Real TwoPi = 2.0 * PI;
  const Real Deg2Rad = PI / 180.0;
  Real       Temp, TUT1;

  TUT1 = (JDUT1 - 2451545.0) / 36525.0;
  Temp = -6.2E-6* TUT1 * TUT1 * TUT1 + 0.093104 * TUT1 * TUT1 + 
          (876600.0*3600 + 8640184.812866) * TUT1 + 67310.54841;  // sec
  Temp = Mod(Temp * Deg2Rad / 240.0, TwoPi); //360/86400 = 1/240, to deg, to rad

  /* ------------------------ Check quadrants --------------------- */
  if (Temp < 0.0)
    Temp += TwoPi;

  if (Show == 'Y')
    printf("%18.12f\n", TUT1);

  return Temp;
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION GSTIM0
|
|  This FUNCTION finds the Greenwich SIDEREAL time at the beginning of a year.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Year        - Year                           1998, 1999, etc.
|
|  OutPuts       :
|    GSTIM0      - Greenwich SIDEREAL Time        0 to 2Pi rad
|
|  Locals        :
|    JDUT1       - Julian Date in UT1             days from 4713 BC
|    Temp        - Temporary variable for Reals   rad
|    TUT1        - Julian Centuries from the
|                  Jan 1, 2000 12 h epoch (UT1)
|
|  Coupling      :
|    REALMOD       MOD FUNCTION for Real variables
|
|  References    :
|    Vallado       2001, 191, Eq 3-43
|
 -----------------------------------------------------------------------------*/
Real GSTim0(Integer Year)
{
  const Real TwoPi = 2.0 * PI;
  Real       JDUT1, Temp, TUT1;

  JDUT1 = 367.0 * Year - 
          (Integer((7 * (Year + Integer(10.0 / 12.0))) * 0.25) ) +
          (Integer(275.0 / 9.0)) + 1721014.5;
  TUT1 = (JDUT1 - 2451545.0) / 36525.0;
  Temp  = 1.75336855923327 + 628.331970688841 * TUT1 + 
          6.77071394490334E-06 * TUT1 * TUT1 -
          4.50876723431868E-10 * TUT1 * TUT1 * TUT1;

  /* ------------------------ Check quadrants --------------------- */
  Temp = Mod(Temp, TwoPi);
  if (Temp < 0.0)
    Temp += TwoPi;

  return Temp;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE HMS_RAD
|
|  This PROCEDURE converts Hours, minutes and seconds into radians.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Hr          - Hours                          0 .. 24
|    MIN         - minutes                        0 .. 59
|    SEC         - Seconds                        0.0 .. 59.99
|    Direction   - Which set of vars to output    FROM  TOO
|
|  OutPuts       :
|    HMS         - Result                         rad
|
|  Locals        :
|    Temp        - Conversion from hours to rad   0.261799
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001, 200, Alg 19 Alg 20, Ex 3-9
|
 -----------------------------------------------------------------------------*/
void HMS_Rad(Integer& Hr, Integer& Min, Real& Sec, Direction dir, Real& HMS)
{
  Real Temp;

  Temp = 15.0 * PI / 180.0;
  if (dir == FROM)
  {
    Temp = HMS / Temp;
    Hr   = Integer(Temp);
    Min  = Integer((Temp - Hr) * 60.0);
    Sec  = (Temp - Hr - Min / 60.0) * 3600.0;
  }
  else
    HMS = (Hr + Min / 60.0 + Sec / 3600.0) * Temp;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE HMS_SEC
|
|  This PROCEDURE converts Hours, minutes and Seconds into seconds from the
|    beginning of the day.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Hr          - Hours                          0 .. 24
|    MIN         - minutes                        0 .. 59
|    SEC         - Seconds                        0.0 .. 59.99
|    Direction   - Which set of vars to output    FROM  TOO
|
|  OutPuts       :
|    SEC         - Seconds                        0.0 .. 86400.0
|
|  Locals        :
|    Temp        - Temporary variable
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
void HMS_Sec(Integer& Hr, Integer& Min, Real& Sec, Direction dir, Real& UTSec)
{
  Real Temp;

  if (dir == FROM)
  {
    Temp = UTSec / 3600.0; // hr
    Hr   = Integer(Temp);
    Min  = Integer((Temp - Hr) * 60.0);
    Sec  = (Temp - Hr - Min / 60.0) * 3600.0;
  }
  else
    UTSec = Hr * 3600.0 + Min * 60.0 + Sec;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE HMS_UT
|
|  This PROCEDURE converts Hours, minutes and Seconds into Universal Time.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Hr          - Hours                          0 .. 24
|    MIN         - minutes                        0 .. 59
|    SEC         - Seconds                        0.0 .. 59.99
|    Direction   - Which set of vars to output    FROM  TOO
|
|  OutPuts       :
|    UT          - Universal Time                 HrMin.SEC
|
|  Locals        :
|    None.
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001, 201, Alg 21, Ex 3-10
|
 -----------------------------------------------------------------------------*/
void HMS_UT(Integer& Hr, Integer& Min, Real& Sec, Direction dir, Real& UT)
{
  if (dir == FROM)
  {
    Hr  = Integer(UT * 0.01);
    Min = Integer(UT - Hr * 100.0);
    Sec = Mod(UT, Integer(UT)) * 100.0;
  }
  else
    UT = Hr * 100.0 + Min + Sec * 0.01;
}

void InitTime(void)
{
  LMonth[ 0] = 31;
  LMonth[ 1] = 28;
  LMonth[ 2] = 31;
  LMonth[ 3] = 30;
  LMonth[ 4] = 31;
  LMonth[ 5] = 30;
  LMonth[ 6] = 31;
  LMonth[ 7] = 31;
  LMonth[ 8] = 30;
  LMonth[ 9] = 31;
  LMonth[10] = 30;
  LMonth[11] = 31;

  MonthTitle[ 0] = "Jan";
  MonthTitle[ 1] = "Feb";
  MonthTitle[ 2] = "Mar";
  MonthTitle[ 3] = "Apr";
  MonthTitle[ 4] = "May";
  MonthTitle[ 5] = "Jun";
  MonthTitle[ 6] = "Jul";
  MonthTitle[ 7] = "Aug";
  MonthTitle[ 8] = "Sep";
  MonthTitle[ 9] = "Oct";
  MonthTitle[10] = "Nov";
  MonthTitle[11] = "Dec";

  DayTitle[0] = "Sun";
  DayTitle[1] = "Mon";
  DayTitle[2] = "Tue";
  DayTitle[3] = "Wed";
  DayTitle[4] = "Thu";
  DayTitle[5] = "Fri";
  DayTitle[6] = "Sat";
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE INVJULIANDAY
|
|  This PROCEDURE finds the Year, month, day, hour, Minute and second
|  given the Julian date. TU can be UT1, TDT, TDB, etc.
|
|  Algorithm     : Set up starting values
|                  Find Leap Year - Use 1900 because 2000 is a leap year
|                  Find the elapsed days through the year in a loop
|                  Call routine to find each individual value
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    JD          - Julian Date                    days from 4713 BC
|
|  OutPuts       :
|    Year        - Year                           1900 .. 2100
|    Mon         - Month                          1 .. 12
|    Day         - Day                            1 .. 28,29,30,31
|    Hr          - Hour                           0 .. 23
|    MIN         - Minute                         0 .. 59
|    SEC         - Second                         0.0 .. 59.999
|
|  Locals        :
|    Days        - Day of year plus fractional
|                  portion of a day               days
|    Tu          - Julian Centuries from 0 h
|                  Jan 0, 1900
|    Temp        - Temporary real values
|    LeapYrs     - Number of Leap years from 1900
|
|  Coupling      :
|    DAYS2MDHMS  - Finds Month, day, hour, Minute and second given Days and Year
|
|  References    :
|    Vallado       2001, 203-205, Alg 22, Ex 3-13
|
 -----------------------------------------------------------------------------*/
void InvJulianDay
    (
      Real JD, 
      Integer& Year, Integer& Mon, Integer& Day, 
      Integer& Hr, Integer& Min, Real& Sec
    )
{
  Integer LeapYrs;
  Real    Days, Tu, Temp;

  /* --------------- Find Year and Days of the year --------------- */
  Temp    = JD - 2415019.5;
  Tu      = Temp / 365.25;
  Year    = 1900 + Integer(Tu);
  LeapYrs = Integer((Year - 1901) * 0.25);
  // Nudge by 8.64x10-7 Sec to get even outputs
  Days    = Temp - ((Year - 1900) * 365.0 + LeapYrs) + 0.00000000001; 

  /* ------------ Check for case of beginning of a year ----------- */
  if (Days < 1.0)
  {
    Year    = Year - 1;
    LeapYrs = Integer((Year - 1901) * 0.25);
    Days    = Temp - ((Year - 1900) * 365.0 + LeapYrs);
  }

  /* ----------------- Find remaing data  ------------------------- */
  Days2MDHMS(Year, Days, Mon, Day, Hr, Min, Sec);
  Sec = Sec - 0.00000086400;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE JULIANDAY
|
|  This PROCEDURE finds the Julian date given the Year, Month, Day, and Time.
|    The Julian date is defined by each elapsed day since noon, Jan 1, 4713 BC.
|
|  Algorithm     : Calculate the answer in one step for efficiency
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Year        - Year                           1900 .. 2100
|    Mon         - Month                          1 .. 12
|    Day         - Day                            1 .. 28,29,30,31
|    Hr          - Universal Time Hour            0 .. 23
|    MIN         - Universal Time MIN             0 .. 59
|    SEC         - Universal Time SEC             0.0 .. 59.999
|    WhichType   - Julian or Gregorian calender   'J' or 'G'
|
|  Outputs       :
|    JD          - Julian Date                    days from 4713 BC
|
|  Locals        :
|    B           - Var to aid Gregorian dates
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001, 186-188, Alg 14, Ex 3-14
|
 -----------------------------------------------------------------------------*/
void JulianDay
    (
      Integer Year, Integer Mon, Integer Day, Integer Hr, Integer Min, Real Sec,
      Real& JD
    )
{
  JD = 367.0 * Year -
       Integer((7 * (Year + Integer((Mon + 9) / 12))) * 0.25) +
       Integer( 275 * Mon / 9 ) +
       Day + 1721013.5 +
       ((Sec / 60.0 + Min) / 60.0 + Hr) / 24.0;  // UT in days 
       // - 0.5*SGN(100.0*Year + Mon - 190002.5) + 0.5;
}

void JulianDayAll
    (
      Integer Year, Integer Mon, Integer Day, Integer Hr, Integer Min, Real Sec,
      char WhichType, Real& JD
    )
{
  Real B;

  if (Mon <= 2)
  {
    Year = Year - 1;
    Mon  = Mon + 12;
  }
  if (WhichType == 'J')
    /* --------- Use for Julian calender, every 4 years --------- */
    B = 0.0;
  else
  {
    /* ---------------------- Use for Gregorian ----------------- */
    B  = 2 - Integer(Year * 0.01) + Integer(Integer(Year * 0.01) * 0.25);
    JD = Integer(365.25 * (Year + 4716)) +
         Integer(30.6001 * (Mon + 1)) +
         Day + B - 1524.5 +
         ((Sec / 60.0 + Min) / 60.0 + Hr) / 24.0;  // UT in days
  }
}

/*------------------------------------------------------------------------------|
|                           PROCEDURE LSTIME
|
|  This PROCEDURE finds the Local SIDEREAL time at a given location.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
|    JDUT1       - Julian Date in UT1             days from 4713 BC
|
|  OutPuts       :
|    LST         - Local SIDEREAL Time            0.0 to 2Pi rad
|    GST         - Greenwich SIDEREAL Time        0.0 to 2Pi rad
|
|  Locals        :
|    None.
|
|  Coupling      :
|    REALMOD       MOD FUNCTION for REAL variables
|    GSTIME        Finds the Greenwich SIDEREAL Time
|
|  References    :
|    Vallado       2001, 189-193. Alg 15, Ex 3-5
|
 -----------------------------------------------------------------------------*/
void LSTime(Real Lon, Real JDUT1, Real& LST, Real& GST)
{
  const Real TwoPi = 2.0 * PI;

  GST = GSTime(JDUT1);
  LST = Lon + GST;

  /* ------------------------ Check quadrants --------------------- */
  LST = Mod(LST, TwoPi);
  if (LST < 0.0)
    LST = LST + TwoPi;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE MOONRISESET
|
|  This PROCEDURE finds the Universal time for Moonrise and Moonset given the
|    day and SITE location.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    JD          - Julian Date                    days from 4713 BC
|    Latgd       - SITE latitude (SOUTH -)        -65ø to 65ø rad
|    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
|
|  OutPuts       :
|    UTMoonRise  - Universal time of Moonrise     hrs
|    UTMoonSet   - Universal time of Moonset      hrs
|    MoonPhaseAng- Phase angle of the Moon        deg
|    Error       - Error Parameter
|
|  Locals        :
|    MoonAngle   - ANGLE between the Moon vector
|                  and a point on the Earth       rad
|    JDTemp      - Julian date for Moonrise/set   days from 4713 BC
|    UTTemp      - Temporary UT time              days
|    TUT1        - Julian Centuries from the
|                  Jan 1, 2000 12 h epoch (UT1)
|    RtAsc       - Right ascension                rad
|    Decl        - Declination                    rad
|    MeanLonMoon -                                rad
|    MeanAnomaly -                                rad
|    EclpLong    - Longitude of the ecliptic      rad
|    Obliquity   - Obliquity of the ecliptic      rad
|    Try
|    l, m, n     - Direction cosines
|    EclpLat
|    MoonGHA, MoonGHAn
|    DGHA
|    LHAn
|    LST
|    DeltaUT
|    t, tn
|    LonEclSun
|    LonEclMoon
|    TTDB
|    x
|    GST         - for 0 h UTC of each day        rad
|    LHA         - Local hour ANGLE               rad
|    Year        - Year                           1900 .. 2100
|    Mon         - Month                          1 .. 12
|    Day         - Day                            1 .. 28,29,30,31
|    Hr          - Hour                           0 .. 23
|    MIN         - Minute                         0 .. 59
|    SEC         - Second                         0.0 .. 59.999
|    Opt         - Idx to do rise and set calc    1,2
|
|  Coupling      :
|    INVJULIANDAY- Finds the Year day mon hr MIN SEC from the Julian Date
|    JULIANDAY   - Finds the Julian date given Year, mon day, hr, MIN, SEC
|    ARCSIN      - Arc sine FUNCTION
|    ARCCOS      - Arc cosine FUNCTION
|    REALMOD     - MOD FUNCTION for REAL variables
|
|  References    :
|    Vallado       2001, 276-280, Alg 32, Ex 5-4
|
-----------------------------------------------------------------------------*/
void MoonRiseSet
    (
      Real JD, Real Latgd, Real Lon, 
      Real& UTMoonRise, Real& UTMoonSet, Real& MoonPhaseAng, char* Error
    )
{
  const Real TwoPi = 2.0 * PI;
  const Real Deg2Rad = PI / 180.0;

  Byte    Try;
  Real    l, m, n, EclpLong, EclpLat, Obliquity, MoonGHA, DGHA, LHAn, LHA, LST,
          MoonGHAn, DeltaUT, tn, GST, t, LonEclSun, LonEclMoon, MeanAnomaly,
          MeanLong, TTDB, x, Sec, JDTemp, UTTemp, RtAsc, Decl;
  Integer Opt, i, Year, Month, Day, Hr, Min;

  strcpy(Error, "ok");
  /* ----------- Do once for MoonRise (1), then set (2) ----------- */
  /* --------------- Make sure lon is within +- 180 deg ----------- */
  if (Lon > PI)
    Lon = Lon - TwoPi;
  if (Lon < -PI)
    Lon = Lon + TwoPi;

  Try = 1;
  Opt = 1;
  while(Opt <= 2)
  {
    InvJulianDay(JD, Year, Month, Day, Hr, Min, Sec);
    JulianDay(Year, Month, Day, 0, 0, 0.0, JDTemp);
    if (Show == 'Y')
      fprintf(FileOut, " JD at 0 hr %18.12f\n", JD);
    UTTemp = 0.5;   // days

    if (Try == 2)
      if (Opt == 1)
        UTTemp = 0.25;  // days
      else
        UTTemp = 0.75;  // days
        
    tn = UTTemp;        // Initial Guess
    JDTemp = JDTemp + UTTemp;
    i = 0;

    while (1 == 1)
    {
      TTDB = (JDTemp - 2451545.0) / 36525.0;
      if (Show == 'Y')
        fprintf(FileOut, "Optin %3d JDTemp %18.12f TTDB %14.10f\n",
                         Opt, JDTemp, TTDB);
      EclpLong = 218.32 + 481267.883 * TTDB +
                   6.29 * sin((134.9 + 477198.85 * TTDB) * Deg2Rad) -
                   1.27 * sin((259.2 - 413335.38 * TTDB) * Deg2Rad) +
                   0.66 * sin((235.7 + 890534.23 * TTDB) * Deg2Rad) +
                   0.21 * sin((269.9 + 954397.70 * TTDB) * Deg2Rad) -
                   0.19 * sin((357.5 +  35999.05 * TTDB) * Deg2Rad) -
                   0.11 * sin((186.6 + 966404.05 * TTDB) * Deg2Rad);   // Deg

      EclpLat =    5.13 * sin(( 93.3 + 483202.03 * TTDB) * Deg2Rad) +
                   0.28 * sin((228.2 + 960400.87 * TTDB) * Deg2Rad) -
                   0.28 * sin((318.3 +   6003.18 * TTDB) * Deg2Rad) -
                   0.17 * sin((217.6 - 407332.20 * TTDB) * Deg2Rad);   // Deg

      EclpLong  = Mod(EclpLong * Deg2Rad, TwoPi);
      EclpLat   = Mod(EclpLat  * Deg2Rad, TwoPi);
      Obliquity = 23.439291 - 0.0130042 * TTDB;  // deg
      Obliquity = Obliquity * Deg2Rad;
      if (Show == 'Y')
        fprintf(FileOut, "Ecl lat,lon,obl %11.7f %11.7f %11.7f\n",
                  EclpLat / Deg2Rad, EclpLong / Deg2Rad, Obliquity / Deg2Rad);

      /* ------------- Find the geocentric direction cosines ---------- */
      l = cos(EclpLat) * cos(EclpLong);
      m = cos(Obliquity) * cos(EclpLat) * sin(EclpLong) -
          sin(Obliquity) * sin(EclpLat);
      n = sin(Obliquity) * cos(EclpLat) * sin(EclpLong) +
          cos(Obliquity) * sin(EclpLat);

      RtAsc = Atan2(m, l);  // more accurate than COS() TAN()
      if (Show == 'Y')
        fprintf(FileOut, "dir cosines lmn %11.7f %11.7f %11.7f rtasc %11.7f\n",
                         l, m, n, RtAsc / Deg2Rad);

      /* ---- Check that RtAsc is in the same quadrant as EclpLong ---- */
      if (EclpLong < 0.0)
        EclpLong = EclpLong + TwoPi;   // make sure it's in 0 to 2Pi range

      if (fabs(EclpLong - RtAsc) > 0.5 * PI)
        RtAsc = RtAsc + 0.5 * PI * 
                Integer((EclpLong - RtAsc + 0.5) / (0.5 * PI));

      Decl = asin(n);
      if (Show == 'Y')
        fprintf(FileOut, "rtasc decl %11.7f %11.7f\n", 
                          RtAsc / Deg2Rad, Decl / Deg2Rad);

      LSTime(Lon, JDTemp, LST, GST);
      MoonGHAn = LST - Lon - RtAsc;
      if (Show == 'Y')
        fprintf(FileOut, "LST GHAMn %11.7f %11.7f\n", 
                          LST / Deg2Rad, MoonGHAn / Deg2Rad);

      if (i == 0)
      {
        LHA  = MoonGHAn + Lon;
        DGHA = 347.8 * Deg2Rad;
      }
      else
        DGHA = (MoonGHAn - MoonGHA) / DeltaUT;
      if (Show == 'Y')
        fprintf(FileOut, "%3d  LHA DGHA %11.7f %11.7f\n",
                         LHA / Deg2Rad, DGHA / Deg2Rad);

      if (DGHA < 0.0)
        DGHA = DGHA + TwoPi / fabs(DeltaUT);
      if (Show == 'Y')
        fprintf(FileOut, " DGHA  %11.7f same as before\n", DGHA / Deg2Rad);

      LHAn = 0.00233 - (sin(Latgd) * sin(Decl)) / (cos(Latgd) * cos(Decl));
      if (LHAn > 1.0)
        LHAn = 0.0;
      if (LHAn < -1.0)
        LHAn = -1.0;
      LHAn = acos(LHAn);
      if (Show == 'Y')
        fprintf(FileOut, " LHAn  %11.7f  opt %3d\n", LHAn / Deg2Rad, Opt);
      if (Opt == 1)
        LHAn = TwoPi - LHAn;

      if (fabs(DGHA) > 0.0001)
        DeltaUT = (LHAn - LHA) / DGHA;
      else
      {
        DeltaUT = (LHAn - LHA);
        DeltaUT = 1.0;  // day 
        fprintf(FileOut, "x\n");
      }
      if (Show == 'Y')
        fprintf(FileOut, " deltaut  %11.7f\n", DeltaUT);
      t = tn;
      if (fabs(DeltaUT)  > 0.5)
      {
        if (fabs(DGHA) > 0.001)
          if (DeltaUT < 0.0)
          {
            DeltaUT = DeltaUT + TwoPi / DGHA;
            if (fabs(DeltaUT) > 0.51)
              i = 6;
          }
          else
          {
            // day after
            DeltaUT = DeltaUT - TwoPi / DGHA;
            if (fabs(DeltaUT) > 0.51)
              i = 6;
          }
        else
        {
          DeltaUT = DeltaUT;  // Take another try
          fprintf(FileOut, "y");
        }
      }
      if (Show == 'Y')
        fprintf(FileOut, " deltaut  %11.7f\n", DeltaUT);
      tn     = UTTemp + DeltaUT;
      JDTemp = JDTemp - UTTemp + tn;
      i++;
      MoonGHA = MoonGHAn;
      if (Show == 'Y')
        fprintf(FileOut, "end %3d  tn t %11.7f %11.7f %18.8f\n",
                         tn, t, JDTemp);
      if ((fabs(tn - t)  < 0.008) || (i > 5))
        break;
    }

    UTTemp = tn * 24.0;
    if (i > 5)
      UTTemp = 9999.99;
    if (UTTemp < 9999.0)
      UTTemp = Mod( UTTemp, 24.0);
    if (UTTemp < 0.0)
      UTTemp = UTTemp + 24.0;
    if (UTTemp > 900)
      UTTemp = 24.0;

    if (Show == 'Y')
      fprintf(FileOut, "opt %3d  uttemp %11.7f  try %3d  ------- \n",
                       Opt, UTTemp, Try);

    switch (Opt)
    {
      case 1:
        UTMoonRise = UTTemp;
        break;
      case 2:
        UTMoonSet  = UTTemp;
    }

    Try++;
    if ((i > 5) && (Try < 3))   // only make 2 guesses
      printf("try #2 for %2d", Opt);
    else
    {
      // go on to next calc
      if ((i > 5) && (Try > 2))
      {
        if (Opt == 1)
          strcpy(Error, "No Rise");
        if (Opt == 2)
          strcpy(Error, "No Set");
      }
      Opt++;
      Try = 1;
    }
  }

  /* ------------- determine phase ANGLE of the MOON -------------- */
  MeanLong = 280.4606184 + 36000.77005361 * TTDB;
  MeanLong = Mod(MeanLong, 360.0);  // deg

  MeanAnomaly = 357.5277233 + 35999.05034 * TTDB;
  MeanAnomaly = Mod(MeanAnomaly * Deg2Rad, TwoPi);  // rad
  if (MeanAnomaly < 0.0)
    MeanAnomaly = TwoPi + MeanAnomaly;

  LonEclSun = MeanLong + 1.914666471 * sin(MeanAnomaly) +
              0.019994643 * sin(2.0 * MeanAnomaly);  // deg

  LonEclMoon = 218.32 + 481267.883 * TTDB +
                 6.29 * sin((134.9+ 477198.85 * TTDB) * Deg2Rad) -
                 1.27 * sin((259.2- 413335.38 * TTDB) * Deg2Rad) +
                 0.66 * sin((235.7+ 890534.23 * TTDB) * Deg2Rad) +
                 0.21 * sin((269.9+ 954397.70 * TTDB) * Deg2Rad) -
                 0.19 * sin((357.5+  35999.05 * TTDB) * Deg2Rad) -
                 0.11 * sin((186.6+ 966404.05 * TTDB) * Deg2Rad);    // Deg
  LonEclMoon = Mod(LonEclMoon, 360.0);

  MoonPhaseAng = LonEclMoon - LonEclSun;

  if (MoonPhaseAng < 0.0)
    MoonPhaseAng = 360.0 + MoonPhaseAng;   // deg
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE SUNRISESET
|
|  This PROCEDURE finds the Universal time for Sunrise and Sunset given the
|    day and SITE location.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    JD          - Julian Date                    days from 4713 BC
|    Latgd       - SITE latitude (SOUTH -)        -65ø to 65ø rad
|    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
|    WhichKind   - Character for which rise/set   'S' 'C' 'N' 'A'
|
|  OutPuts       :
|    UTSunRise   - Universal time of sunrise      hrs
|    UTSunSet    - Universal time of sunset       hrs
|    Error       - Error Parameter
|
|  Locals        :
|    SunAngle    - ANGLE between the SUN vector
|                  and a point on the Earth       rad
|    JDTemp      - Julian date for sunrise/set    days from 4713 BC
|    UTTemp      - Temporary UT time              days
|    TUT1        - Julian Centuries from the
|                  Jan 1, 2000 12 h epoch (UT1)
|    Ra          - Right ascension                rad
|    Decl        - Declination                    rad
|    MeanLonSun  -                                rad
|    MeanAnomalySun                               rad
|    LonEcliptic - Longitude of the ecliptic      rad
|    Obliquity   - Obliquity of the ecliptic      rad
|    GST         - for 0 h UTC of each day        rad
|    LHA         - Local hour ANGLE               rad
|    Year        - Year                           1900 .. 2100
|    Mon         - Month                          1 .. 12
|    Day         - Day                            1 .. 28,29,30,31
|    Hr          - Hour                           0 .. 23
|    MIN         - Minute                         0 .. 59
|    SEC         - Second                         0.0 .. 59.999
|    Opt         - Idx to do rise and set calc    1,2
|
|  Coupling      :
|    INVJULIANDAY- Finds the Year day mon hr MIN SEC from the Julian Date
|    JULIANDAY   - Finds the Julian date given Year, mon day, hr, MIN, SEC
|    ARCSIN      - Arc sine FUNCTION
|    ARCCOS      - Arc cosine FUNCTION
|    REALMOD     - MOD FUNCTION for REAL variables
|
|  References    :
|    Vallado       2001, 267-271, Alg 30, Ex 5-2
|
 -----------------------------------------------------------------------------*/
void SunRiseSet
    (
      Real JD, Real Latgd, Real Lon, char WhichKind, 
      Real& UTSunRise, Real& UTSunSet, char* Error
    )
{
  const Real TwoPi   = 2.0 * PI;
  const Real Rad2Deg = 180.0 / PI;
  const Real Deg2Rad = PI / 180.0;

  Integer Opt, Year, Month, Day, Hr, Min;
  Real JDTemp, UTTemp, SunAngle, TUT1, Ra, Sec, MeanLonSun, MeanAnomalySun,
       LonEcliptic, Decl, Obliquity, GST, LHA;

  strcpy(Error, "ok");
  /* --------------- Make sure lon is within +- 180 deg ----------- */
  if (Lon > PI)
    Lon = Lon - 2.0 * PI;
  if (Lon < -PI)
    Lon = Lon + 2.0 * PI;
  
  switch (WhichKind)
  {
    case 'S':
      SunAngle = (90.0 + 50.0 / 60.0) * Deg2Rad; // Sunrise / set
      break;
    case 'C':
      SunAngle = 96.0 * Deg2Rad;                 // Civil
      break;
    case 'N':
      SunAngle = 102.0 * Deg2Rad;                // Nautical
      break;
    case 'A':
      SunAngle = 108.0 * Deg2Rad;                // Astronomical
      break;
    default:
      SunAngle = (90.0 + 50.0 / 60.0) * Deg2Rad; // Sunrise / set
      break;
  }
  InvJulianDay(JD, Year, Month, Day, Hr, Min, Sec);

  for (Opt = 1; Opt <= 2; Opt++)
  {
    if (Opt == 1)
      JulianDay(Year, Month, Day, 6, 0, 0.0, JDTemp);   // Sunrise
    else
      JulianDay(Year, Month, Day, 18, 0, 0.0, JDTemp);  // Sunset
    JDTemp = JDTemp - Lon * Rad2Deg / 15.0 / 24.0;

    TUT1 = (JDTemp - 2451545.0) / 36525.0;
    MeanLonSun = 280.4606184 + 36000.77005361 * TUT1;   // deg

    MeanAnomalySun = 357.5277233 + 35999.05034 * TUT1;  // deg
    MeanAnomalySun = Mod(MeanAnomalySun * Deg2Rad, TwoPi);
    if (MeanAnomalySun < 0.0)
      MeanAnomalySun = MeanAnomalySun + TwoPi;

    LonEcliptic = MeanLonSun + 1.914666471 * sin(MeanAnomalySun) +
                  0.019994643 * sin(2.0 * MeanAnomalySun); // deg
    LonEcliptic = Mod(LonEcliptic * Deg2Rad,TwoPi);
    if (LonEcliptic < 0.0)
      LonEcliptic = LonEcliptic + TwoPi;

    Obliquity = 23.439291 - 0.0130042 * TUT1; // deg
    Obliquity = Obliquity * Deg2Rad;

    Ra   = atan(cos(Obliquity) * tan(LonEcliptic));
    Decl = asin(sin(Obliquity) * sin(LonEcliptic));

    if (Ra < 0.0)
      Ra = Ra + TwoPi;
    if (LonEcliptic > PI)
      Ra = Ra + PI;
    if ((LonEcliptic < PI) && (Ra > PI))
      Ra = Ra - PI;

    LHA = (cos(SunAngle) - sin(Decl) * sin(Latgd)) / 
          (cos(Decl) * cos(Latgd)); // rad

    if (fabs(LHA) <= 1.0)
      LHA = acos(LHA);      // rad
    else
      strcpy(Error, "Not ok");
    if (strcmp(Error, "ok") == 0)
    {
      if (Opt == 1)
        LHA = TwoPi - LHA;  // Sunrise only
      GST = 1.75336855923327 + 628.331970688841 * TUT1 +
            6.77071394490334E-06 * TUT1 * TUT1 -
            4.50876723431868E-10 * TUT1 * TUT1 * TUT1;
      GST = Mod(GST, TwoPi);
      if (GST < 0.0)
        GST = GST + TwoPi;
      UTTemp = LHA + Ra  - GST;
      UTTemp = UTTemp * Rad2Deg / 15.0;  // hrs
      UTTemp = Mod(UTTemp, 24.0);
      UTTemp = UTTemp - Lon * Rad2Deg / 15.0; // hrs
      if (UTTemp < 0.0)
      {
        UTTemp = UTTemp + 24.0;
        strcpy(Error, "Day before");
      }
      if (UTTemp > 24.0)
      {
        UTTemp = UTTemp - 24.0;
        strcpy(Error, "Day after");
      }
    }
    else
      UTTemp = 99.99;

    if (Opt == 1)
      UTSunRise = UTTemp;
    else
      UTSunSet  = UTTemp;
  }
}
