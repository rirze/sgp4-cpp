#ifndef _ASTMATH_H_
#define _ASTMATH_H_

/*   *****************************************************************

                           Module - AstMATH.H

  This file contains most of the math procedures and functions.

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

     ***************************************************************** */

#include <math.h>
#include <stdio.h>

/* General object definitions */
typedef char          Byte;
typedef double        Real; 
typedef int           Integer;
typedef long int      LINT;
typedef int           SINT;
typedef unsigned int  UINT;
typedef unsigned long ULINT;

extern char Show;
extern FILE *FileOut;

class Matrix;

/* error codes for vector accessing */
typedef enum {VOK, EMPTY_VECTOR, BAD_INDEX} VStatus;

/*
 * Class Definitions
 */

/*****************************************************************************
 *                    Define the Vector Class
*****************************************************************************/
class Vector
{
  private:
    #define VDIM 3
    Real vec[VDIM];
    UINT  dim;
    Real  mag;

  protected:
    void xmag(void);
  public:
    /* Constructor/Destructor */
    Vector(void);
    Vector(UINT d);
   ~Vector(void);

    /* Overloaded Operators */
    Real&   operator()(UINT);
    Vector& operator=(Vector&);
    Vector& operator+(Vector&);
    Vector& operator-(Vector&);
    Real    operator*(Vector&);
    Vector& operator*(Real);
    Vector& operator*(Matrix&);
    Vector& operator/(Real);

    /* User Functions */
    Vector& Add(Vector&);
    VStatus Add(Vector&, Vector&);
    Real    Angle(Vector&);
    VStatus Angle(Vector&, Real&);
    void    Clear(void);
    Vector& Cross(Vector&);
    VStatus Cross(Vector&, Vector&);
    UINT    Dim(void);
    void    Display(char*, UINT);
    Real    Dot(Vector&);
    Real    Get(UINT);
    VStatus Get(Real&, UINT);
    Real    Mag(void);
    void    Mag(Real);
    Vector& Norm(void);
    VStatus Norm(Vector&);
    Vector& Rot1(Real);
    VStatus Rot1(Real, Vector&);
    Vector& Rot2(Real);
    VStatus Rot2(Real, Vector&);
    Vector& Rot3(Real);
    VStatus Rot3(Real, Vector&);
    VStatus Set(Real, UINT);

    /* Friend Functions */
    friend Vector& operator*(Real, Vector&);
    friend VStatus Mult(Vector&, Matrix&);
};

/* error codes for matrix accessing */
typedef enum {MOK, EMPTY_MATRIX, BAD_INDEX1, BAD_INDEX2, NOT_SQUARE} MStatus;

/*****************************************************************************
 *                    Define the Matrix Class
*****************************************************************************/
class Matrix
{
  private:
    typedef Real *mvtype;
    #define MDIMR 3
    #define MDIMC 3
    Real  mat[MDIMR][MDIMC];
    UINT   rdim;
    UINT   cdim;
  public:
    Matrix(void);
    Matrix(UINT, UINT);
   ~Matrix(void);

    /* Overloaded Operators */
    Real&   operator()(UINT, UINT);
    Matrix& operator+(Real);
    Matrix& operator+(Matrix&);
    Matrix& operator-(Real);
    Matrix& operator-(Matrix&);
    Matrix& operator*(Real);
    Matrix& operator*(Matrix&);
    Vector& operator*(Vector&);
    Matrix& operator/(Real);

    /* User Functions */
    void Clear(void);
    Real Determinant(void);
    MStatus Determinant(Real&);
    UINT DimC(void);
    UINT DimR(void);
    MStatus Display(char*, Integer);
    Real Get(UINT, UINT);
    MStatus Get(Real&, UINT, UINT);
    MStatus LUBkSubstitute(Matrix&, Integer*);
    MStatus LUDeCompose(Matrix&, Integer*);
    Matrix& Identity(void);
    MStatus Identity(Matrix&);
    MStatus Inverse(Matrix&);
    Matrix& Inverse(void);
    Matrix& Scale(Real);
    MStatus Scale(Matrix&, Real);
    MStatus Set(Real, UINT, UINT);
    bool    Square(void);
    Matrix& Transpose(void);
};

/* error codes for complex accessing */
typedef enum {COK, EMPTY_COMPLEX} CStatus;

/*****************************************************************************
 *                    Define the Complex Class
*****************************************************************************/
class Complex
{
  private:
    Real val[2];
    Real mag;
  protected:
    void xmag(void);
  public:
    Complex(void);
    Complex(Real, Real);
   ~Complex(void);

    /* Overloaded Operators */
    Complex& operator=(Complex&);
    Complex& operator+(Complex&);
    Complex& operator+(Real);
    Complex& operator-(Complex&);
    Complex& operator-(Real);
    Complex& operator*(Complex&);
    Complex& operator*(Real);
    Complex& operator/(Complex&);
    Complex& operator/(Real);

    /* User Functions */
    Complex& Add(Complex&);
    CStatus  Add(Complex&, Complex&);
    Real     Angle(Complex&);
    CStatus  Angle(Complex&, Real&);
    void     Clear(void);
    void     Display(void);
    Complex& Div(Complex&);
    Complex& Div(Real);
    Complex& Dot(Complex&);
    Real     GetI(void);
    CStatus  GetI(Real&);
    Real     GetR(void);
    CStatus  GetR(Real&);
    Real     Mag(void);
    void     Mag(Real&);
    Complex& Mult(Complex&);
    Complex& Mult(Real);
    CStatus  Set(Real, Real);
    CStatus  SetI(Real);
    CStatus  SetR(Real);
    Real     Sqrt(void);
    CStatus  Sqrt(Real&);
    Complex& Sub(Complex&);
    CStatus  Sub(Complex&, Complex&);

    /* Friend Functions */
    friend Complex& operator+(Real, Complex&);
    friend Complex& operator-(Real, Complex&);
    friend Complex& operator*(Real, Complex&);
    friend Complex& Div(Real, Complex&);
    friend Complex& Mult(Real, Complex&);
};

/*****************************************************************************
 *                    Non Class Mathematics
*****************************************************************************/
Real Arccosh(Real);
Real Arcsinh(Real);
Real Arctanh(Real);
Real Atan2(Real, Real);
Real Binomial(ULINT, ULINT);
Real ComplexSqrt(Real);
Real Cot(Real);
void Cubic(Real, Real, Real, Real, Real&, Real&, Real&, Real&, Real&, Real&);
Real Csc(Real);
#define EVEN(a) (((abs(a)%2) == 0) ? true : false)
void Factor(Real *, UINT, Real **);
Real Factorial(ULINT);
Real Fraction(Real);
bool IsInt(Real);

#define LOG(a) (log10(a))
Real Max(Real, Real);
#define MaxLongInt 2^64-1
Real Min(Real, Real);
#define ODD(a) (((abs(a)%2) == 1) ? true : false)
Real Mod(Real, Real);
#define PI M_PI
void Plane(Real, Real, Real, Real, Real, Real, Real, Real, Real,
           Real&, Real&, Real&,  Real&);
void Polyfit(UINT, UINT, Matrix, Matrix, Real&, Real&);
Real Power(Real, Real);
void Quadratic(Real, Real, Real, Real&, Real&, Real&, Real&);
void Quartic
    (
      Real, Real, Real, Real, Real, 
      Real&, Real&, Real&, Real&, Real&, Real&, Real&, Real&
    );
Real Round(Real);
Real Sec(Real);
Real Sgn(Real);

/*****************************************************************************
 *                    General Utilities
*****************************************************************************/
void FilePrint(Vector&, char*, Integer, FILE*);
void Print(Vector&, char*, Integer);
void FilePrint(Matrix&, char*, Integer, FILE*);
void Print(Matrix&, char*, Integer);

#endif
