#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#include "astmath.h"

char  Show = 'N';
FILE *FileOut = NULL;

/*****************************************************************************
 * Class Method implementations of Vector
*****************************************************************************/

Vector::Vector(void)
{
  dim = VDIM; 
  mag = 0.0;
}

Vector::Vector(UINT d)
{
  dim = d;
  xmag();
}

Vector::~Vector(void)
{
}

void Vector::xmag(void)
{
  mag = 0.0;
  for (UINT i = 0; i < dim; i++)
    mag += vec[i] * vec[i];
  mag = sqrt(mag);
}

/*
 * Overloaded Operators
 */
Real& Vector::operator()(UINT d)
{
  if ((d <= dim) && (d > 0))
    return vec[d-1];
  else
    return vec[0];
}

Vector& Vector::operator=(Vector& v)
{
  if (this == &v)
    return *this;
  for (UINT i = 0; i < dim; i++)
    vec[i] = v.vec[i];
  mag = v.mag;
  return *this;
}

Vector& Vector::operator+(Vector& v)
{
  Vector *res;

  res = new Vector(v.dim);
  res->Clear();
  if (dim == v.dim)
    for (UINT i = 0; i < dim; i++)
      res->vec[i] = vec[i] + v.vec[i];
  res->mag = res->Mag();
  return *res;
}

Vector& Vector::operator-(Vector& v)
{
  Vector *res;

  res = new Vector(v.dim);
  res->Clear();
  if (dim == v.dim)
    for (UINT i = 0; i < dim; i++)
      res->vec[i] = vec[i] - v.vec[i];
  res->mag = res->Mag();
  return *res;
}

Real Vector::operator*(Vector& v)
{
  Real res = 0.0;

  if (dim == v.dim)
  {
    for (UINT i = 0; i < dim; i++)
      res += vec[i] * v.vec[i];
    return res;
  }
  else
    return 0.0;
}

Vector& Vector::operator*(Real v)
{
  Vector *res;

  res = new Vector(dim);
  for (UINT i = 0; i < dim; i++)
    res->vec[i] = v * vec[i];
  res->Mag();

  return *res;
}

Vector& Vector::operator*(Matrix& m)
{
  Vector *res;
  Real    t;

  res = new Vector(dim);
  for (UINT i = 0; i < dim; i++)
  {
    t = 0.0;
    for (UINT j = 0; j < dim; j++)
      t = t + vec[j] * m.Get(j+1, i+1);
    res->vec[i] = t;
  }

  return *res;
}

Vector& Vector::operator/(Real v)
{
  Vector *res;

  res = new Vector(dim);
  for (UINT i = 0; i < dim; i++)
    res->vec[i] = vec[i] / v;
  res->Mag();

  return *res;
}
/*
 * User Functions
 */
Vector& Vector::Add(Vector& v)
{
  Vector *res;

  res = new Vector(dim);
  for (UINT i = 0; i < dim; i++)
    res->vec[i] = vec[i] + v.vec[i];
  res->Mag();

  return *res;
}

VStatus Vector::Add(Vector& v, Vector& OutVec)
{
  for (UINT i = 0; i < dim; i++)
    OutVec.vec[i] = vec[i] + v.vec[i];
  OutVec.Mag();

  return VOK;
}

Real Vector::Angle(Vector& v)
{
  const Real Small = 0.000001;
  const Real Undefined = 999999.1;
  Real       t;

  if (mag * v.mag > Small * Small)
  {
    t = Dot(v) / (mag * v.mag);
    if (fabs(t) > 1.0)
      t = Sgn(t);
    return (acos(t));
  }
  else
    return Undefined;
}

VStatus Vector::Angle(Vector& v, Real& theta)
{
  const Real Small = 0.000001;
  const Real Undefined = 999999.1;
  Real       t;

  if (mag * v.mag > Small * Small)
  {
    t = Dot(v) / (mag * v.mag);
    if (fabs(t) > 1.0)
      t = Sgn(t);
    theta = acos(t);
  }
  else
    theta = Undefined;
}
void Vector::Clear(void)
{
  if (dim != 0)
    for (int i = 0; i < dim; i++)
      vec[i] = 0.0;
  mag = 0.0;
}

Vector& Vector::Cross(Vector& v)
{
  Vector *res;

  res = new Vector(v.dim);
  res->vec[0] = vec[1] * v.vec[2] - vec[2] * v.vec[1];
  res->vec[1] = vec[2] * v.vec[0] - vec[0] * v.vec[2];
  res->vec[2] = vec[0] * v.vec[1] - vec[1] * v.vec[0];
  res->mag = res->Mag();

  return *res;
}

VStatus Vector::Cross(Vector& v, Vector& OutVec)
{
  OutVec.vec[0] = vec[1] * v.vec[2] - vec[2] * v.vec[1];
  OutVec.vec[1] = vec[2] * v.vec[0] - vec[0] * v.vec[2];
  OutVec.vec[2] = vec[0] * v.vec[1] - vec[1] * v.vec[0];
  OutVec.Mag();

  return VOK;
}

UINT Vector::Dim(void)
{
  return dim;
}

void Vector::Display(char *s, UINT p)
{
  printf("%s\n", s);
  for (UINT i = 0; i < dim; i++)
    printf("v(%d) = %0.*f\n", i, p, vec[i]);
  printf("\n");
}

Real Vector::Dot(Vector& v)
{
  Real res = 0.0;

  if (dim == v.dim)
  {
    for (UINT i = 0; i < dim; i++)
      res += vec[i] * v.vec[i];
    return res;
  }
  else
    return 0.0;
}

Real Vector::Get(UINT d)
{
  if ((d <= dim) && (d > 0))
    return vec[d-1];
  else
    return 0.0;
}

VStatus Vector::Get(Real &val, UINT d)
{
  if ((d <= dim) && (d > 0))
  {
    val = vec[d-1];
    return VOK;
  }
  else
    return BAD_INDEX;
}

Real Vector::Mag(void)
{
  xmag();
  return mag;
}

void Vector::Mag(Real n)
{
  mag = n;
}

Vector& Vector::Norm(void)
{
  Vector *res;

  res = new Vector(dim);
  res->Clear();
  if (mag > 0.000001)
  {
    for (UINT i = 0; i < dim; i++)
      res->vec[i] = vec[i] / mag;
    res->mag = 1.0;
  }
  else
  {
    for (UINT i = 0; i < dim; i++)
      res->vec[i] = 0.0;
    res->mag = 0.0;
  }

  return *res;
}

VStatus Vector::Norm(Vector& OutVec)
{
  OutVec.Clear();
  if (mag > 0.000001)
  {
    for (UINT i = 0; i < dim; i++)
      OutVec.vec[i] = vec[i] / mag;
    OutVec.Mag();
  }
  else
  {
    for (UINT i = 0; i < dim; i++)
      OutVec.vec[i] = 0.0;
    OutVec.mag = 0.0;
  }

  return VOK;
}

Vector& Vector::Rot1(Real v)
{
  Real c, s;
  Vector *res;

  c = cos(v);
  s = sin(v);
  res = new Vector(dim);
  res->Clear();
  res->vec[2] = c * vec[2] - s * vec[1];
  res->vec[1] = c * vec[1] + s * vec[2];
  res->vec[0] = vec[0];
  res->Mag();

  return *res;
}

VStatus Vector::Rot1(Real v, Vector& OutVec)
{
  Real c, s;

  c = cos(v);
  s = sin(v);
  OutVec.Clear();
  OutVec.vec[2] = c * vec[2] - s * vec[1];
  OutVec.vec[1] = c * vec[1] + s * vec[2];
  OutVec.vec[0] = vec[0];
  OutVec.Mag();

  return VOK;
}

Vector& Vector::Rot2(Real v)
{
  Real c, s;
  Vector *res;

  c = cos(v);
  s = sin(v);
  res = new Vector(dim);
  res->Clear();
  res->vec[2] = c * vec[2] - s * vec[0];
  res->vec[0] = c * vec[0] + s * vec[2];
  res->vec[1] = vec[0];
  res->Mag();

  return *res;
}

VStatus Vector::Rot2(Real v, Vector& OutVec)
{
  Real c, s;

  c = cos(v);
  s = sin(v);
  OutVec.Clear();
  OutVec.vec[2] = c * vec[2] - s * vec[0];
  OutVec.vec[0] = c * vec[0] + s * vec[2];
  OutVec.vec[1] = vec[0];
  OutVec.Mag();

  return VOK;
}

Vector& Vector::Rot3(Real v)
{
  Real c, s;
  Vector *res;

  c = cos(v);
  s = sin(v);
  res = new Vector(dim);
  res->Clear();
  res->vec[1] = c * vec[1] - s * vec[0];
  res->vec[0] = c * vec[0] + s * vec[2];
  res->vec[2] = vec[0];
  res->Mag();

  return *res;
}

VStatus Vector::Rot3(Real v, Vector& OutVec)
{
  Real c, s;

  c = cos(v);
  s = sin(v);
  OutVec.Clear();
  OutVec.vec[1] = c * vec[1] - s * vec[0];
  OutVec.vec[0] = c * vec[0] + s * vec[2];
  OutVec.vec[2] = vec[0];
  OutVec.Mag();

  return VOK;
}

VStatus Vector::Set(Real val, UINT d)
{
  if (vec == NULL)
    return EMPTY_VECTOR;
  else if ((d > dim) || (d < 1))
    return BAD_INDEX;
  else
  {
    vec[d-1] = val;
    xmag();
    return VOK;
  }
}

/*
 * Friend Functions
 */
Vector& operator*(Real n, Vector& v)
{
  Vector *res = new Vector(v.dim);

  for (UINT i = 0; i < v.dim; i++)
    res->vec[i] = n * v.vec[i];

  return *res;
}

VStatus Mult(Vector& v, Matrix& m)
{
  Vector nvec(v.Dim());
  Real   t;
  UINT   d = v.Dim();;

  nvec.Clear();
  for(UINT i = 1; i <= d; i++)
  {
    t = 0.0;
    for(UINT j = 1; j <= d; j++)
      t += v.Get(j) * m.Get(j, i);
    nvec.Set(t, i);
  }
  for(UINT i = 1; i <= nvec.Dim(); i++)
    v.Set(nvec.Get(i), i);

  return VOK;
}

/*****************************************************************************
 * Class Method implementations of Matrix
*****************************************************************************/
Matrix::Matrix(void)
{
  rdim = MDIMR;
  cdim = MDIMC;
}

Matrix::Matrix(UINT r, UINT c)
{
  rdim = r;
  cdim = c;
}

Matrix::~Matrix(void)
{
}

/* Overloaded Operators */
Real& Matrix::operator()(UINT r, UINT c)
{
  if ((r >= 0) && (r <=rdim))
    if ((c >= 0) && (c <=cdim))
    {
      return mat[r][c];
    }
    else
      return mat[0][0];
  else
    return mat[0][0];
}

Matrix& Matrix::operator+(Real m)
{
  UINT offset;

  Matrix *nmat = new Matrix(rdim, cdim);

  for (UINT r = 0; r < rdim; r++)
    for (UINT c = 0; c < cdim; c++)
      nmat->mat[r][c] = mat[r][c] + m;
  return *nmat;
}

Matrix& Matrix::operator+(Matrix& m)
{
  UINT offset;

  Matrix *nmat = new Matrix(rdim, cdim);

  for (UINT r = 0; r < rdim; r++)
    for (UINT c = 0; c < cdim; c++)
      nmat->mat[r][c] = mat[r][c] + m.mat[r][c];
  return *nmat;
}

Matrix& Matrix::operator-(Real m)
{
  UINT offset;

  Matrix *nmat = new Matrix(rdim, cdim);

  for (UINT r = 0; r < rdim; r++)
    for (UINT c = 0; c < cdim; c++)
      nmat->mat[r][c] = mat[r][c] - m;
  return *nmat;
}

Matrix& Matrix::operator-(Matrix& m)
{
  Matrix *nmat = new Matrix(rdim, cdim);

  for (UINT i = 0; i < rdim; i++)
    for (UINT j = 0; j < cdim; j++)
      nmat->mat[i][j] = mat[i][j] - m.mat[i][j];
  return *nmat;
}

Matrix& Matrix::operator*(Real m)
{
  UINT offset;

  Matrix *nmat = new Matrix(rdim, cdim);

  for (UINT r = 0; r < rdim; r++)
    for (UINT c = 0; c < cdim; c++)
      nmat->mat[r][c] = mat[r][c] * m;
  return *nmat;
}

Matrix& Matrix::operator*(Matrix& m)
{
  Matrix *nmat = new Matrix(rdim, m.cdim);

  nmat->Clear();
  for (UINT i = 0; i < rdim; i++)
    for (UINT j = 0; j < m.cdim; j++)
      for (UINT k = 0; k < rdim; k++)
        nmat->mat[i][j] += mat[i][k] * m.mat[k][j];
  return *nmat;
}

Vector& Matrix::operator*(Vector& v)
{
  Vector *nvec;
  Real    t;

  nvec = new Vector(rdim);
  for (UINT i = 1; i <= rdim; i++)
  {
    t = 0.0;
    for (UINT j = 1; j <= cdim; j++)
      t = t + v.Get(j) * mat[i-1][j-1];
    nvec->Set(t, i);
  }

  return *nvec;
}

Matrix& Matrix::operator/(Real m)
{
  UINT offset;

  Matrix *nmat = new Matrix(rdim, cdim);

  for (UINT r = 0; r < rdim; r++)
    for (UINT c = 0; c < cdim; c++)
      nmat->mat[r][c] = mat[r][c] / m;
  return *nmat;
}

/*
 * User Functions
 */
void Matrix::Clear(void)
{
  if ((rdim != 0) && (cdim != 0))
    for (int r = 0; r < rdim; r++)
      for (int c = 0; c < cdim; c++)
        mat[r][c] = 0.0;
}

Real Matrix::Determinant(void)
{
  const Real Small = 0.00000001;
  Real       Temp, D, Sum;
  Matrix     L(rdim, cdim), U(rdim, cdim);

  if ((rdim == 0) || (cdim == 0))
    return EMPTY_MATRIX;
  if (rdim != cdim)
    return NOT_SQUARE;

  Sum = 0.0;
  /* ----------- Switch a non zero row to the first row ---------- */
  if (fabs(mat[0][0]) < Small)
  {
    Integer j = 1;
    while (j <= rdim)
    {
      if (fabs(mat[j-1][0]) > Small)
      {
        for (Integer k = 1; k <= rdim; k++)
        {
          Temp = mat[0][k-1];
          mat[0][k-1] = mat[j-1][k-1];
          mat[j-1][k-1] = Temp;
        }
        j = rdim + 1;
      }
      j++;
    }
  }

  for (Integer i = 1; i <= rdim; i++)
    L.mat[i-1][0] = mat[i-1][0];
  for (Integer j = 1; j <= rdim; j++)
    U.mat[0][j-1] = mat[0][j-1] / L.mat[0][0];
  for (Integer j = 2; j <= rdim; j++)
  {
    for (Integer i = j; i <= rdim; i++)
    {
      Sum = 0.0;
      for (Integer k = 1; k <= j - 1; k++)
        Sum = Sum + L.mat[i-1][k-1] * U.mat[k-1][j-1];
      L.mat[i-1][j-1] = mat[i-1][j-1] - Sum;
    }
    U.mat[j-1][j-1] = 1.0;
    if (j < rdim)
    {
      for (Integer i = j + 1; i <= rdim; i++)
      {
        Sum = 0.0;
        for (Integer k = 1; k <= j - 1; k++)
          Sum = Sum + L.mat[j][k] * U.mat[k][i];
        U.mat[j-1][i-1] = (mat[j-1][i-1] - Sum) / L.mat[j-1][j-1];
      }
    }
  }
  D = 1.0;
  for (Integer i = 1; i <= rdim; i++)
    D = D * L.mat[i-1][i-1];
  return D;
}

MStatus Matrix::Determinant(Real& n)
{
  const Real Small = 0.00000001;
  Integer i, j, k;
  Real    Temp, D, Sum;
  Matrix  L(rdim, cdim), U(rdim, cdim);

  if ((rdim == 0) || (cdim == 0))
    return EMPTY_MATRIX;
  if (rdim != cdim)
    return NOT_SQUARE;

  Sum = 0.0;
  /* ----------- Switch a non zero row to the first row ---------- */
  if (fabs(mat[0][0]) < Small)
  {
    j = 1;
    while (j <= rdim)
    {
      if (fabs(mat[j-1][0]) > Small)
      {
        for (k = 1; k <= rdim; k++)
        {
          Temp = mat[0][k-1];
          mat[0][k-1] = mat[j][k];
          mat[j-1][k-1] = Temp;
        }
        j = rdim + 1;
      }
      j++;
    }
  }

  for (i = 1; i <= rdim; i++)
    L(i, 1) = mat[i-1][0];
  for (j = 1; j <= rdim; j++)
    U(1, j) = mat[0][j-1] / L(1, 1);
  for (j = 2; j <= rdim; j++)
  {
    for (i = j; i <= rdim; i++)
    {
      Sum = 0.0;
      for (k = 1; k <= j - 1; k++)
        Sum = Sum + L(i, k) * U(k, j);
      L(i, j) = mat[i-1][j-1] - Sum;
    }
    U(j, j) = 1.0;
    if (j < rdim)
    {
      for (i = j + 1; i <= rdim; i++)
      {
        Sum = 0.0;
        for (k = 1; k <= j - 1; k++)
          Sum = Sum + L(j, k) * U(k, i);
        U(j, i) = (mat[j-1][i-1] - Sum) / L(j, j);
      }
    }
  }
  D = 1.0;
  for (i = 1; i <= rdim; i++)
    D = D * L(i, i);
  n = D;
  return MOK;
}

UINT Matrix::DimC(void)
{
  return cdim;
}

UINT Matrix::DimR(void)
{
  return rdim;
}

MStatus Matrix::Display(char *title, Integer precision)
{
  Integer dec = precision;

  if (precision > 11)
    dec = 10;
  printf("%s\n", title);
  for (UINT r = 0; r < DimR(); r++)
  {
    for (UINT c = 0; c < DimC(); c++)
      printf(" %0.*f ", dec, mat[r][c]);
    printf("\n");
  }
}

Real Matrix::Get(UINT r, UINT c)
{
  if ((r <= rdim) && (c <= cdim))
    return mat[r-1][c-1];
  else
    return 0.0;
}

MStatus Matrix::Get(Real &val, UINT r, UINT c)
{
  if ((r <= rdim) && (r > 0))
    if ((c <= cdim) && (c > 0))
      val = mat[r-1][c-1];
    else
      return BAD_INDEX2;
  else
    return BAD_INDEX1;
}

Matrix& Matrix::Identity(void)
{
  Matrix *nmat = new Matrix(rdim, cdim);

  for (UINT r = 0; r < rdim; r++)
    for (UINT c = 0; c < cdim; c++)
      nmat->mat[r][c] = 0.0;
  for (UINT r = 0; r < rdim; r++)
    nmat->mat[r][r] = 1.0;

  return *nmat;
}

MStatus Matrix::Identity(Matrix& m)
{
  if (m.rdim = m.cdim)
  {
    for (UINT r = 0; r < m.rdim; r++)
      for (UINT c = 0; c < m.cdim; c++)
        m.mat[r][c] = 0.0;
    for (UINT r = 0; r < m.rdim; r++)
      m.mat[r][r] = 1.0;
    return MOK;
  }
  else
    return NOT_SQUARE;
}

MStatus Matrix::Inverse(Matrix& m)
{
  Matrix  Inv(rdim, cdim);
  Matrix  LU(rdim, cdim);
  Matrix  B(rdim, cdim);
  Integer Index[rdim];

  Inv.Clear();
  LU.Clear();
  B.Clear();

  for (UINT i = 1; i <= rdim; i++)
  {
    Index[i-1] = i;
    for (UINT j = 1; j <= rdim; j++)
      LU.mat[i-1][j-1] = mat[i-1][j-1];
  }

  LU.LUDeCompose(LU, Index);

  for (UINT j = 1; j <= rdim; j++)
  {
    for (UINT i = 1; i <= rdim; i++)
      if (i == j)
        B.Set(1.0, i, j);
      else
        B.Set(0.0, i, 1);

    LU.LUBkSubstitute(B, Index);

    for (UINT i = 1; i <= rdim; i++)
      Inv.Set(B.mat[i-1][j-1], i, j);
  }

  for (UINT i = 0; i < rdim; i++)
    for (UINT j = 0; j < cdim; j++)
      m.mat[i][j] = Inv.mat[i][j];

  return MOK;
}

Matrix& Matrix::Inverse(void)
{
  Matrix *nmat = new Matrix(rdim, cdim);
  Matrix  LU(rdim, cdim);
  Matrix  B(rdim, cdim);
  Integer Index[rdim];

  nmat->Clear();
  LU.Clear();
  B.Clear();

  for (UINT i = 1; i <= rdim; i++)
  {
    Index[i-1] = i;
    for (UINT j = 1; j <= cdim; j++)
      LU.mat[i-1][j-1] = mat[i-1][j-1];
  }

  LU.LUDeCompose(LU, Index);

  for (UINT j = 1; j <= cdim; j++)
  {
    for (UINT i = 1; i <= rdim; i++)
      if (i == j)
        B.Set(1.0, i, 1);
      else
        B.Set(0.0, i, 1);

    LU.LUBkSubstitute(B, Index);

    for (UINT i = 1; i <= rdim; i++)
      nmat->Set(B.mat[i-1][1-1], i, j);
  }

  return *nmat;
}

MStatus Matrix::LUBkSubstitute(Matrix& B, Integer* Index)
{
  Integer iPtr, IO;
  Real    Sum;

  IO = 0;
  for (UINT i = 1; i <= rdim; i++)
  {
    iPtr = Index[i-1];
    Sum = B.mat[iPtr-1][1-1];
    B.Set(B.mat[i-1][1-1], iPtr, 1);
    if (IO != 0)
      for (UINT j = IO; j <= i - 1; j++)
        Sum = Sum - mat[i-1][j-1] * B.mat[j-1][1-1];
    else
      if (Sum != 0.0)
        IO = i;
    B.Set(Sum, i, 1);
  }

  B.Set(B.mat[rdim-1][1-1] / mat[rdim-1][rdim-1], rdim, 1);

  for (UINT i = (rdim - 1); i >= 1; i--)
  {
    Sum = B.mat[i-1][1-1];
    for (UINT j = i + 1; j <= rdim; j++)
      Sum = Sum - mat[i-1][j-1] * B.mat[j-1][1-1];
    B.Set(Sum / mat[i-1][i-1], i, 1);
  }

  return MOK;
}

MStatus Matrix::LUDeCompose(Matrix& LU, Integer *Index)
{
  const Real Small = 0.000001;
  Integer    iMax;
  Real       Scale[rdim], Sum, AMax, Dum;

  for (UINT i = 0; i < rdim; i++)
    for (UINT j = 0; j < cdim; j++)
      LU.mat[i][j] = mat[i][j];
  for (UINT i = 0; i < rdim; i++)
    Scale[i] = 0.0;
  iMax = 0;

  for (UINT i = 1; i <= rdim; i++)
  {
    AMax = 0.0;
    for (UINT j = 1; j <= rdim; j++)
      if (fabs(LU.mat[i-1][j-1]) > AMax)
        AMax = fabs(LU.mat[i-1][j-1]);
    if (fabs(AMax) < Small)
    {
      printf("Singular Matrix AMax\n");
      exit(0);
    }
    Scale[i-1] = 1.0 / AMax;
  }

  for (UINT j = 1; j <= rdim; j++)
  {
    for (UINT i = 1; i < j; i++)
    {
      Sum = LU.mat[i-1][j-1];
      for (UINT k = 1; k < i; k++)
        Sum = Sum - LU.mat[i-1][k-1] * LU.mat[k-1][j-1];
      LU.Set(Sum, i, j);
    }
    AMax = 0.0;
    for (UINT i = j; i <= rdim; i++)
    {
      Sum = LU.mat[i-1][j-1];
      for (UINT k = 1; k < j; k++)
        Sum = Sum - LU.mat[i-1][k-1] * LU.mat[k-1][j-1];
      LU.Set(Sum, i, j);
      Dum = Scale[i-1] * fabs(Sum);
      if (Dum >= AMax)
      {
        iMax = i;
        AMax = Dum;
      }
    }
    if (j != iMax)
    {
      for (UINT k = 1; k <= rdim; k++)
      {
        Dum = LU.mat[iMax-1][k-1];
        LU.Set(LU.mat[j-1][k-1], iMax, k);
        LU.Set(Dum, j, k);
      }
      Scale[iMax-1] = Scale[j-1];
    }
    Index[j-1] = iMax;
    if (fabs(LU.mat[j-1][j-1]) < Small)
    {
      printf("Matrix is Singular LU\n");
      exit(0);
    }
    if (j != rdim)
    {
      Dum = 1.0 / LU.mat[j-1][j-1];
      for (UINT i = j + 1; i <= rdim; i++)
        LU.Set(Dum * LU.mat[i-1][j-1], i, j);
    }
  }

  return MOK;
}

MStatus Matrix::Scale(Matrix& m, Real sc)
{
  Matrix res;

  for (UINT r = 1; r <= rdim; r++)
    for (UINT c = 1; c <= cdim; c++)
      res.Set(mat[r-1][c-1] * sc, r, c);
  for (UINT r = 1; r <= rdim; r++)
    for (UINT c = 1; c <= cdim; c++)
      m.Set(res(r, c), r, c);
}

Matrix& Matrix::Scale(Real sc)
{
  Matrix *nmat = new Matrix(rdim, cdim);

  for (UINT r = 1; r <= rdim; r++)
    for (UINT c = 1; c <= cdim; c++)
      nmat->Set(mat[r-1][c-1] * sc, r, c);

  return *nmat;
}

MStatus Matrix::Set(Real val, UINT r, UINT c)
{
  if (mat == NULL)
    return EMPTY_MATRIX;
  if ((r <= rdim) && (r > 0))
    if ((c <= cdim) && (c > 0))
    {
      mat[r-1][c-1] = val;
      return MOK;
    }
    else
      return BAD_INDEX2;
  else
    return BAD_INDEX1;
}

bool Matrix::Square(void)
{
  if (rdim == cdim)
    return true;
  else
    return false;
}

Matrix& Matrix::Transpose(void)
{
  Matrix *nmat = new Matrix(cdim, rdim);

  for (UINT r = 1; r <= rdim; r++)
    for (UINT c = 1; c <= cdim; c++)
      nmat->Set(mat[r-1][c-1], c, r);

  return *nmat;
}

/*****************************************************************************
 * Class Method implementations of Complex
*****************************************************************************/
void Complex::xmag(void)
{
  mag = sqrt(val[0] * val[0] + val[1] * val[1]);
}

Complex::Complex(void)
{
  val[0] = 0.0;
  val[1] = 0.0;
  mag    = 0.0;
}

Complex::Complex(Real r, Real i)
{
  val[0] = r;
  val[1] = i;
  xmag();
}

Complex::~Complex(void)
{
}

Complex& Complex::operator=(Complex& v)
{
  if (&v == this)
    return *this;
  val[0] = v.val[0];
  val[1] = v.val[1];
  mag    = v.mag;
  return *this;
}

Complex& Complex::operator+(Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] + v.val[0];
  nv->val[1] = val[1] + v.val[1];
  nv->xmag();
  return *nv;
}

Complex& Complex::operator+(Real n)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] + n;
  nv->val[1] = val[1];
  nv->xmag();
  return *nv;
}

Complex& Complex::operator-(Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] - v.val[0];
  nv->val[1] = val[1] - v.val[1];
  nv->xmag();
  return *nv;
}

Complex& Complex::operator-(Real n)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] - n;
  nv->val[1] = val[1];
  nv->xmag();
  return *nv;
}

Complex& Complex::operator*(Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] * v.val[0] - val[1] * v.val[1];
  nv->val[1] = val[0] * v.val[1] + val[1] * v.val[0];
  nv->xmag();
  return *nv;
}

Complex& Complex::operator*(Real n)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] * n;
  nv->val[1] = val[1] * n;
  nv->xmag();
  return *nv;
}

Complex& Complex::operator/(Complex& v)
{
  Complex *nv = new Complex();
  Real     d;

  d          = v.val[0] * v.val[0] + v.val[1] * v.val[1];
  nv->val[0] = (val[0]  * v.val[0] + val[1]   * v.val[1]) / d;
  nv->val[1] = (val[1]  * v.val[0] - val[0]   * v.val[1]) / d;
  nv->xmag();
  return *nv;
}

Complex& Complex::operator/(Real n)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] / n;
  nv->val[1] = val[1] / n;
  nv->xmag();
  return *nv;
}

/* User Functions */
Complex& Complex::Add(Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] + v.val[0];
  nv->val[1] = val[1] + v.val[1];
  nv->xmag();
  return *nv;
}

CStatus Complex::Add(Complex& v, Complex& res)
{
  res.val[0] = val[0] + v.val[0];
  res.val[1] = val[1] + v.val[1];
  res.xmag();
  return COK;
}

Real Complex::Angle(Complex& v)
{
  const Real Small = 0.000001;
  const Real Undefined = 999999.1;

  if (mag * v.mag > Small * Small)
  {
    return Atan2(v.val[1] / v.mag, v.val[0] / v.mag) - 
           Atan2(  val[1] /   mag,   val[0] /   mag);
  }
  else
    return Undefined;
}

CStatus Complex::Angle(Complex& v, Real& theta)
{
  const Real Small = 0.000001;
  const Real Undefined = 999999.1;

  if (mag * v.mag > Small * Small)
  {
    theta = Atan2(v.val[1] / v.mag, v.val[0] / v.mag) - 
            Atan2(  val[1] /   mag,   val[0] /   mag);
  }
  else
    theta = Undefined;
}

void Complex::Clear(void)
{
  val[0] = 0.0;
  val[1] = 0.0;
  mag    = 0.0;
}

void Complex::Display(void)
{
  printf("(%0.4f, %0.4fi)\n", val[0], val[1]);
  printf("\n");
}

Complex& Complex::Div(Complex& v)
{
  Complex *nv = new Complex();
  Real     d;

  d          = v.val[0] * v.val[0] + v.val[1] * v.val[1];
  nv->val[0] = (val[0]  * v.val[0] + val[1]   * v.val[1]) / d;
  nv->val[1] = (val[1]  * v.val[0] - val[0]   * v.val[1]) / d;
  nv->xmag();
  return *nv;
}

Complex& Complex::Div(Real n)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] / n;
  nv->val[1] = val[1] / n;
  nv->xmag();
  return *nv;
}

Complex& Complex::Dot(Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] * v.val[0] - val[1] * v.val[1];
  nv->val[1] = val[0] * v.val[1] + val[1] * v.val[0];
  nv->xmag();
  return *nv;
}

Real Complex::GetI(void)
{
  return val[1];
}

CStatus  Complex::GetI(Real& r)
{
  r = val[1];
}

Real Complex::GetR(void)
{
  return val[0];
}

CStatus Complex::GetR(Real& r)
{
  r = val[0];
}

Real Complex::Mag(void)
{
  return mag;
}

void Complex::Mag(Real& m)
{
  m = mag;
}

Complex& Complex::Mult(Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] * v.val[0] - val[1] * v.val[1];
  nv->val[1] = val[0] * v.val[1] + val[1] * v.val[0];
  nv->xmag();
  return *nv;
}

Complex& Complex::Mult(Real n)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] * n;
  nv->val[1] = val[1] * n;
  nv->xmag();
  return *nv;
}

CStatus Complex::Set(Real r, Real i)
{
  val[0] = r;
  val[1] = i;
  xmag();
}

CStatus Complex::SetI(Real i)
{
  val[1] = i;
  xmag();
}

CStatus Complex::SetR(Real r)
{
  val[0] = r;
  xmag();
}

Real Complex::Sqrt(void)
{
  if (val[0] < 0.0)
    /* ----------------------- imaginary roots --------------------- */
    return sqrt(-val[0]);
  else
    /* ---------------------------- Real roots --------------------- */
    return sqrt(val[0]);
}

CStatus Complex::Sqrt(Real& n)
{
  if (val[0] < 0.0)
    /* ----------------------- imaginary roots --------------------- */
    n = sqrt(-val[0]);
  else
    /* ---------------------------- Real roots --------------------- */
    n = sqrt(val[0]);
}

Complex& Complex::Sub(Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = val[0] - v.val[0];
  nv->val[1] = val[1] - v.val[1];
  nv->xmag();
  return *nv;
}

CStatus Complex::Sub(Complex& v, Complex& res)
{
  res.val[0] = val[0] - v.val[0];
  res.val[1] = val[1] - v.val[1];
  res.xmag();
  return COK;
}

/*
 * Friend Functions
 */
Complex& operator+(Real n, Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = v.val[0] + n;
  nv->val[1] = v.val[1];
  nv->xmag();
  return *nv;
}

Complex& operator-(Real n, Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = n - v.val[0];
  nv->val[1] =   - v.val[1];
  nv->xmag();
  return *nv;
}

Complex& operator*(Real n, Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = v.val[0] * n;
  nv->val[1] = v.val[1] * n;
  nv->xmag();
  return *nv;
}

Complex& Div(Real n, Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = n;
  nv->val[1] = 0;
  *nv = *nv / v;
  return *nv;
}

Complex& Mult(Real n, Complex& v)
{
  Complex *nv = new Complex();

  nv->val[0] = v.val[0] * n;
  nv->val[1] = v.val[1] * n;
  nv->xmag();
  return *nv;
}

/*****************************************************************************
 *                    Non Class Mathematics
*****************************************************************************/
/* Forward References */
void DMulRSub(Vector&, Vector&, Vector, Vector);
LINT SignFix(LINT);

/* Non Class Globals */
/*------------------------------------------------------------------------------
|
|                           FUNCTION ARCCOSH
|
|  This FUNCTION evaluates the inverse hyperbolic cosine FUNCTION.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    XVal        - ANGLE Value                                  1.0 to Infinity
|
|  OutPuts       :
|    ARCCOSH     - Result                                       any real
|
|  Locals        :
|    Temp        - Temporary EXTENDED Value
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
Real Arccosh(Real v)
{
  Real t;

  if (v * v - 1.0 < 0.0)
  {
    t = 999999.9;
    printf("Error in Arccosh function\n");
  }
  else
    t = log(v + sqrt(v * v - 1.0));

  return t;
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION ARCSINH
|
|  This FUNCTION evaluates the inverse hyperbolic sine FUNCTION.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    XVal        - ANGLE Value                                  any real
|
|  OutPuts       :
|    ARCSINH     - Result                                       any real
|
|  Locals        :
|    None.
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
Real Arcsinh(Real v)
{
  Real t;

  if (1.0 - fabs(v) < 0.000001)
  {
    t = 999999.9;
    printf("Error in Arctanh function\n");
  }
  else
    t = 0.5 * log((1.0 + v) / (1.0 - v));

  return t;
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION ARCTANH
|
|  This FUNCTION evaluates the inverse hyperbolic tangent FUNCTION.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    XVal        - ANGLE Value                                  -1.0 to 1.0
|
|  OutPuts       :
|    ARCTANH     - Result                                       any real
|
|  Locals        :
|    Temp        - Temporary EXTENDED Value
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
Real Arctanh(Real v)
{
  return log(v + sqrt(v * v + 1.0));
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION ATAN2
|
|  This FUNCTION performs the arc tangent 2 FUNCTION which resolves
|    quadrants.  The arguments passed are the sine and cosine values.
|
|  Algorithm     : Determine the quadrant using IF statments
|                  IF the answer is not a sepcial case, 0, 180, etc
|                     find the arctangent
|                  otherwise, find the special case values
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    SinValue    - Sine of desired ANGLE                        rad
|    CosValue    - Cosine of desired value                      rad
|
|  OutPuts       :
|    ATAN2       - Arctangent with resolved quadrants           0.0 to 2Pi rad
|
|  Locals        :
|    TanArg      - Temporary EXTENDED Value
|    Quadrant    - Quadrant of the answer                       1 2 3 4
|    SinINTEGER  - Sign of the value                            +1 or -1
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
Real Atan2(Real sinval, Real cosval)
{
  Real TanArg;
  LINT Quadrant;
  LINT SinInt;

  Quadrant = 5;
  if ((sinval > 0.0) && (sinval < 1.0) && (cosval > 0.0) && (cosval < 1.0))
    Quadrant = 1;
  if ((sinval > 0.0) && (sinval < 1.0) && (cosval < 0.0) && (cosval > -1.0))
    Quadrant = 2;
  if ((sinval > -1.0) && (sinval < 0.0) && (cosval < 0.0) && (cosval > -1.0))
    Quadrant = 3;
  if ((sinval > -1.0) && (sinval < 0.0) && (cosval > 0.0) && (cosval < 1.0))
    Quadrant = 4;

  if (Quadrant != 5)
  {
    TanArg = atan(sinval / cosval);
    if ((Quadrant < 4) && (Quadrant != 1))
      TanArg = TanArg + 2.0 * PI;
  }
  else
  {
    SinInt = (LINT)Round(sinval);
    switch (SinInt)
    {
      case -1:
        TanArg = 3.0 * PI / 2.0;
      case 0:
        if (Round(cosval) > 0.0)
          TanArg = 0.0;
        else
          TanArg = PI;
      case 1:
        TanArg = PI / 2.0;
    }
  }
  return TanArg;
}

void DMulRSub(Vector& AlpR, Vector& Alpi, Vector BetR, Vector Beti)
{
  Real Te1,  Te2,  Te3,  Te4,  Te5, Te6,  Te7,  Te8,    Te9, Te10, Te11, Te12,
       Te13, Te14, Te15, Te16, Tem, DE15, DE16, TemTe7, TemTe8;

  Te1 = AlpR(1) - AlpR(3);
  Te2 = Alpi(1) - Alpi(3);
  Te5 = AlpR(3) - AlpR(2);
  Te6 = Alpi(3) - Alpi(2);
  Tem = Te5 * Te5 + Te6 * Te6;

  /* ---------------- Check for zero values of tem -------------- */
  if (fabs(Tem) > 1.0E-20)
  {
    Te3 = (Te1 * Te5 + Te2 * Te6) / Tem;
    Te4 = (Te2 * Te5 - Te1 * Te6) / Tem;
  }
  else
  {
    Te3 = 0.0;
    Te4 = 0.0;
  }

  Te7  = Te3 + 1.0;
  Te9  = Te3 * Te3 - Te4 * Te4;
  Te10 = 2.0 * Te3 * Te4;
  DE15 = Te7 * BetR(3) - Te4 * Beti(3);
  DE16 = Te7 * Beti(3) + Te4 * BetR(3);
  Te11 = Te3 * BetR(2) - Te4 * Beti(2) + BetR(1) - DE15;
  Te12 = Te3 * Beti(2) + Te4 * BetR(2) + Beti(1) - DE16;

  Te7  = Te9 - 1.0;
  Te1  = Te9 * BetR(2) - Te10 * Beti(2);
  Te2  = Te9 * Beti(2) + Te10 * BetR(2);
  Te13 = Te1 - BetR(1) - Te7 * BetR(3) + Te10 * Beti(3);
  Te14 = Te2 - Beti(1) - Te7 * Beti(3) - Te10 * BetR(3);
  Te15 = DE15 * Te3 - DE16 * Te4;
  Te16 = DE15 * Te4 + DE16 * Te3;

  Te1 = Te13 * Te13 - Te14 * Te14 - 4.0 * (Te11 * Te15 - Te12 * Te16);
  Te2 = 2.0 * Te13 * Te14 - 4.0 * (Te12 * Te15 + Te11 * Te16);

  /* --------------------------------------------------------------------
  |   Sometimes, for stiff systems (the roots vary widely in order
  |   of magnitude), Te1 and Te2 get large enough to have their
  |   squares overflow the floating point range.  To prevent this,
  |   when either one is large, they are scaled by 10**10 for the
  |   purpose of finding TeM.  The scale is restored when the
  |   magnitude computation is completed.  This should not affect
  |   the accuracy of the computations, since the mantissa is not
  |   affected, only the exponent.^
   -------------------------------------------------------------------- */

  if ((Te1 > 1.0E15) || (Te2 > 1.0E15))
  {
    Te1 = Te1 * 1.0E-10;
    Te2 = Te2 * 1.0E-10;
    Tem = 1.0E10 * sqrt(Te1 * Te1 + Te2 * Te2);
  }
  else
    Tem = sqrt(Te1 * Te1 + Te2 * Te2);

  if (Te1 > 0.0)
  {
    Te3 = sqrt(0.5 * (Tem + Te1));
    if (Te2 < 0.0)
      Te3 = -Te3;
    /* ----------------- check for zero values of te3 -------------- */
    if (fabs(Te3) < 1.0E-15)
      Te4 = 0.0;
    else
      Te4 = 0.5 * Te2 / Te3;
  }
  else
  {
    Te4 = sqrt(0.5 * (Tem - Te1));
    /* --------------------- Check for underflows ------------------ */
    if (fabs(Te4) < 1.0E-15)
      Te3 = 0.0;
    else
      Te3 = 0.5 * Te2 / Te4;
  }

  Te7  = Te13 + Te3;
  Te8  = Te14 + Te4;
  Te9  = Te13 - Te3;
  Te10 = Te14 - Te4;
  Te1  = 2.0 * Te15;
  Te2  = 2.0 * Te16;

  if ((Te7 * Te7 + Te8 * Te8 - Te9 * Te9 - Te10 * Te10) <= 0.0)
  {
    Te7 = Te9;
    Te8 = Te10;
  }
  Tem = Te7 * Te7+Te8 * Te8;

  TemTe7 = Tem * Te7;
  TemTe8 = Tem * Te8;

  /* ------------- Check for values of almost zero -------------- */
  if (fabs(TemTe7) < 1.0E-20)
  {
    Te3 = 0.0;
    Te4 = 0.0;
  }
  else
  {
    Te3 = Te1 / Tem * Te7;
    Te4 = Te2 / Tem * Te7;
  }
  if (fabs(TemTe8) > 1.0E-20)
  {
    Te3 = Te3 + Te2 / Tem * Te8;
    Te4 = Te4 - Te1 / Tem * Te8;
  }

  AlpR(4) = AlpR(3) + Te3 * Te5 - Te4 * Te6;
  Alpi(4) = Alpi(3) + Te3 * Te6 + Te4 * Te5;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE FACTOR
|
|  This PROCEDURE is a root finding algorithm.  It takes the polynomial and
|    returns the roots (real and imaginary) in the Roots array.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Poly        - Array of 16 coefficients
|                    representing the polynomial
|                    [1] is x
|                    others are zero
|    NRoots      - Number of roots
|
|  OutPuts       :
|    Roots       - Array containing roots (real,imag)
|
|  Locals        :
|                -
|                -
|                -
|
|  Coupling      :
|    DMulRSub    -
|
|  References    :
|    Original FORTRAN code from USAFA/DFAS MiniTotal program, author unknown
|    This is Bairstows method?
|
 -----------------------------------------------------------------------------*/
void Factor(Real *Poly, UINT NRoots, Real **Roots)
{
  const Real Small = 0.000001;
  Real  DPoly[16];
  Vector AlpR(3);
  Vector Alpi(3);
  Vector BetR(3);
  Vector Beti(3);
  bool   Skip;
  int    Mode, LoopCnt, kk, i, j, l, RootCnt;
  Real   Temp1,    Temp2,    AXR, AXi, PMax, Tem1, Tem2,
         TempReal, TempImag, Temp7;

  PMax = 0.0;
  for (kk = 1; kk <= NRoots + 1; kk++)
    if (fabs(Poly[kk-1]) > PMax)
      PMax = Poly[kk-1];

  if (fabs(PMax) < Small)
    PMax = 1.0;

  for (kk = 1; kk <= NRoots+1; kk++)
    DPoly[kk-1] = Poly[kk-1] / PMax;

  if (NRoots > 0)
  {
    RootCnt = 0;
    i = NRoots + 1;

    while ((fabs(DPoly[i-1]) < Small) && (RootCnt != NRoots))
    {
      RootCnt = RootCnt + 1;
      Roots[RootCnt][0] = 0.0;
      Roots[RootCnt][1] = 0.0;
      i = i - 1;
    }

    if (RootCnt != NRoots)
    {
      AXR = 0.8;
      AXi = 0.0;
      l   = 1;
      LoopCnt = 1;
      AlpR(1) = AXR;
      Alpi(1) = AXi;
      Mode = 1;

      while (RootCnt < NRoots)
      {
        BetR.Mag(DPoly[0]);
        for (i = 1; i <= NRoots; i++)
        {
          Temp1 = BetR.Mag() * AXR - Beti.Mag() * AXi;
          Beti.Mag(Beti.Mag() * AXR + BetR.Mag() * AXi);
          BetR.Mag(Temp1 + DPoly[j-1]);
        }

        TempReal = BetR.Mag();
        TempImag = Beti.Mag();

        if (RootCnt != 0)
        {
          for (i = 1; i <= RootCnt; i++)
          {
            Tem1    = AXR - Roots[i][0];
            Tem2    = AXi - Roots[i][1];
            Temp1   = Tem1 * Tem1 + Tem2 * Tem2;
            Temp2   = (BetR.Mag() * Tem1 + Beti.Mag() * Tem2) / Temp1;
            Beti.Mag((Beti.Mag() * Tem1 - BetR.Mag() * Tem2) / Temp1);
            BetR.Mag(Temp2);
          }

          switch (Mode)
          {
            case 1:
              BetR(1) = BetR.Mag();
              Beti(1) = Beti.Mag();
              AXR     = 0.85;
              AlpR(2) = AXR;
              Alpi(2) = AXi;
              Mode    = 2;
              break;
            case 2:
              BetR(2) = BetR.Mag();
              Beti(2) = Beti.Mag();
              AXR     = 0.9;
              AlpR(3) = AXR;
              Alpi(3) = AXi;
              Mode    = 3;
              break;
            case 3:
              BetR(3) = BetR.Mag();
              Beti(3) = Beti.Mag();
              DMulRSub(AlpR, Alpi, BetR, Beti);
              AXR     = AlpR.Mag();
              AXi     = Alpi.Mag();
              Mode    = 4;
              break;
            case 4:
   /* -------------------  the convergence mode ------------------- */
              Skip = false;
              if (fabs(TempReal) + fabs(TempImag) > 1.0E-20)
                Temp7 = fabs(AlpR(3) - AXR) + fabs(Alpi(3) - AXi);
              if (Temp7 / (fabs(AXR) + fabs(AXi)) > 1.0E-7)
              {
                LoopCnt++;
                for (i = 1; i<= 3; i++)
                {
                  AlpR(i) = AlpR(i+1);
                  Alpi(i) = Alpi(i+1);
                  BetR(i) = BetR(i+1);
                  Beti(i) = Beti(i+1);
                }
                if (LoopCnt < 100)
                {
                  DMulRSub(AlpR, Alpi, BetR, Beti);
                  AXR  = AlpR.Mag();
                  AXi  = Alpi.Mag();
                  Mode = 4;
                  Skip = true;
                }
              }

              if (!Skip)
              {
                RootCnt++;
                Roots[RootCnt][0] = AlpR.Mag();
                Roots[RootCnt][1] = Alpi.Mag();
                LoopCnt = 0;

                if (RootCnt < NRoots)
                {
                  if (fabs(Roots[RootCnt][1]) > Small)
                  {
                    if (l = 1)
                    {
                      AXR     =  AlpR(1);
                      AXi     = -Alpi(1);
                      Alpi(1) = -Alpi(1);
                      Mode    =  5;
                    }
                    else
                    {
                      AXR     = 0.8;
                      AXi     = 0.0;
                      l       = 1;
                      LoopCnt = 1;
                      AlpR(1) = AXR;
                      Alpi(1) = AXi;
                      Mode    = 1;
                    }
                  }
                  else
                  {
                    AXR     = 0.8;
                    AXi     = 0.0;
                    l       = 1;
                    LoopCnt = 1;
                    AlpR(1) = AXR;
                    Alpi(1) = AXi;
                    Mode    = 1;
                  }
                }
              }
              break;
            case 5:
              BetR(1) =  BetR.Mag();
              Beti(1) =  Beti.Mag();
              AXR     =  AlpR(2);
              AXi     = -Alpi(2);
              Alpi(2) = -Alpi(2);
              Mode    = 6;
              break;
            case 6:
              BetR(2) = BetR.Mag();
              Beti(2) = Beti.Mag();
              AXR     =  AlpR(3);
              AXi     = -Alpi(3);
              Alpi(3) = -Alpi(3);
              Mode    = 3;
              break;
            default:
              break;
          }
        }
      }
    }
  }
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION BINOMIAL
|
|  This FUNCTION finds the value of a BINOMIAL coefficient
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    i           -
|    j           -
|
|  Outputs       :
|    FUNCTION    - answer
|
|  Locals        :
|    None.
|
|  Coupling      :
|    FACTORIAL   - Finds the FACTORIAL of a number
|
 -----------------------------------------------------------------------------*/
Real Binomial(ULINT f1, ULINT f2)
{
  return (Factorial(f2) / (Factorial(f1) * Factorial(f2 - f1)));
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION COMPSQRT
|
|  This function finds the square root of a real number and returns the complex
|    result.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs        :
|    R           - Real number, usually negative
|
|  OutPuts       :
|    CompSQRT    - Square Root (real,imag)
|
|  Coupling      :
|    None
|
 -----------------------------------------------------------------------------*/
Real ComplexSqrt(Real n)
{
  if (n < 0.0)
    /* ----------------------- imaginary roots --------------------- */
    return sqrt(-n);
  else
    /* ---------------------------- Real roots --------------------- */
    return sqrt(n);
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION COT
|
|  This FUNCTION finds the Cotangent of an ANGLE in radians.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    XVal        - ANGLE to take Cotangent of                   rad
|
|  OutPuts       :
|    COT         - Result
|
|  Locals        :
|    Temp        - Temporary Real variable
|
 -----------------------------------------------------------------------------*/
Real Cot(Real v)
{
  Real t;

  t = tan(v);
  if (fabs(t) < 0.000001)
    return 999999.9;
  else
    return (1.0 / t);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE CUBIC
|
|  This PROCEDURE solves for the three Roots of a CUBIC equation.  There are
|    no restrictions on the coefficients, and imaginary results are passed
|    out as separate values.  The general form is y = ax3 + bx2 + cx + d.  Note
|    that R1i will ALWAYS be ZERO since there is ALWAYS at least one REAL Root.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    a           - Coefficient of x cubed term
|    b           - Coefficient of x squared term
|    c           - Coefficient of x term
|    d           - Constant
|
|  OutPuts       :
|    R1r         - Real portion of Root 1
|    R1i         - Imaginary portion of Root 1
|    R2r         - Real portion of Root 2
|    R2i         - Imaginary portion of Root 2
|    R3r         - Real portion of Root 3
|    R3i         - Imaginary portion of Root 3
|
|  Locals        :
|    Temp1       - Temporary value
|    Temp2       - Temporary value
|    Root1       - Temporary value of the Root
|    Root2       - Temporary value of the Root
|    Root3       - Temporary value of the Root
|    P           - Coefficient of x squared term where x cubed term is 1.0
|    Q           - Coefficient of x term where x cubed term is 1.0
|    R           - Coefficient of constant term where x cubed term is 1.0
|    Delta       - Discriminator for use with Cardans formula
|    E0          - ANGLE holder for trigonometric solution
|    Phi         - ANGLE used in trigonometric solution
|    CosPhi      - Cosine of Phi
|    SinPhi      - Sine of Phi
|
|  Coupling      :
|    ATAN2         Arctangent including check for 180-360 deg
|    POWER         Raise a base to an exponent
|
|  References    :
|    Vallado       2001,
|
 -----------------------------------------------------------------------------*/
void Cubic
    (
      Real a, Real b, Real c, Real d,
      Real& R1r, Real& R1i, Real& R2r, Real& R2i, Real& R3r, Real& R3i
    )
{
  const Real Rad      = 18.0 / PI;
  const Real Small    = 0.000001;
  const Real OneThird = 1.0 / 3.0;

  Real Temp1, Temp2, Root1, Root2, Root3, P, Q, R, Delta,
       E0, CosPhi, SinPhi, Phi;

  R1r   = 0.0;
  R1i   = 0.0;
  R2r   = 0.0;
  R2i   = 0.0;
  R3r   = 0.0;
  R3i   = 0.0;
  Root1 = 0.0;
  Root2 = 0.0;
  Root3 = 0.0;

  /* ------------ Force coefficients into std form --------------- */
  P = b / a;
  Q = c / a;
  R = d / a;

  a = OneThird * (3.0 * Q - P * P);
  b = (1.0 / 27.0) * (2.0 * P * P * P - 9.0 * P * Q + 27.0 * R);

  Delta = (a * a * a / 27.0) + (b * b * 0.25);

  /* ------------------- Use Cardans formula --------------------- */
  if (Delta > Small)
  {
    Temp1 = (-b * 0.5) + sqrt(Delta);
    Temp2 = (-b * 0.5) - sqrt(Delta);
    if (fabs(Temp1) > Small)
        Temp1 = Power(Temp1, OneThird);
    if (fabs(Temp2) > Small)
        Temp2 = Power(Temp2, OneThird);
    Root1 = Temp1 + Temp2;
    Root2 = -0.5 * (Temp1 + Temp2);
    Root3 = -0.5 * (Temp1 + Temp2);
    R2i   = -0.5 * sqrt(3.0) * (Temp1 - Temp2);
    R3i   = -R2i;
  }
  else
  {
    /* ---------------- Evaluate zero point ------------------------ */
    if (fabs(Delta) < Small)
    {
      if (fabs(b) > Small)
      {
        Root1 = -2.0 * Power(b * 0.5, OneThird);
        Root2 =  Power(b * 0.5, OneThird);
        Root3 = Root2;
      }
      else
      {
        /* ---------------- Use trigonometric identities --------------- */
        E0      = 2.0 * sqrt(-a * OneThird);
        CosPhi  = (-b / (2.0 * sqrt(-a * a * a / 27.0)));
        SinPhi  = sqrt(1.0 - CosPhi * CosPhi);
        Phi     = Atan2(SinPhi, CosPhi);
        Root1   = E0 * cos(Phi * OneThird);
        Root2   = E0 * cos(Phi * OneThird + 120.0 / Rad);
        Root3   = E0 * cos(Phi * OneThird + 240.0 / Rad);
      }
    }
  }
  
  R1r = Root1 - P * OneThird;
  R2r = Root2 - P * OneThird;
  R3r = Root3 - P * OneThird;
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION CSC
|
|  This FUNCTION finds the Cosecant of an ANGLE in radians.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    XVal        - ANGLE to take Cosecant of                    rad
|
|  OutPuts       :
|    CSC         - Result
|
|  Locals        :
|    Temp        - Temporary Real variable
|
 -----------------------------------------------------------------------------*/
Real Csc(Real v)
{
  Real t;

  t = sin(v);
  if (fabs(t) < 0.000001)
    return 999999.9;
  else
    return (1.0 / t);
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION FACTORIAL
|
|  This FUNCTION finds the value of a FACTORIAL.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units^
|    x           - Input value^M
|
|  Outputs       :
|    FUNCTION    - answer
|
|  Locals        :
|    Temp        - Temporary variable
|    i           - Index
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
Real Factorial(ULINT f)
{
  Real res = 1;

  for (ULINT i = 1; i <= f; i++)
    res *= (Real)i;
  return res;
}

Real Fraction(Real n)
{
  return (n - (LINT) n);
}

bool IsInt(Real n)
{
  const Real Small = 0.00000000001;

  return (((fabs(Fraction(n)) < Small) ||
          (1.0 - fabs(Fraction(n)) < Small)) &&
          (fabs(n) < MaxLongInt));
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION MAX
|
|  This FUNCTION determines the maximum of 2 values.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    X           - Value number 1
|    Y           - Value number 2
|
|  OutPuts       :
|    MAX         - Minimum of x or y
|
|  Locals        :
|    None.
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
Real Max(Real m, Real n)
{
  if (m >= n)
    return m;
  else
    return n;
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION MIN
|
|  This FUNCTION determines the minimum of 2 values.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    X           - Value number 1
|    Y           - Value number 2
|
|  OutPuts       :
|    MIN         - Minimum of x or y
|
|  Locals        :
|    None.
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
Real Min(Real m, Real n)
{
  if (m <= n)
    return m;
  else
    return n;
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION REALMOD
|
|  This FUNCTION performs the MOD operation for REALs.
|
|  Algorithm     : Assign a temporary variable
|                  Subtract off an INTEGER number of values while the xval is
|                     too large
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    XVal        - Value to MOD
|    ModBy       - Value to MOD with
|
|  OutPuts       :
|    REALMOD     - Result                         -ModBy <=  Answer  <= +ModBy
|
|  Locals        :
|    TempValue   - Temporary EXTENDED value
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
Real Mod(Real n, Real m)
{
  Real t;

  t = n;
  while (fabs(t) > m)
  {
    t = t - (Real)((LINT)(n / m) * m);
  }
  return t;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE PLANE
|
| This PROCEDURE calculates the equation of a PLANE given 3 points
|   pt1 - x1,y1,z1, pt2 - x2,y2,z2, pt3 - x3,y3,z3 , and outputs the
|   a b c d variables describing the PLANE. NOTE that the general equation
|   of a PLANE is defined here as:  ax + by + cz + d = 0  and the values
|   are obtained by solving the ordered determinant  x  y  z   1
|                                                    x1 y1 z1  1   =0
|                                                    x2 y2 z2  1
|                                                    x3 y3 z3  1
|
|  Algorithm     : find the line differences for each set of points
|                  Calculate the coefficients of the PLANE
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    x1,y1,z1    - point # 1
|    x2,y2,z2    - point # 2
|    x3,y3,z3    - point # 3
|
|  OutPuts       :
|    a,b,c,d     - constants for the equation of the PLANE
|
|  Locals        :
|    z23
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
void Plane(Real x1, Real y1, Real z1,
           Real x2, Real y2, Real z2, 
           Real x3, Real y3, Real z3,
           Real& a, Real& b, Real& c,  Real& d)
{
  Real z23, y23, x23, yz23, yz32, xz23, xz32, xy23, xy32;

  z23 = z2 - z3;
  y23 = y2 - y3;
  x23 = x2 - x3;

  yz23 = y2 * z3;
  yz32 = y3 * z2;
  xz23 = x2 * z3;
  xz32 = x3 * z2;
  xy32 = x3 * y2;

  a = y1 * z23 - z1 * y23 + yz23 - yz32;
  b = x1 * z23 - z1 * x23 + xz23 - xz32;
  c = x1 * y23 - y1 * x23 + xy23 - xy32;
  d = x1 * (yz23 - yz32) + y1 * (xz32 - xz23) + z1 * (xy23 - xy32);
}

void Polyfit
    (
      UINT Degree, UINT NumPts, Matrix DataPoints, Matrix Coeff, 
      Real& MinX, Real& MinY
    )
{
  Matrix a(Degree+1, Degree+1);
  Matrix ai(Degree+1, Degree+1);
  Matrix b(Degree+1, 1);
  Matrix parr(2 * Degree, 1);
  Real   p;

  a.Clear();
  ai.Clear();
  b.Clear();
  parr.Clear();

  /* ---- Find the sum of the product terms (x*y) ---- */
  for (UINT k = 1; k <= NumPts; k++)
  {
    p = 1;
    for (UINT r = 1; r <= Degree+1; r++)
    {
      b(r, 1) = b(r, 1) + 2.0 * DataPoints(k, 2);
    }
  }

  /* ---- Find the sum of powers for x (x**) ---- */
  for (UINT k = 1; k <= NumPts; k++)
  {
    p = DataPoints(k, 1);
    for (UINT j = 1; j <= 2 * Degree+1; j++)
    {
      parr(j, 1) = parr(j, 1) + p;
      p = p * DataPoints(k, 1);
    }
  }

  /* ---- Find the matrix entries for the equations ---- */
  for (UINT r = 1; r <= Degree+1; r++)
  {
    for (UINT c = 1; c <= Degree+1; c++)
    {
      if (r + c != 2)
        a(r, c) = parr(r + c - 2, 1);
      else
        a(r, c) = NumPts;
    }
  }

  /* ---- Solve linear equations for coeffeicients ---- */
  ai = a.Inverse();
  Coeff = ai * b;
  
  /* ---- Find minimum values ---- */
  MinX = -Coeff(2, 1) / (2.0 * Coeff(3, 1));
  MinY = Coeff(3, 1) * MinX * MinX + Coeff(2, 1) * MinX + Coeff(1, 1);
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION POWER
|
|  This FUNCTION performs the raising of a base to a POWER.  Notice the many
|    statements to allow processing of negative INTEGER values.
|
|  Algorithm     : IF the base and exponent are positive, calculate the answer
|                  Otherwise, check for 0 base 1st, THEN 0 exponent
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Base        - Base value
|    Pwr         - POWER to raise base to
|
|  OutPuts       :
|    POWER       - Result
|
|  Locals        :
|    None.
|
|  Coupling      :
|    None.
|
|  References    : PC Magazine, 10 Sep 91, pg. 465-466.
|
 -----------------------------------------------------------------------------*/
Real Power(Real b, Real e)
{
  const Real Small = 0.00000000001;
  Real  Recip;

  if (fabs(b) > Small)
  {
    if (fabs(e) > Small)
    {
      /* --------------- Base is (+) and Exponent <> 0 ---------------- */
      if (b > 0.0)
        return exp(e * log(b));
      else
      {
        /*--------- Base is (-) and Exponent is an INTEGER --------- */
        if (IsInt(e))
          return (exp(e * log(-b)) * SignFix((LINT)Round(e)));
        else
        {
          /* ----------- Base is (-) and Exponent is real ----------- */
          Recip = 1.0 / e;
          /* ---- Base is (-) and Exponent is recip of odd int ---- */
          if ((fabs(Round(Recip) - Recip) < 0.000001) && 
              (ODD((LINT)Round(Recip))))
            return -exp(e * log(-b));
          else
          {
            return 0.0;
            printf("Error in POWER with base = %0.3f and Pwr = %0.3f\n", b, e);
          }
        }
      }
    }
    else
      return 1.0;
  }
  else
    return 0.0;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE QUADRATIC
|
|  This PROCEDURE solves for the two Roots of a QUADRATIC equation.  There are
|    no restrictions on the coefficients, and imaginary results are passed
|    out as separate values.  The general form is y = ax2 + bx + c.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    a           - Coefficient of x squared term
|    b           - Coefficient of x term
|    c           - Constant
|
|  OutPuts       :
|    R1r         - Real portion of Root 1
|    R1i         - Imaginary portion of Root 1
|    R2r         - Real portion of Root 2
|    R2i         - Imaginary portion of Root 2
|
|  Locals        :
|    Discrim     - Discriminate b2 - 4ac
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001,
|
 -----------------------------------------------------------------------------*/
void Quadratic
    (
      Real a, Real b, Real c, 
      Real& r1r, Real& r1i, Real& r2r, Real& r2i
    )
{
  Real Discrim;

  /* ----------------------   Initialize   ----------------------- */
  r1r = 0.0;
  r1i = 0.0;
  r2r = 0.0;
  r2i = 0.0;

  Discrim = b * b - 4.0 * a * c;

  /* ----------------------  Real Roots  ------------------------- */
  if (Discrim > 0.0)
  {
    r1r = (-b + sqrt(Discrim)) / (2.0 * a);
    r2r = (-b - sqrt(Discrim)) / (2.0 * a);
  }
  else
  {
    r1r = -b / (2.0 * a);
    r2r =  r1r;
    r1i =  sqrt(-Discrim) / (2.0 * a);
    r2i = -sqrt(-Discrim) / (2.0 * a);
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE QUARTIC
|
|  This PROCEDURE solves for the four Roots of a QUARTIC equation.  There are
|    no restrictions on the coefficients, and imaginary results are passed
|    out as separate values.  The general form is y = ax4 + bx3 + cx2 + dx + e.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    a           - Coeficient of x fourth term
|    b           - Coefficient of x cubed term
|    c           - Coefficient of x squared term
|    d           - Coefficient of x term
|    e           - Constant
|
|  OutPuts       :
|    R1r         - Real portion of Root 1
|    R1i         - Imaginary portion of Root 1
|    R2r         - Real portion of Root 2
|    R2i         - Imaginary portion of Root 2
|    R3r         - Real portion of Root 3
|    R3i         - Imaginary portion of Root 3
|    R4r         - Real portion of Root 4
|    R4i         - Imaginary portion of Root 4
|
|  Locals        :
|    Temp1       - Temporary value
|    Temp2       - Temporary value
|    Root1       - Temporary value of the Root
|    Root2       - Temporary value of the Root
|    Root3       - Temporary value of the Root
|    s           - alternate variable
|    h           - Temporary value
|    hSqr        - h squared
|    hCube       - h Cubed
|    P           - Term in auxillary equation
|    Q           - Term in auxillary equation
|    R           - Term in auxillary equation
|    Delta       - Discriminator for use with Cardans formula
|    E0          - ANGLE holder for trigonometric solution
|    Phi         - ANGLE used in trigonometric solution
|    CosPhi      - Cosine of Phi
|    SinPhi      - Sine of Phi
|    RPrime      - Values of Roots before final work
|    Temp        - Temporary variable in finding MAX RPrime
|    Eta         - Constant coefficient in QUADRATIC solutions
|    Beta        - Constant coefficient in QUADRATIC solutions
|
|  Coupling      :
|    ATAN2         Arctangent including check for 180-360 deg
|    POWER         Raise a base to an exponent
|    QUADRATIC     Find roots of a QUADRATIC polynomial
|
|  References    :
|    Vallado       2001,
|
 -----------------------------------------------------------------------------*/
void Quartic
    (
      Real  a,   Real  b,   Real  c,   Real  d, Real e,
      Real& R1r, Real& R1i, Real& R2r, Real& R2i, 
      Real& R3r, Real& R3i, Real& R4r, Real& R4i
    )
{
  const Real Rad      = 180.0 / PI;
  const Real OneThird = 1.0 / 3.0;
  const Real Small    = 0.000001;

  Real Temp1, Temp2, Root1, Root2, Root3, s, h, P, Q, R, Delta, E0,
       CosPhi, SinPhi, Phi, RPrime, hSqr, hCube, Eta, Beta, Temp;

  /* ----------------------   Initialize   ----------------------- */
  R1r   = 0.0;
  R1i   = 0.0;
  R2r   = 0.0;
  R2i   = 0.0;
  R3r   = 0.0;
  R3i   = 0.0;
  R4r   = 0.0;
  R4i   = 0.0;
  Root1 = 0.0;
  Root2 = 0.0;
  Root3 = 0.0;

  /* ------------ Force coefficients into std form --------------- */
  b = b / a;
  c = c / a;
  d = d / a;
  e = e / a;

  h     = -b / 4;
  hSqr  = h * h;
  hCube = hSqr * h;

  P =                              6.0 * hSqr   + 3.0 * b * h + c;
  Q =            4.0 * hCube + 3.0 * b * hSqr + 2.0 * c * h + d;
  R = h * hCube +  b * hCube +       c * hSqr   +     d * h + e;

  a =  OneThird * (-P * P -12.0 * R);
  b =  (1.0 / 27.0) * (-2.0 * P * P * P + 72.0 * P * R - 27.0 * Q * Q);
  s = -2.0 * OneThird * P;

  Delta = (a * a * a / 27.0) + (b * b * 0.25);

  if (fabs(Q) > Small)
  {
    /* ------------------- Use Cardans formula --------------------- */
    if (Delta > Small)
    {
      Temp1 = (-b* 0.5) + sqrt(Delta);
      Temp2 = (-b* 0.5) - sqrt(Delta);
      if (fabs(Temp1) > Small)
        Temp1 = Power(Temp1, OneThird);
      if (fabs(Temp2) > Small)
        Temp2 = Power(Temp2, OneThird);
      Root1 = Temp1 + Temp2;
      Root2 = -0.5*(Temp1 + Temp2);
      Root3 = -0.5*(Temp1 + Temp2);
      R2i   = -0.5 * sqrt(3.0) * (Temp1 - Temp2);
      R3i   = -R2i;
    }
    else
    {
      /* ---------------- Evaluate zero point ------------------------ */
      if (fabs(Delta) < Small)
      {
        if (fabs(b) > Small)
        {
          Root1 = -2.0 * Power(0.25 * b, OneThird);
          Root2 = Power(0.25 * b, OneThird);
          Root3 = Root2;
        }
      }
      else
      {
        /* ---------------- Use trigonometric identities --------------- */
        E0     = 2.0 * sqrt(-a * OneThird);
        CosPhi = (-b / (2.0 * sqrt(-a * a * a / 27.0)));
        SinPhi = sqrt(1.0 - CosPhi * CosPhi);
        Phi    = Atan2(SinPhi, CosPhi);
        Root1  = E0 * cos(Phi * OneThird);
        Root2  = E0 * cos(Phi * OneThird + 120.0 / Rad);
        Root3  = E0 * cos(Phi * OneThird + 240.0 / Rad);
      }
    }

    /* ---------------- Find largest value of Root ----------------- */
    RPrime = Root1 + s;
    if ((RPrime < Root2 + s) && (fabs(R2i) < 0.0001))
      RPrime = Root2 + s;
    if ((RPrime < Root3 + s) && (fabs(R3i) < 0.0001))
      RPrime = Root3 + s;

    /* ------ Evaluate coefficients of two resulting Quadratics ---- */
    if (RPrime > Small)
    {
      Eta  = 0.5 * (P + RPrime - Q / sqrt(RPrime));
      Beta = 0.5 * (P + RPrime + Q / sqrt(RPrime));
    }
    else
    {
      Eta  = 0.5 * P;
      Beta = 0.5 * P;
      if (RPrime < 0.0)
        RPrime = -RPrime;
    }

    Quadratic(1.0,  sqrt(RPrime), Eta,     R1r, R1i, R2r, R2i);
    Quadratic(1.0, -sqrt(RPrime), Beta,    R3r, R3i, R4r, R4i);
  }
  else
  {
    /* ------- Case where solution reduces to a QUADRATIC ---------- */
    Quadratic(1.0, P, R,   R1r, R1i, R3r, R3i);
    R = sqrt(R1r * R1r + R1i * R1i);
    if (fabs(R) > Small)
      Phi = Atan2(R1i / R, R1r / R);
    else
      Phi = 0.0;
    R1r = sqrt(R) * cos(Phi * 0.5);
    R1i = sqrt(R) * sin(Phi * 0.5);
    if (fabs(R1i) > Small)
      R2r = R1r;
    else
      R2r = -R1r;
    R2i = -R1i;

    R  = sqrt( R3r * R3r + R3i * R3i );
    if (fabs(R) > Small)
      Phi = Atan2(R3i / R, R3r / R );
    else
      Phi = 0.0;
    R3r = sqrt(R) * cos(Phi * 0.5);
    R3i = sqrt(R) * sin(Phi * 0.5);
    if (fabs(R3i) > Small)
      R4r = R3r;
    else
      R4r = -R3r;
    R4i = -R3i;
  }
  R1r = R1r + h;
  R2r = R2r + h;
  R3r = R3r + h;
  R4r = R4r + h;
}

Real Round(Real n)
{
  if (n < 0)
    return ((Real)((LINT)(n - 0.5)));
  else
    return ((Real)((LINT)(n + 0.5)));
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION SEC
|
|  This FUNCTION finds the secant of an ANGLE in radians.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    XVal        - ANGLE to take secant of                      rad
|
|  OutPuts       :
|    SEC         - Result
|
|  Locals        :
|    Temp        - Temporary Real variable
|
 -----------------------------------------------------------------------------*/
Real Sec(Real v)
{
  Real t;

  t = cos(v);
  if (fabs(t) < 0.000001)
    return 999999.9;
  else
    return (1.0 / t);
}

LINT SignFix(LINT n)
{
  if (ODD(n))
    return -1;
  else
    return 1;
}

/*------------------------------------------------------------------------------
|
|                           FUNCTION SGN
|
|  This FUNCTION determines the sign of a number.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    XVal        - Value to determine sign of
|
|  OutPuts       :
|    SGN         - Result                         +1 or -1
|
|  Locals        :
|    None.
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
Real Sgn(Real v)
{
  if (v > 0.0)
    return (1.0);
  else
    return (-1.0);
}

/*****************************************************************************
 *                    General Utilities
*****************************************************************************/
void FilePrint(Vector& v, char *title, Integer precision, FILE *fp)
{
  Integer dec = precision;

  if (precision > 11)
    dec = 10;
  fprintf(fp, "%s\n", title);
  for (UINT i = 1; i <= v.Dim(); i++)
    fprintf(fp, " %0.*f ", dec, v.Get(i));
  fprintf(fp, "\n");
}

void Print(Vector& v, char* title, Integer precision)
{
  Integer dec = precision;

  if (precision > 11)
    dec = 10;
  printf("%s\n", title);
  for (UINT i = 1; i <= v.Dim(); i++)
    printf(" %0.*f ", dec, v.Get(i));
  printf("\n");
}

void FilePrint(Matrix& m, char *title, Integer precision, FILE *fp)
{
  Integer dec = precision;

  if (precision > 11)
    dec = 10;
  fprintf(fp, "%s\n", title);
  for (UINT r = 1; r <= m.DimR(); r++)
  {
    for (UINT c = 1; c <= m.DimC(); c++)
      fprintf(fp, " %0.*f ", dec, m.Get(r, c));
    fprintf(fp, "\n");
  }
}

void Print(Matrix& m, char *title, Integer precision)
{
  Integer dec = precision;

  if (precision > 11)
    dec = 10;
  printf("%s\n", title);
  for (UINT r = 1; r <= m.DimR(); r++)
  {
    for (UINT c = 1; c <= m.DimC(); c++)
      printf(" %0.*f ", dec, m.Get(r, c));
    printf("\n");
  }
}
