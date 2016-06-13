#include "myinclude.h"

// #define mymax(A, B) ((A) > (B) ? (A) : (B))
// #define mymin(A, B) ((A) < (B) ? (A) : (B))


ptrdiff_t copyscaledvectovec(double* dy, double da, double* dx, ptrdiff_t n)
{
  EASYDCOPY(n, dx, dy);
  EASYDSCAL(n, da, dy);
  return 0;
}

ptrdiff_t move_in_dir(double* dy, double* dx, double da, double* dir, ptrdiff_t n)
{
  if(dy == dx) EASYDAXPY(n, da, dir, dy);
  else if(dy == dir) {
    EASYDSCAL(n, da, dy);
    EASYDAXPY(n, 1.0, dx, dy);
  }
  else {
    EASYDCOPY(n, dx, dy);
    EASYDAXPY(n, da, dir, dy);
  }
  return 1;
}

ptrdiff_t mydaxpy(ptrdiff_t n, double da, double* dx, ptrdiff_t incx, double* dy, ptrdiff_t incy)
{
  return daxpy_(&n,&da,dx,&incx,dy,&incy);
}

ptrdiff_t mydcopy(ptrdiff_t n, double* dx, ptrdiff_t incx, double* dy, ptrdiff_t incy)
{
  return dcopy_(&n,dx,&incx,dy,&incy);
}

double myddot(ptrdiff_t n, double* dx, ptrdiff_t incx, double* dy, ptrdiff_t incy)
{
  return ddot_(&n, dx, &incx, dy, &incy);
}

double mydnrm2(ptrdiff_t n, double* dx, ptrdiff_t incx)
{
  return dnrm2_(&n, dx, &incx);
}

ptrdiff_t mydscal(ptrdiff_t n, double da, double* dx, ptrdiff_t incx)
{
  return dscal_(&n, &da, dx, &incx);
}


ptrdiff_t createlowrankmat(lowrankmat** passedR, ptrdiff_t ncol, ptrdiff_t nrow)
{
  lowrankmat *R;

  MYCALLOC(R, lowrankmat, 1);
  
  R->ncol = ncol;
  R->nrow = nrow;
  
  MYCALLOC(R->d, double, ncol + 1);
  MYCALLOC(R->ent, double, nrow*ncol + 1);

  *passedR = R;

  return 1;
}

ptrdiff_t destroylowrankmat(lowrankmat* R)
{
  MYFREE(R->d);
  MYFREE(R->ent);
  MYFREE(R);

  return 1;
}

ptrdiff_t createsparsesymmmat(sparsesymmmat** passedS, ptrdiff_t nnz)
{
  sparsesymmmat *S;

  MYCALLOC(S, sparsesymmmat, 1);
  MYCALLOC(S->row, ptrdiff_t, nnz+1);
  MYCALLOC(S->col, ptrdiff_t, nnz+1);
  S->nnz = nnz;
  MYCALLOC(S->ent, double, nnz+1);
  MYCALLOC(S->XS_in, ptrdiff_t, nnz+1);

  *passedS = S;

  return 1;

}

ptrdiff_t destroysparsesymmmat(sparsesymmmat* S)
{
  MYFREE(S->row);
  MYFREE(S->col);
  MYFREE(S->ent);
  MYFREE(S->XS_in);
  MYFREE(S);

  return 1;
}

ptrdiff_t creatediagmat(diagmat** passedD, ptrdiff_t nnz)
{
  diagmat *D;

  MYCALLOC(D, diagmat, 1);
  MYCALLOC(D->ind, ptrdiff_t, nnz+1);
  D->nnz = nnz;
  MYCALLOC(D->ent, double, nnz+1);
  MYCALLOC(D->XS_in, ptrdiff_t, nnz+1);

  *passedD = D;

  return 1;

}

ptrdiff_t destroydiagmat(diagmat* D)
{
  MYFREE(D->ind);
  MYFREE(D->ent);
  MYFREE(D->XS_in);
  MYFREE(D);

  return 1;
}


ptrdiff_t createdatamat(datamat** passedA, char type, ptrdiff_t ncol_or_nnz, ptrdiff_t dim, char* label)
{
  datamat *A;

  MYCALLOC(A, datamat, 1);
  A->type = type;
  MYCALLOC(A->label, char, 30);
  strcpy(A->label, label);

  if(type == 'l')
    createlowrankmat(&(A->lr), ncol_or_nnz, dim);
  
  if(type == 's')
    createsparsesymmmat(&(A->sp), ncol_or_nnz);

  if(type == 'd')
    creatediagmat(&(A->diag), ncol_or_nnz);

  // if type = 'u' then do nothing

  *passedA = A;

  return 1;
}

ptrdiff_t destroydatamat(datamat* A)
{
  if(A->type == 'l')
    destroylowrankmat(A->lr);
  
  if(A->type == 's')
    destroysparsesymmmat(A->sp);

  if(A->type == 'd')
    destroydiagmat(A->diag);

  MYFREE(A->label);
  MYFREE(A);

  return 1;
}

