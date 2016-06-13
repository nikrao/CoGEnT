#include "myinclude.h"

// Parameters and macros just for this file
#define RANDOM 0
// #define DATABLOCKIND(DATA,BLOCK,NUMBLOCK) ((DATA+1)-1)*NUMBLOCK + BLOCK - 1

ptrdiff_t main(ptrdiff_t argc, char *argv[])
{
  ptrdiff_t h, k, ret, n, nr, *maxranks, *ranks;
  double *R, *lambda, pieces[8];
  FILE *fid;

  ptrdiff_t m, numblk, *blksz;
  char *blktype;
  double *b, *CAent;
  ptrdiff_t *CArow, *CAcol, *CAinfo_entptr, *CAinfo_rowcolptr;
  char *CAinfo_type, *CAinfo_storage;

  ptrdiff_t numbfgsvecs;
  double rho_f;
  double rho_c;
  double sigmafac;
  ptrdiff_t rankreduce;
  double gaptol;
  ptrdiff_t checkbd;
  ptrdiff_t typebd;
  ptrdiff_t dthresh_dim;
  double dthresh_dens;
  ptrdiff_t timelim;
  double rankredtol;
  ptrdiff_t inputtype;
  ptrdiff_t printlevel;

  char *filein;
  char *fileout;

  // Begin: Perform some simple checks
  if(argc != 2 && argc != 3 && argc != 4 && argc != 5) {
    printf("Usage #1: %s <input_file> [params_file] [soln_in] [soln_out]\n", argv[0]);
    printf("Usage #2: %s gen_params\n", argv[0]);
    exit(0);
  }
  // End: Perform some simple checks

  if(argc == 2 && strcmp("gen_params", argv[1]) == 0) {

    generate_params();

    return 0;

  }

  if(argc == 2) {

    ret = getparams(NULL, &inputtype, &rho_f, &rho_c, &sigmafac,
                    &rankreduce, &timelim, &printlevel, &dthresh_dim,
                    &dthresh_dens, &numbfgsvecs, &rankredtol, &gaptol,
                    &checkbd, &typebd);

  }

  if(argc != 2) {

    ret = getparams(argv[2], &inputtype, &rho_f, &rho_c, &sigmafac,
                    &rankreduce, &timelim, &printlevel, &dthresh_dim,
                    &dthresh_dens, &numbfgsvecs, &rankredtol, &gaptol,
                    &checkbd, &typebd);

  }

  if(ret == 0) printf("Warning (main): Some problems getting parameters.\n");
  if(ret == -1) {
    printf("Error (main): Problem getting parameters.\n");
    exit(0);
  }


       if(argc == 4) { filein = argv[3]; fileout = NULL;    }
  else if(argc == 5) { filein = argv[3]; fileout = argv[4]; }
  else               { filein = NULL;    fileout = NULL;    }

  if(inputtype == 1) {

    readdata_sdpa(argv[1], &m, &numblk, &blksz, &blktype, &b, &CAent, &CArow, &CAcol,
                  &CAinfo_entptr, &CAinfo_rowcolptr, &CAinfo_type, &CAinfo_storage);

  }
  if(inputtype == 2) {

    readdata_sdplr(argv[1], &m, &numblk, &blksz, &blktype, &b, &CAent, &CArow, &CAcol,
                  &CAinfo_entptr, &CAinfo_rowcolptr, &CAinfo_type, &CAinfo_storage);

  }
  if(inputtype == 1000) {

     readdata_raw(argv[1], &m, &numblk, &blksz, &blktype, &b, &CAent, &CArow, &CAcol,
                  &CAinfo_entptr, &CAinfo_rowcolptr, &CAinfo_type, &CAinfo_storage);

  }

  // Quick error check

  for(k = 0; k < numblk; k++)
    if(blksz[k] == 0) {
      printf("Error (main): Block %d has size 0.\n", k);
      exit(0);
    }


  // Allocate space and setup info

  // First get dimensions and ranks
  MYCALLOC(maxranks, ptrdiff_t, numblk);
  MYCALLOC(ranks, ptrdiff_t, numblk);
  getstorage(m, numblk, blksz, blktype, CAinfo_entptr, &n, &nr, maxranks);
  for(k = 0; k < numblk; k++) ranks[k] = maxranks[k];

  // Allocate space
  MYCALLOC(R, double, nr);
  MYCALLOC(lambda, double, m);

  // Setup initial point; either default or from file
  if(filein != NULL && (fid = fopen(filein, "r")) != NULL) {
    readin(m,numblk,blksz,blktype,R,lambda,maxranks,ranks,pieces,fid);
    fclose(fid);
  }
  else {
    if(RANDOM) srand( (size_t)time( NULL ) );
    else       srand(925);
    for(h = 0; h < nr; h++) {
      R[h]  = (double)rand()/RAND_MAX;
      R[h] -= (double)rand()/RAND_MAX;
    }
    pieces[0] = (double)0;
    pieces[1] = (double)0;
    pieces[2] = (double)0;
    pieces[3] = (double)0;
    pieces[4] = (double)0;
    pieces[5] = (double)0.0;
    pieces[6] = (double)1.0/n;
    pieces[7] = (double)1.0;
  }

//    writedata_sdpa("temp.dat-s", m, numblk, blksz, blktype, b, CAent, CArow, CAcol, CAinfo_entptr, CAinfo_type);
//   exit(0);
  
  sdplrlib(m, numblk, blksz, blktype, b, CAent, CArow, CAcol,
           CAinfo_entptr, CAinfo_type,
           numbfgsvecs, rho_f, rho_c, sigmafac, rankreduce, gaptol,
           checkbd, typebd, dthresh_dim, dthresh_dens, timelim,
           rankredtol, printlevel, R-1, lambda-1, maxranks, ranks,
           pieces);

  if(fileout != NULL && (fid = fopen(fileout, "w")) != NULL) {
    writeout(m,numblk,blksz,blktype,R,lambda,maxranks,ranks,pieces,fid);
    fclose(fid);
  }

  MYFREE(R);
  MYFREE(lambda);
  MYFREE(maxranks);
  MYFREE(ranks);

  if(inputtype == 1 || inputtype == 2) {
    MYFREE(blksz);
    MYFREE(blktype);
    MYFREE(b);
    MYFREE(CAent);
    MYFREE(CArow);
    MYFREE(CAcol);
    MYFREE(CAinfo_entptr);
    MYFREE(CAinfo_rowcolptr);
    MYFREE(CAinfo_type);
    MYFREE(CAinfo_storage);
  }

  return 0;
}


ptrdiff_t getstorage(ptrdiff_t m, ptrdiff_t numblk, ptrdiff_t* blksz, char* blktype, 
               ptrdiff_t* CAinfo_entptr, ptrdiff_t* passedn, ptrdiff_t* passednr,
               ptrdiff_t* maxranks)
{
  ptrdiff_t i, k, *cons, n, nr, ind;

  MYCALLOC(cons, ptrdiff_t, m+1);

  n = 0;
  nr = 0;

  for(k = 1; k <= numblk; k++) {

    if(blktype[k-1] == 's') {

      for(i = 1; i <= m; i++) {

        ind = DATABLOCKIND(i,k,numblk);
        if(CAinfo_entptr[ind+1] > CAinfo_entptr[ind]) cons[i] = 1;
        else cons[i] = 0;

      }

      cons[0] = 0;
      for(i = 1; i <= m; i++) cons[0] += cons[i];
      // Note: this formula MUST match exactly the one found in initialize.c
      maxranks[k-1] = mymin((ptrdiff_t)sqrt(2*cons[0]) + 1, blksz[k-1]);
      nr += blksz[k-1]*maxranks[k-1];
      n += blksz[k-1];

    }
    else if(blktype[k-1] == 'd') {
      maxranks[k-1] = 1;
      nr += blksz[k-1];
      n += blksz[k-1];
    }
  }

  *passedn  = n;
  *passednr = nr;

  MYFREE(cons);

  return 0;

}

ptrdiff_t readin(ptrdiff_t m, ptrdiff_t numblk, ptrdiff_t* blksz, char* blktype, double* R, double* lambda, ptrdiff_t* maxranks, ptrdiff_t* ranks, double* pieces, FILE* fid)
{
  ptrdiff_t h, k, base, ti1, ti2, ti3;
  char tc1;

  fscanf(fid, "dual variable %d\n", &ti1);

  if(ti1 != m) {
    printf("Error (readin): Input solution and problem file do not match.\n");
    exit(0);
  }

  for(h = 0; h < m; h++) fscanf(fid, "%lf\n", &(lambda[h]));

  base = 0;
  for(k = 0; k < numblk; k++) {

    fscanf(fid, "primal variable %d %c %d %d %d\n", &ti1, &tc1, &ti2, &ti3, &(ranks[k])); ti1--;

    if(ti1 != k || tc1 != blktype[k] || ti2 != blksz[k] || ti3 != maxranks[k]) {
      printf("Error (readin): Input solution and problem file do not match.\n");
      exit(0);
    }

    for(h = 0; h < blksz[k]*ranks[k]; h++) fscanf(fid, "%lf\n", &(R[base + h]));
    base += blksz[k]*ranks[k];
  }

  fscanf(fid, "special majiter ");
  fscanf(fid, "%lf\n", &(pieces[0]));
  fscanf(fid, "special iter ");
  fscanf(fid, "%lf\n", &(pieces[1]));
  fscanf(fid, "special lambdaupdate ");
  fscanf(fid, "%lf\n", &(pieces[2]));
  fscanf(fid, "special CG ");
  fscanf(fid, "%lf\n", &(pieces[3]));
  fscanf(fid, "special curr_CG ");
  fscanf(fid, "%lf\n", &(pieces[4]));
  fscanf(fid, "special totaltime ");
  fscanf(fid, "%lf\n", &(pieces[5]));
  fscanf(fid, "special sigma ");
  fscanf(fid, "%lf\n", &(pieces[6]));
  fscanf(fid, "special scale ");
  fscanf(fid, "%lf\n", &(pieces[7]));

  return 0;
}

ptrdiff_t writeout(ptrdiff_t m, ptrdiff_t numblk, ptrdiff_t* blksz, char* blktype, double* R, double* lambda, ptrdiff_t* maxranks, ptrdiff_t* ranks, double* pieces, FILE* fid)
{
  ptrdiff_t h, k, base;

  fprintf(fid, "dual variable %d\n", m);
  for(h = 0; h < m; h++) fprintf(fid, "%.16e\n", lambda[h]*pieces[7]);
  base = 0;
  for(k = 0; k < numblk; k++) {
    fprintf(fid, "primal variable %d %c %d %d %d\n", k+1, blktype[k], blksz[k], maxranks[k], ranks[k]);
    for(h = 0; h < blksz[k]*ranks[k]; h++)
      fprintf(fid, "%.16e\n", R[base + h]);
    base += blksz[k]*ranks[k];
  }
  fprintf(fid, "special majiter ");
  fprintf(fid, "%d\n", (ptrdiff_t)pieces[0]);
  fprintf(fid, "special iter ");
  fprintf(fid, "%d\n", (ptrdiff_t)pieces[1]);
  fprintf(fid, "special lambdaupdate ");
  fprintf(fid, "%d\n", (ptrdiff_t)pieces[2]);
  fprintf(fid, "special CG ");
  fprintf(fid, "%d\n", (ptrdiff_t)pieces[3]);
  fprintf(fid, "special curr_CG ");
  fprintf(fid, "%d\n", (ptrdiff_t)pieces[4]);
  fprintf(fid, "special totaltime ");
  fprintf(fid, "%.16e\n", (double)pieces[5]);
  fprintf(fid, "special sigma ");
  fprintf(fid, "%.16e\n", (double)pieces[6]);
  fprintf(fid, "special scale ");
  fprintf(fid, "%.16e\n", (double)pieces[7]);

  return 0;
}


