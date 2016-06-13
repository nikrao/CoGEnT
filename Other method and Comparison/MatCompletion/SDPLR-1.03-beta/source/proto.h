/* copystructures.c */
ptrdiff_t copystructures(problemdata *data, ptrdiff_t m, ptrdiff_t numblk, ptrdiff_t *blksz, char *blktype, double *b, double *CAent, ptrdiff_t *CArow, ptrdiff_t *CAcol, ptrdiff_t *CAinfo_entptr, char *CAinfo_type);
/* dataoper.c */
double function(problemdata *data, double *R);
ptrdiff_t Aoper(problemdata *d, double *U, double *V, double *UVt, ptrdiff_t same, ptrdiff_t obj, double *results);
ptrdiff_t Aoper_formUVt(problemdata *data, double *passedUVt, double *U, double *V, ptrdiff_t same);
ptrdiff_t gradient(problemdata *data, double *R);
ptrdiff_t StimesR(problemdata *data, double *S, double *y, double *R, double *result);
ptrdiff_t Stimesmat(problemdata *data, double *S, double *y, double *vec, double *result, ptrdiff_t n, ptrdiff_t m, ptrdiff_t k);
ptrdiff_t AToper(problemdata *data, double *y, double *S, ptrdiff_t obj);
double C_normdatamat(problemdata *data);
double normdatamat(problemdata *data, ptrdiff_t matnum);
ptrdiff_t essential_calcs(problemdata *data, double *R, double normC, double normb, double *val, double *rho_c_val, double *rho_f_val);
/* eigen.c */
ptrdiff_t my_arpack(ptrdiff_t (*matvec)(double *, double *), ptrdiff_t n, ptrdiff_t nev, ptrdiff_t ncv, char *which, ptrdiff_t maxitr, ptrdiff_t printlevel, double *evals, double *evecs, ptrdiff_t *nconv, ptrdiff_t *nummatvec);
ptrdiff_t simple_Stimesvec_block(double *out, double *in);
ptrdiff_t Smineval(problemdata *d, double *eval);
ptrdiff_t Sblockmineval(problemdata *d, double *blkevals);
/* initialize.c */
ptrdiff_t initialize(problemdata *data, ptrdiff_t *maxranks);
ptrdiff_t deinitialize(problemdata *data);
ptrdiff_t locatetype(problemdata *data, ptrdiff_t blk, char type, ptrdiff_t **passed_ind, ptrdiff_t *passed_nind);
/* lbfgs.c */
ptrdiff_t dirlbfgs(problemdata *data, lbfgsvec *vecs, double *dir, double *grad, ptrdiff_t oldest, ptrdiff_t numlbfgsvecs, ptrdiff_t scale);
ptrdiff_t updatelbfgs1(problemdata *data, lbfgsvec *vecs, double *grad, ptrdiff_t oldest);
ptrdiff_t updatelbfgs2(problemdata *data, lbfgsvec *vecs, double *dir, double *grad, double stepsize, ptrdiff_t *oldest, ptrdiff_t numlbfgsvecs);
/* linesearch.c */
double linesearch(problemdata *data, double *R, double *D, double max, double *funcval, ptrdiff_t update);
/* main.c */
ptrdiff_t main(ptrdiff_t argc, char *argv[]);
ptrdiff_t getstorage(ptrdiff_t m, ptrdiff_t numblk, ptrdiff_t *blksz, char *blktype, ptrdiff_t *CAinfo_entptr, ptrdiff_t *passedn, ptrdiff_t *passednr, ptrdiff_t *maxranks);
ptrdiff_t readin(ptrdiff_t m, ptrdiff_t numblk, ptrdiff_t *blksz, char *blktype, double *R, double *lambda, ptrdiff_t *maxranks, ptrdiff_t *ranks, double *pieces, FILE *fid);
ptrdiff_t writeout(ptrdiff_t m, ptrdiff_t numblk, ptrdiff_t *blksz, char *blktype, double *R, double *lambda, ptrdiff_t *maxranks, ptrdiff_t *ranks, double *pieces, FILE *fid);
/* mexsdplr.c */
ptrdiff_t classifyApq(ptrdiff_t *blksz, char *blktype, ptrdiff_t *inblk, ptrdiff_t *inblkptr, ptrdiff_t p, ptrdiff_t q, ptrdiff_t A_tr, ptrdiff_t *cons, ptrdiff_t *blk, ptrdiff_t *i, ptrdiff_t *j);
/* misc.c */
ptrdiff_t print_notes(void);
ptrdiff_t printparams(problemdata *data);
ptrdiff_t print_dimacs_errors(problemdata *data, double *R);
ptrdiff_t printheading(ptrdiff_t start);
/* params.c */
ptrdiff_t getparams(char *paramfile, ptrdiff_t *inputtype, double *rho_f, double *rho_c, double *sigmafac, ptrdiff_t *rankreduce, ptrdiff_t *timelim, ptrdiff_t *printlevel, ptrdiff_t *dthresh_dim, double *dthresh_dens, ptrdiff_t *numbfgsvecs, double *rankredtol, double *gaptol, ptrdiff_t *checkbd, ptrdiff_t *typebd);
ptrdiff_t getparams_maxlinelength(FILE *datafile);
ptrdiff_t getparams_getline(FILE *datafile, char *buffer, ptrdiff_t bufsiz);
ptrdiff_t getparams_tolower(char *buff, ptrdiff_t buffsz);
ptrdiff_t getparams_isnumber(char *str);
ptrdiff_t generate_params(void);
ptrdiff_t generate_params_info(ptrdiff_t num);
/* rankreduce.c */
ptrdiff_t dorankreduce(problemdata *d, double *R);
/* readdata.c */
ptrdiff_t readdata_sdpa(char *datafilename, ptrdiff_t *passed_m, ptrdiff_t *passed_numblk, ptrdiff_t **passed_blksz, char **passed_blktype, double **passed_b, double **passed_CAent, ptrdiff_t **passed_CArow, ptrdiff_t **passed_CAcol, ptrdiff_t **passed_CAinfo_entptr, ptrdiff_t **passed_CAinfo_rowcolptr, char **passed_CAinfo_type, char **passed_CAinfo_storage);
ptrdiff_t quicksort5(ptrdiff_t *A1, ptrdiff_t *A2, ptrdiff_t *A3, ptrdiff_t *A4, double *A5, ptrdiff_t p, ptrdiff_t r);
ptrdiff_t partition5(ptrdiff_t *A1, ptrdiff_t *A2, ptrdiff_t *A3, ptrdiff_t *A4, double *A5, ptrdiff_t p, ptrdiff_t r);
void skip_to_end_of_line(FILE *datafile);
ptrdiff_t get_line(FILE *datafile, char *buffer, ptrdiff_t bufsiz);
ptrdiff_t max_line_length(FILE *datafile);
ptrdiff_t readdata_sdplr(char *datafilename, ptrdiff_t *passed_m, ptrdiff_t *passed_numblk, ptrdiff_t **passed_blksz, char **passed_blktype, double **passed_b, double **passed_CAent, ptrdiff_t **passed_CArow, ptrdiff_t **passed_CAcol, ptrdiff_t **passed_CAinfo_entptr, ptrdiff_t **passed_CAinfo_rowcolptr, char **passed_CAinfo_type, char **passed_CAinfo_storage);
ptrdiff_t readdata_raw(char *datafilename, ptrdiff_t *passed_m, ptrdiff_t *passed_numblk, ptrdiff_t **passed_blksz, char **passed_blktype, double **passed_b, double **passed_CAent, ptrdiff_t **passed_CArow, ptrdiff_t **passed_CAcol, ptrdiff_t **passed_CAinfo_entptr, ptrdiff_t **passed_CAinfo_rowcolptr, char **passed_CAinfo_type, char **passed_CAinfo_storage);
ptrdiff_t writedata_raw(char *datafilename, ptrdiff_t m, ptrdiff_t numblk, ptrdiff_t *blksz, char *blktype, double *b, double *CAent, ptrdiff_t *CArow, ptrdiff_t *CAcol, ptrdiff_t *CAinfo_entptr, ptrdiff_t *CAinfo_rowcolptr, char *CAinfo_type, char *CAinfo_storage);
ptrdiff_t writedata_sdplr(char *datafilename, ptrdiff_t m, ptrdiff_t numblk, ptrdiff_t *blksz, char *blktype, double *b, double *CAent, ptrdiff_t *CArow, ptrdiff_t *CAcol, ptrdiff_t *CAinfo_entptr, char *CAinfo_type);
ptrdiff_t writedata_sdpa(char *datafilename, ptrdiff_t m, ptrdiff_t numblk, ptrdiff_t *blksz, char *blktype, double *b, double *CAent, ptrdiff_t *CArow, ptrdiff_t *CAcol, ptrdiff_t *CAinfo_entptr, char *CAinfo_type);
/* sdplrlib.c */
ptrdiff_t sdplrlib(ptrdiff_t m, ptrdiff_t numblk, ptrdiff_t *blksz, char *blktype, double *b, double *CAent, ptrdiff_t *CArow, ptrdiff_t *CAcol, ptrdiff_t *CAinfo_entptr, char *CAinfo_type, ptrdiff_t numbfgsvecs, double rho_f, double rho_c, double sigmafac, ptrdiff_t rankreduce, double gaptol, ptrdiff_t checkbd, ptrdiff_t typebd, ptrdiff_t dthresh_dim, double dthresh_dens, ptrdiff_t timelim, double rankredtol, ptrdiff_t printlevel, double *R, double *lambda, ptrdiff_t *maxranks, ptrdiff_t *ranks, double *pieces);
#ifdef Kaushik_addition
void myprint(ptrdiff_t majiter, ptrdiff_t iter, double val, double rho_f_val, 
ptrdiff_t CG, double totaltime, int rank, double sigma);
#else
void myprint(ptrdiff_t majiter, ptrdiff_t iter, double val, double rho_f_val, 
    ptrdiff_t CG, double totaltime);
#endif
ptrdiff_t do_scaling(problemdata *data, double value, double *norm);
/* timefuncs.c */
#ifndef __WIN32
double current_time(double timeorig);
#else
double current_time(clock_t timeorig);
#endif
/* util.c */
ptrdiff_t copyscaledvectovec(double *dy, double da, double *dx, ptrdiff_t n);
ptrdiff_t move_in_dir(double *dy, double *dx, double da, double *dir, ptrdiff_t n);
ptrdiff_t mydaxpy(ptrdiff_t n, double da, double *dx, ptrdiff_t incx, double *dy, ptrdiff_t incy);
ptrdiff_t mydcopy(ptrdiff_t n, double *dx, ptrdiff_t incx, double *dy, ptrdiff_t incy);
double myddot(ptrdiff_t n, double *dx, ptrdiff_t incx, double *dy, ptrdiff_t incy);
double mydnrm2(ptrdiff_t n, double *dx, ptrdiff_t incx);
ptrdiff_t mydscal(ptrdiff_t n, double da, double *dx, ptrdiff_t incx);
ptrdiff_t createlowrankmat(lowrankmat **passedR, ptrdiff_t ncol, ptrdiff_t nrow);
ptrdiff_t destroylowrankmat(lowrankmat *R);
ptrdiff_t createsparsesymmmat(sparsesymmmat **passedS, ptrdiff_t nnz);
ptrdiff_t destroysparsesymmmat(sparsesymmmat *S);
ptrdiff_t creatediagmat(diagmat **passedD, ptrdiff_t nnz);
ptrdiff_t destroydiagmat(diagmat *D);
ptrdiff_t createdatamat(datamat **passedA, char type, ptrdiff_t ncol_or_nnz, ptrdiff_t dim, char *label);
ptrdiff_t destroydatamat(datamat *A);
