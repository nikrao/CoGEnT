//! \brief Core iteration information for implementation of L-BFGS.
//!
//! SDPLR uses the limited memory BFGS (L-BFGS) method as its core
//! unconstrained minimization procedure. This method is based on
//! storing several types of vectors and scalars for each iteration of
//! the algorithm. These vectors and scalars are more or less linear
//! combinations of points and previous points. This structure is the
//! basic collection of information for one iteration. In SDPLR, an
//! array of structures of this type will be maintained and updated.
//! This size of this array is \ref NUMBFGSVECS, which is the so-called
//! "memory" of the L-BFGS procedure and is a user-specified integer
//! parameter.

typedef struct {

  //! \brief Difference between iterates
  //!
  //! This stores the difference between successive iterates, i.e.,
  //! \f$ x^{k+1} - x^k \f$. (Here, \f$ x \f$ is generic.)

  double *s; 

  //! \brief Difference between gradients
  //!
  //! This stores the difference between successive gradients, i.e.,
  //! \f$ \nabla f(x^{k+1}) - \nabla f(x^k) \f$. (Here, \f$ x \f$ and
  //! \f$ f \f$ are generic.)
 
  double *y;

  //! \brief Reciprocal of inner product between \ref s and \ref y
  //!
  //! This stores \f$ (y^T s)^{-1} \f$. It is significant only as
  //! a temporary storage places.

  double rho;

  //! \brief Temporary storage.
  //!
  //! Significant only as a temporary storage place related to \ref s,
  //! \ref y, and \ref rho.
  
  double a;

} lbfgsvec;

typedef struct {
  double*  d;
  double*  ent;
  ptrdiff_t      nrow;
  ptrdiff_t      ncol;
} lowrankmat;

typedef struct {
  ptrdiff_t*    row;
  ptrdiff_t*    col;
  ptrdiff_t     nnz;
  double* ent;
  ptrdiff_t*    XS_in;
} sparsesymmmat;

typedef struct {
  ptrdiff_t*    ind;
  ptrdiff_t     nnz;
  double* ent;
  ptrdiff_t*    XS_in;
} diagmat;

typedef struct {
  lowrankmat*    lr;
  sparsesymmmat* sp;
  diagmat*       diag;
  char           type;
  double*        multval;
  char*          label;
} datamat;

typedef struct {

  // user options
  double     rho_f;
  double     rho_c;
  double     sigmafac;
  ptrdiff_t        rankreduce;
  ptrdiff_t        timelim;
  ptrdiff_t        printlevel;
  ptrdiff_t        dthresh_dim;
  double     dthresh_dens;
  ptrdiff_t        numbfgsvecs;
  double     rankredtol;
  double     gaptol;
  ptrdiff_t        checkbd;
  ptrdiff_t        typebd;

  // very basic data (both SDPA and lowrank formats)
  ptrdiff_t        m;
  ptrdiff_t        numblk;
  ptrdiff_t*       blksz;
  char*      blktype;
  datamat*** A;
  double*    b;
  datamat**  C;

  // auxiliary data, applies to lowrank format only
  double     evaladj;

  //
  // algorithm data, doesn't change once it is setup
  //
  // global dimension
  ptrdiff_t        n;
  ptrdiff_t        nr;
  // number of limited memory vecs to store
  ptrdiff_t        M;
  // pointers to different types of matrices
  ptrdiff_t**      lrind;
  ptrdiff_t*       nlrind;
  ptrdiff_t**      usind;
  ptrdiff_t*       nusind;
  // rank information
  ptrdiff_t*       rank;
  ptrdiff_t*       maxrank;
  // norm of C for internal stopping criterion
  double     normC;


  //
  // algorithm data, changes per iteration
  //
  // basic structures
  double*      lambda;
  double       sigma;
  double*      vio;
  double*      G;

  // timing data structures
  double totaltime;
  clock_t start;
  clock_t finish;

  // What are these?
  double   *S;
  double   *X;
  double   *y;
  ptrdiff_t      *XS_blkptr;
  char     *XS_blksto;
  ptrdiff_t     **XS_colptr;
  ptrdiff_t     **XS_rowind;
  ptrdiff_t      *AA_rowptr;
  ptrdiff_t      *AA_colind;
  double   *AA_colval_one;
  double   *AA_colval_two;
  ptrdiff_t      *lr_mat;
  ptrdiff_t      *lr_blk;
  ptrdiff_t       lr_num;

} problemdata;



