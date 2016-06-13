#include "myinclude.h"

#ifdef Kaushik_addition

int int_cmp (const int *n1, const int *n2)
{
  return n1[0]-n2[0];
}

int is_nnz(int q_i, int num_ind, int *nnz_row_ind)
{
  int i, result;
  result = 0;
  for (i=1; i <= num_ind; i++){
    if (nnz_row_ind[i] == q_i){
      result = 1;
    }
  }
  return result;
}

int agg_id(int k, int q_i, int q_j, int num_ind, int *nnz_row_ind, int 
*nnz_row_cnt, int *nnz_rows, int **nnz_row_id, int **nnz_row_ct)
{
  int i, id;
  id = 0;
  for (i=1; i <= num_ind; i++){
    if (nnz_row_ind[i] == q_i){
      id = nnz_row_cnt[i];
    }
  }
  if (id == 0){
    printf("%d %d %d %d\n", id, k, q_i, q_j);
  }
  return id;
}
#endif


// ptrdiff_t locatetype(problemdata *data, ptrdiff_t blk, char type, ptrdiff_t **passed_ind, ptrdiff_t *passed_nind);

ptrdiff_t initialize (problemdata * data, ptrdiff_t * maxranks)
{
  ptrdiff_t                 h, i, j, k, ct;
  /* ptrdiff_t  p, q; */
  ptrdiff_t               **AGGMAT;
  ptrdiff_t                *colptr;
  sparsesymmmat      *A;
  /* ptrdiff_t                *xadj, *adjncy, numedges, dim, options[10]; */
  /* ptrdiff_t                *perm, *invp; */
  ptrdiff_t                *spind, nspind;
  /* ptrdiff_t  *temprow, *tempcol; */
  /* char                type; */
  /* lowrankmat         *LR; */
  // Temp variables
  datamat* dA;
  ptrdiff_t col, row, ctr;
  
#ifdef Kaushik_addition
  int **nnz_row_ind, **nnz_row_cnt;
  int *CA_col, *CA_row, *nnz_rows;
  int total_nnz, cnt, num_cols;
#endif

  // The following three sections of code are very similar
  // and could be simplified.

  // Setup pointers to find lowrank data matrices. Used in several routines.
  // Right now, we do not set this up for diagonal blocks.
  MYCALLOC (data->lrind, ptrdiff_t *, data->numblk + 1);
  MYCALLOC (data->nlrind, ptrdiff_t, data->numblk + 1);
  for (k = 1; k <= data->numblk; k++)
    if (data->blktype[k] == 's') 
      locatetype(data, k, 'l', &(data->lrind[k]), &(data->nlrind[k]));

  // BEGIN: Setup new data structures
  // BEGIN: Setup new data structures
  // BEGIN: Setup new data structures
  // BEGIN: Setup new data structures
  
  // Setup structure of X- and S-type matrices.
  
  // For non-diag blocks, this uses somewhat exorbitant allocation of
  // an entire blksz[k] x blksz[k] lower triangular integer matrix,
  // but this is the easiest thing right now.
 
  // Allocate space for block related pointers
  
  MYCALLOC(data->XS_blkptr, ptrdiff_t , data->numblk+2);
  MYCALLOC(data->XS_blksto, char, data->numblk+1);
  MYCALLOC(data->XS_colptr, ptrdiff_t*, data->numblk+1);
  MYCALLOC(data->XS_rowind, ptrdiff_t*, data->numblk+1);

  // First pointer is to position 1

  data->XS_blkptr[1] = 1;

  // Loop over blocks

  for (k = 1; k <= data->numblk; k++) {

    // First find location of sparse data matrices (i.e., not low-rank
    // or user) for this block.
   
    locatetype(data, k, 's', &spind, &nspind);

    // Condition on block type

    if(data->blktype[k] == DIAGBLK) {

      // If block type is diagonal, there is no need for colptr or
      // rowind. Just set blkptr and blksto.

      data->XS_blkptr[k+1] = data->XS_blkptr[k] + data->blksz[k];
      data->XS_blksto[k] = DENSE;

      // Determine absolute position of nonzero ntries of k-th block of
      // Ai in k-th block of S.

      for(i = 0; i <= data->m; i++) {
        if(i == 0) dA = data->C[k];
        else       dA = data->A[i][k];
        if(dA->type == 'd')
          for(j = 1; j <= dA->diag->nnz; j++)
            dA->diag->XS_in[j] = data->XS_blkptr[k] - 1 + dA->diag->ind[j];
            // printf("A(i=%d,k=%d,r=%d,c=%d) --> %d\n", i, k, dA->diag->ind[j], dA->diag->ind[j], data->XS_blkptr[k] - 1 + dA->diag->ind[j]);
      }


    }
    else if(data->blktype[k] == SDPBLK) {

      // If block type is SDP...

      // First thing is to determine whether storage will be sparse
      // or dense. So we figure out positions of aggregate nonzeros
      // (lower triangular part only) and count them.
#ifndef Kaushik_addition
      MYCALLOC(AGGMAT, ptrdiff_t*, data->blksz[k] + 1);
      for(i = 1; i <= data->blksz[k]; i++) MYCALLOC(AGGMAT[i], ptrdiff_t, i + 1);
      for(i = 1; i <= data->blksz[k]; i++) for(j = 1; j <= i; j++) AGGMAT[i][j] = 0;

      for (i = 1; i <= nspind; i++) {
        if(spind[i] == 0) A = data->C          [k]->sp;
        else              A = data->A[spind[i]][k]->sp;
        for(j = 1; j <= A->nnz; j++) AGGMAT[A->row[j]][A->col[j]] = 1;
      }
#else
      //Collecting all the nnz row and col indices in CA_row and CA_col
      total_nnz = 0;
      for (i = 1; i <= nspind; i++) {
        if(spind[i] == 0) A = data->C          [k]->sp;
        else              A = data->A[spind[i]][k]->sp;
        total_nnz = total_nnz + A->nnz; }
      cnt = 1;
      CA_col = (int*) calloc(total_nnz+1, sizeof(int));
      CA_row = (int*) calloc(total_nnz+1, sizeof(int));
      for (i = 1; i <= nspind; i++) {
        if(spind[i] == 0) A = data->C          [k]->sp;
        else              A = data->A[spind[i]][k]->sp;
        for(j = 1; j <= A->nnz; j++) {
          CA_col[cnt] = A->col[j]; CA_row[cnt] = A->row[j]; cnt++;
        }
      }

      num_cols = data->blksz[k];
      //row count for each column
      nnz_rows = (int*) calloc(num_cols+1, sizeof(int));
      for(col = 1; col <= num_cols; col++){
        nnz_rows[col] = 0;
        for(j = 1; j<=total_nnz ; j++) if(col == CA_col[j]) nnz_rows[col]++;
      }
      
      //non-zero row indices and its aggregate count for each column
      nnz_row_ind = (int**) calloc(num_cols+1, sizeof(int*));
      nnz_row_cnt = (int**) calloc(num_cols+1, sizeof(int*));
      for(col = 1; col <= num_cols; col++){
        nnz_row_ind[col] = (int*) calloc(nnz_rows[col]+1, sizeof(int));
        nnz_row_cnt[col] = (int*) calloc(nnz_rows[col]+1, sizeof(int));
      }

      cnt = 1;
      for(col = 1; col <= num_cols; col++){
    
        //total nnz row indices per col
        i = 1;
        for(j = 1; j<=total_nnz; j++){
          if(col == CA_col[j]){
            nnz_row_ind[col][i] = CA_row[j];
            i++;
          }
        }
        
        //Sorting and finding unique row indices
        qsort (nnz_row_ind[col]+1, nnz_rows[col], sizeof(int), int_cmp);
        for(i=1 ; i<=nnz_rows[col]; i++)
        {
          if(i == 1)
          {
            nnz_row_cnt[col][i] = cnt;
            cnt++;
          }else {
            if( nnz_row_ind[col][i] == nnz_row_ind[col][i-1] ){
              for(j = i; j < nnz_rows[col]; j++) 
                nnz_row_ind[col][j] = nnz_row_ind[col][j+1];
              nnz_rows[col]--;
              i--;
            }else {
              nnz_row_cnt[col][i] = cnt;
              cnt++;
            }
          }  
        }
      }
      #endif
            
      ct = 0;
#ifndef Kaushik_addition
      for(j = 1; j <= data->blksz[k]; j++) for(i = j; i <= data->blksz[k]; 
        i++) if(AGGMAT[i][j] == 1) ct++;
#else
      for(j = 1; j <= data->blksz[k]; j++) for(i = j; i <= data->blksz[k]; 
        i++) if(is_nnz(i, nnz_rows[j], nnz_row_ind[j]) == 1) ct++;
#endif
      

      // Make decision about storage type

      if ((double) 2.0 * ct / (data->blksz[k] * (data->blksz[k] + 1)) >= data->dthresh_dens || data->blksz[k] <= data->dthresh_dim) {

        // Storage is dense. Just set blkptr and blksto.
        
        data->XS_blkptr[k+1] = data->XS_blkptr[k] + data->blksz[k]*data->blksz[k];
        data->XS_blksto[k] = DENSE;

        // Determine absolute position of nonzero ntries of k-th block
        // of Ai in k-th block of S.

        for(i = 0; i <= data->m; i++) {
          if(i == 0) dA = data->C[k];
          else       dA = data->A[i][k];
          if(dA->type == 's')
            for(j = 1; j <= dA->sp->nnz; j++) {
              row = dA->sp->row[j];
              col = dA->sp->col[j];
              dA->sp->XS_in[j] = data->XS_blkptr[k] - 1 + SMATIND(row,col,data->blksz[k]); // code assumes row >= col
            }
        }
 
      }
      else {

        // Storage is sparse, so setup colptr

        MYCALLOC(data->XS_colptr[k], ptrdiff_t, data->blksz[k]+2);

        // For convenience in coding, save colptr locally
        colptr = data->XS_colptr[k];

        // Next determine colptr values
        colptr[1] = 1;
        for(j = 1; j <= data->blksz[k]; j++) {
          ct = 0;
#ifndef Kaushik_addition
          for(i = j; i <= data->blksz[k]; i++) if(AGGMAT[i][j] == 1) ct++;
#else
          for(i = j; i <= data->blksz[k]; i++) 
            if( is_nnz(i, nnz_rows[j], nnz_row_ind[j]) == 1) ct++;
#endif
          colptr[j+1] = colptr[j] + ct;
        }

        // Now we can allocate rowind
        
        MYCALLOC(data->XS_rowind[k], ptrdiff_t, colptr[data->blksz[k]+1]-1 + 1); 

        // Setup rowind
         
#ifndef Kaushik_addition
        for(j = 1; j <= data->blksz[k]; j++) {
          ct = 0;
          for(i = j; i <= data->blksz[k]; i++) if(AGGMAT[i][j] == 1) {
            data->XS_rowind[k][colptr[j] + ct] = i;
            ct++;
          }
        }
#else
        for(j = 1; j <= data->blksz[k]; j++) {
          ct = 0;
          for(i = j; i <= data->blksz[k]; i++) 
            if( is_nnz(i, nnz_rows[j], nnz_row_ind[j]) == 1) {
              data->XS_rowind[k][colptr[j] + ct] = i;
              ct++;
          }
        }

#endif

        // Finally, don't forget to set blkptr and blksto. Based on
        // total number of nonzeros in this block.

        data->XS_blkptr[k+1] = data->XS_blkptr[k] + colptr[data->blksz[k]+1]-1;
        data->XS_blksto[k] = SPARSE;

        // Experimental
        // Determine absolute position of nonzero ntries of k-th block of
        // Ai in k-th block of S
        
#ifndef Kaushik_addition
        ct = 0;
        for(j = 1; j <= data->blksz[k]; j++) 
          for(i = j; i <= data->blksz[k]; i++) if(AGGMAT[i][j] == 1) 
            AGGMAT[i][j] = ++ct;
#endif

        for(i = 0; i <= data->m; i++) {
          if(i == 0) dA = data->C[k];
          else       dA = data->A[i][k];
          if(dA->type == 's')
            for(j = 1; j <= dA->sp->nnz; j++) {
              row = dA->sp->row[j];
              col = dA->sp->col[j];
#ifndef Kaushik_addition
              dA->sp->XS_in[j] = data->XS_blkptr[k] - 1 + AGGMAT[row][col];
              // printf("A(i=%d,k=%d,r=%d,c=%d) --> %d\n", i, k, row, col, 
              // data->XS_blkptr[k] - 1 + AGGMAT[row][col]);
#else
              dA->sp->XS_in[j] = data->XS_blkptr[k] - 1 + agg_id(k, row, col,  
               nnz_rows[col], nnz_row_ind[col], nnz_row_cnt[col], nnz_rows,             
               nnz_row_ind, nnz_row_cnt);
#endif
            }
        }

      }

#ifndef Kaushik_addition
      // Free AGGMAT for next round.
      for(i = 1; i <= data->blksz[k]; i++) { MYFREE(AGGMAT[i]); 
        AGGMAT[i] =   NULL; }
      MYFREE(AGGMAT); AGGMAT = NULL;
#else
      for(col = 0; col<=num_cols; col++){
        free(nnz_row_ind[col]); nnz_row_ind[col]= NULL;
        free(nnz_row_cnt[col]); nnz_row_cnt[col]=NULL;
      }
      free(nnz_row_ind); nnz_row_ind=NULL;
      free(nnz_row_cnt); nnz_row_cnt=NULL;
      free(nnz_rows); nnz_rows=NULL;
      free(CA_col); CA_col=NULL;
      free(CA_row); CA_row=NULL;
#endif

    }

    // Free location of sparse data matrices
    MYFREE(spind);

  }

  // Can now setup X and S

  MYCALLOC(data->S, double, data->XS_blkptr[data->numblk+1]-1 + 1);
  MYCALLOC(data->X, double, data->XS_blkptr[data->numblk+1]-1 + 1);

  // Also setup y

  MYCALLOC(data->y, double, data->m + 1);

  // Construct big, sparse A matrix (including C in 0-th row)
  // Watch out! A potentially weird mixture of Fortran and C indexing.
  
  // First allocate AA_rowptr

  MYCALLOC(data->AA_rowptr, ptrdiff_t, data->m+2); 

  // Now determine AA_rowptr

  data->AA_rowptr[0] = 1;

  for(i = 0; i <= data->m; i++) {

    ctr = 0;

    for(k = 1; k <= data->numblk; k++) {

      if(i == 0) dA = data->C[k];
      else       dA = data->A[i][k];

           if(dA->type == 's') ctr += dA->sp  ->nnz;
      else if(dA->type == 'd') ctr += dA->diag->nnz; 

    }

    data->AA_rowptr[i+1] = data->AA_rowptr[i] + ctr;

  }

  // Next job is to allocate AA_colind and AA_colval_one and _two
  
  MYCALLOC(data->AA_colind, ptrdiff_t, data->AA_rowptr[data->m+1]-1 + 1);
  MYCALLOC(data->AA_colval_one, double, data->AA_rowptr[data->m+1]-1 + 1);
  MYCALLOC(data->AA_colval_two, double, data->AA_rowptr[data->m+1]-1 + 1);

  // Then we fill them (will eventually want to sort this)

  for(i = 0; i <= data->m; i++) {

    ctr = 0;

    for(k = 1; k <= data->numblk; k++) {

      if(i == 0) dA = data->C[k];
      else       dA = data->A[i][k];

      if(dA->type == 's') {
        for(j = 1; j <= dA->sp->nnz; j++) {
          data->AA_colind    [data->AA_rowptr[i] + ctr] = dA->sp->XS_in[j];
          data->AA_colval_one[data->AA_rowptr[i] + ctr] = dA->sp->ent[j];
          if(dA->sp->row[j] == dA->sp->col[j])
            data->AA_colval_two[data->AA_rowptr[i] + ctr] =     dA->sp->ent[j];
          else
            data->AA_colval_two[data->AA_rowptr[i] + ctr] = 2.0*dA->sp->ent[j];
          ctr++;
        }
      }
      else if(dA->type == 'd') {
        for(j = 1; j <= dA->diag->nnz; j++) {
          data->AA_colind    [data->AA_rowptr[i] + ctr] = dA->diag->XS_in[j];
          data->AA_colval_one[data->AA_rowptr[i] + ctr] = dA->diag->ent[j];
          data->AA_colval_two[data->AA_rowptr[i] + ctr] = dA->diag->ent[j];
          ctr++;
        }
      }

    }

  }

  // Setup pointers to find lowrank data matrices. Used in several
  // routines. Right now, we do not set this up for diagonal blocks
 
  ct = 0;
  for (k = 1; k <= data->numblk; k++)
    if (data->blktype[k] == 's') {
                                     if (data->C   [k]->type == 'l') ct++;
      for (i = 1; i <= data->m; i++) if (data->A[i][k]->type == 'l') ct++;
    }
 
  data->lr_num = ct;

  MYCALLOC(data->lr_mat, ptrdiff_t, data->lr_num + 1);
  MYCALLOC(data->lr_blk, ptrdiff_t, data->lr_num + 1);

  ct = 0;
  for (k = 1; k <= data->numblk; k++) 
    if (data->blktype[k] == 's') {
                                     if (data->C   [k]->type == 'l') { data->lr_mat[++ct] = 0; data->lr_blk[ct] = k; }
      for (i = 1; i <= data->m; i++) if (data->A[i][k]->type == 'l') { data->lr_mat[++ct] = i; data->lr_blk[ct] = k; }
    }

  if(ct != data->lr_num) { printf("Problem getting lowrank matrices. (%d != %d)\n", ct, data->lr_num); exit(0); }

  // END: Setup new data structures
  // END: Setup new data structures
  // END: Setup new data structures
  // END: Setup new data structures

  // setup other data structures
  MYCALLOC (data->vio, double, data->m + 1);
  MYCALLOC (data->rank, ptrdiff_t, data->numblk + 1);
  MYCALLOC (data->maxrank, ptrdiff_t, data->numblk + 1);

  // Calculate ranks for nondiag blocks

  for (k = 1; k <= data->numblk; k++)
    data->rank[k] = data->maxrank[k] = maxranks[k-1];

  // Setup global dimension n
  data->n = 0;
  for (k = 1; k <= data->numblk; k++)
    data->n += data->blksz[k];

  // Setup global dimension nr
  data->nr = 0;
  for (k = 1; k <= data->numblk; k++)
    data->nr += data->blksz[k] * data->rank[k];

  // Now can setup data->G for gradient
  MYCALLOC (data->G, double, data->nr + 1);

  // Create global structures
  i = -1;
  for(k = 1; k <= data->numblk; k++)
    if(data->blktype[k] == SDPBLK && data->XS_blksto[k] == DENSE)
      i = mymax(i, data->blksz[k]);

  // Next setup global temp space for evaluation of function with lowrank matrices.

  ct = 0;
  for(h = 1; h <= data->lr_num; h++) {
    k = data->lr_blk[h];
    i = data->lr_mat[h];
    if(i == 0) ct = mymax(ct, data->rank[k]*data->C   [k]->lr->ncol);
    else       ct = mymax(ct, data->rank[k]*data->A[i][k]->lr->ncol);
  }

  MYCALLOC(global_UtB, double, ct+1);
  MYCALLOC(global_VtB, double, ct+1);

  return 1;
}


ptrdiff_t deinitialize (problemdata * data)
{
  ptrdiff_t                 i, j, k;

  MYFREE(global_UtB);
  MYFREE(global_VtB);

  MYFREE (data->vio);
  MYFREE (data->G);
  for (j = 1; j <= data->numblk; j++)
    {
      MYFREE (data->lrind[j]);
      destroydatamat (data->C[j]);
    }
  for (i = 1; i <= data->m; i++) {
    for (j = 1; j <= data->numblk; j++)
      destroydatamat (data->A[i][j]);
    MYFREE (data->A[i]);
  }

  for (k = 1; k <= data->numblk; k++) 
    if(data->blktype[k] == SDPBLK && data->XS_blksto[k] == SPARSE) {
      MYFREE(data->XS_colptr[k]);
      MYFREE(data->XS_rowind[k]); 
    }

  MYFREE(data->XS_blkptr);
  MYFREE(data->XS_blksto);
  MYFREE(data->XS_colptr);
  MYFREE(data->XS_rowind);

  MYFREE(data->S);
  MYFREE(data->X);

  MYFREE(data->y);

  MYFREE(data->AA_rowptr);
  MYFREE(data->AA_colind);
  MYFREE(data->AA_colval_one);
  MYFREE(data->AA_colval_two);

  MYFREE(data->lr_mat);
  MYFREE(data->lr_blk);

  MYFREE (data->lrind);
  MYFREE (data->nlrind);
  MYFREE (data->usind);
  MYFREE (data->nusind);
  MYFREE (data->rank);
  MYFREE (data->maxrank);
  MYFREE (data->C);
  MYFREE (data->A);

  return 1;
}


ptrdiff_t locatetype(problemdata *data, ptrdiff_t blk, char type, ptrdiff_t **passed_ind, ptrdiff_t *passed_nind)
{
  ptrdiff_t i, k;
  ptrdiff_t ct;
  ptrdiff_t *ind, nind;

  k = blk;

  ct = 0;                        if (data->C   [k]->type == type) ct++;
  for (i = 1; i <= data->m; i++) if (data->A[i][k]->type == type) ct++;

  nind = ct;
  MYCALLOC (ind, ptrdiff_t, nind + 1);

  ct = 0;                        if (data->C   [k]->type == type) ind[++ct] = 0;
  for (i = 1; i <= data->m; i++) if (data->A[i][k]->type == type) ind[++ct] = i;
  if (ct != nind) { printf ("locatetype: problem with setting up ind\n"); exit (0); }

  *passed_ind = ind;
  *passed_nind = nind;

  return 0;
}
