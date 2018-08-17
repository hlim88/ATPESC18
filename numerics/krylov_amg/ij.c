/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision$
 ***********************************************************************EHEADER*/

/*--------------------------------------------------------------------------
 * Test driver for unstructured matrix interface (IJ_matrix interface).
 * Do `driver -help' for usage info.
 * This driver started from the driver for parcsr_linear_solvers, and it
 * works by first building a parcsr matrix as before and then "copying"
 * that matrix row-by-row into the IJMatrix interface. AJC 7/99.
 *--------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "_hypre_utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_mv.h"

#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"
#include "HYPRE_krylov.h"

#ifdef HAVE_DSUPERLU
#include "superlu_ddefs.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

HYPRE_Int BuildParLaplacian (HYPRE_Int argc , char *argv [], HYPRE_Int arg_index , HYPRE_ParCSRMatrix *A_ptr );
HYPRE_Int BuildParDifConv (HYPRE_Int argc , char *argv [], HYPRE_Int arg_index , HYPRE_ParCSRMatrix *A_ptr);
HYPRE_Int BuildParFromOneFile (HYPRE_Int argc , char *argv [], HYPRE_Int arg_index , HYPRE_Int num_functions , HYPRE_ParCSRMatrix *A_ptr );
HYPRE_Int BuildRhsParFromOneFile (HYPRE_Int argc , char *argv [], HYPRE_Int arg_index , HYPRE_Int *partitioning , HYPRE_ParVector *b_ptr );
HYPRE_Int BuildParLaplacian9pt (HYPRE_Int argc , char *argv [], HYPRE_Int arg_index , HYPRE_ParCSRMatrix *A_ptr );
HYPRE_Int BuildParLaplacian27pt (HYPRE_Int argc , char *argv [], HYPRE_Int arg_index , HYPRE_ParCSRMatrix *A_ptr );
HYPRE_Int BuildParRotate7pt (HYPRE_Int argc , char *argv [], HYPRE_Int arg_index , HYPRE_ParCSRMatrix *A_ptr );
HYPRE_Int BuildParVarDifConv (HYPRE_Int argc , char *argv [], HYPRE_Int arg_index , HYPRE_ParCSRMatrix *A_ptr , HYPRE_ParVector *rhs_ptr );
HYPRE_Int
hypre_PrintTiming2( const char *heading, HYPRE_Real *wall_time_ptr, MPI_Comm comm );

#ifdef __cplusplus
}
#endif
 
hypre_int
main( hypre_int argc,
      char *argv[] )
{
   HYPRE_Int                 arg_index;
   HYPRE_Int                 print_usage;
   HYPRE_Int                 build_matrix_type;
   HYPRE_Int                 build_matrix_arg_index;
   HYPRE_Int                 build_rhs_type;
   HYPRE_Int                 solver_id;
   HYPRE_Int                 ioutdat;
   HYPRE_Int                 i; 
   HYPRE_Int                 max_levels = 25;
   HYPRE_Int                 num_iterations;
   HYPRE_Int                 max_iter = 1000;
   HYPRE_Int                 mg_max_iter = 100;
   HYPRE_Real          norm;
   HYPRE_Real          final_res_norm;
   HYPRE_Real          setup_time, solve_time;
   void               *object;

   HYPRE_IJMatrix      ij_A = NULL; 
   HYPRE_IJVector      ij_b = NULL;
   HYPRE_IJVector      ij_x = NULL;

   HYPRE_ParCSRMatrix  parcsr_A = NULL;
   HYPRE_ParVector     b = NULL;
   HYPRE_ParVector     x;

   HYPRE_Solver        amg_solver;
   HYPRE_Solver        pcg_solver;
   HYPRE_Solver        pcg_precond=NULL, pcg_precond_gotten;

   HYPRE_Int                 num_procs, myid;
   HYPRE_Int		       num_paths = 1;
   HYPRE_Int		       agg_num_levels = 0;
   HYPRE_Int		       ns_coarse = 1, ns_down = -1, ns_up = -1;

   HYPRE_Int		       time_index;
   HYPRE_Int first_local_row, last_local_row, local_num_rows;
   HYPRE_Int first_local_col, last_local_col, local_num_cols;
   HYPRE_Real *values;

   /* parameters for BoomerAMG */
   HYPRE_Real   strong_threshold;
   HYPRE_Real   trunc_factor;
   HYPRE_Int      P_max_elmts = 4;
   HYPRE_Int      cycle_type;
   HYPRE_Int      coarsen_type = 10;
   HYPRE_Int      num_sweeps = 1;  
   HYPRE_Int      relax_type = -1;   
   HYPRE_Int      relax_coarse = -1;   
   HYPRE_Int      relax_up = -1;   
   HYPRE_Int      relax_down = -1;   
   HYPRE_Int      relax_order = 0;   
   HYPRE_Int      level_w = -1;
   HYPRE_Int      level_ow = -1;
   HYPRE_Int      coarse_threshold = 9;
   HYPRE_Int      min_coarse_size = 0;
   HYPRE_Real   relax_wt = 1.0; 
   HYPRE_Real   relax_wt_level = 1.0; 
   HYPRE_Real   outer_wt = 1.0;
   HYPRE_Real   outer_wt_level = 1.0;
   HYPRE_Real   tol = 1.e-8, pc_tol = 0.;
   HYPRE_Real   max_row_sum = 1.;

   /* parameters for GMRES */
   HYPRE_Int	    k_dim;
   /* interpolation */
   HYPRE_Int      interp_type  = 6; /* default value */
   /* aggressive coarsening */
   HYPRE_Int      agg_P_max_elmts  = 0; /* default value */
   HYPRE_Real   agg_trunc_factor  = 0; /* default value */

   HYPRE_Int      print_system = 0;

   HYPRE_Int rel_change = 0;

   /*-----------------------------------------------------------
    * Initialize some stuff
    *-----------------------------------------------------------*/

   /* Initialize MPI */
   hypre_MPI_Init(&argc, &argv);

   hypre_MPI_Comm_size(hypre_MPI_COMM_WORLD, &num_procs );
   hypre_MPI_Comm_rank(hypre_MPI_COMM_WORLD, &myid );
   
   /*-----------------------------------------------------------
    * Set defaults
    *-----------------------------------------------------------*/
 
   build_matrix_type = 2;
   build_matrix_arg_index = argc;
   build_rhs_type = 2;

   solver_id = 0;

   ioutdat = 0;

   /*-----------------------------------------------------------
    * Parse command line
    *-----------------------------------------------------------*/
 
   print_usage = 0;
   arg_index = 1;

   while ( (arg_index < argc) && (!print_usage) )
   {
      if ( strcmp(argv[arg_index], "-laplacian") == 0 )
      {
         arg_index++;
         build_matrix_type      = 2;
         build_matrix_arg_index = arg_index;
      }
      else if ( strcmp(argv[arg_index], "-9pt") == 0 )
      {
         arg_index++;
         build_matrix_type      = 3;
         build_matrix_arg_index = arg_index;
      }
      else if ( strcmp(argv[arg_index], "-27pt") == 0 )
      {
         arg_index++;
         build_matrix_type      = 4;
         build_matrix_arg_index = arg_index;
      }
      else if ( strcmp(argv[arg_index], "-difconv") == 0 )
      {
         arg_index++;
         build_matrix_type      = 5;
         build_matrix_arg_index = arg_index;
      }
      else if ( strcmp(argv[arg_index], "-vardifconv") == 0 )
      {
         arg_index++;
         build_matrix_type      = 6;
         build_matrix_arg_index = arg_index;
      }
      else if ( strcmp(argv[arg_index], "-rotate") == 0 )
      {
         arg_index++;
         build_matrix_type      = 7;
         build_matrix_arg_index = arg_index;
      } 
      else if ( strcmp(argv[arg_index], "-amg") == 0 )
      {
         arg_index++;
         solver_id = 0;
      }
      else if ( strcmp(argv[arg_index], "-amgpcg") == 0 )
      {
         arg_index++;
         solver_id = 1;
      }
      else if ( strcmp(argv[arg_index], "-pcg") == 0 )
      {
         arg_index++;
         solver_id = 2;
      }
      else if ( strcmp(argv[arg_index], "-amggmres") == 0 )
      {
         arg_index++;
         solver_id = 3;
      }
      else if ( strcmp(argv[arg_index], "-gmres") == 0 )
      {
         arg_index++;
         solver_id = 4;
      }
      else if ( strcmp(argv[arg_index], "-amgbicgstab") == 0 )
      {
         arg_index++;
         solver_id = 9;
      }
      else if ( strcmp(argv[arg_index], "-bicgstab") == 0 )
      {
         arg_index++;
         solver_id = 10;
      }
      else if ( strcmp(argv[arg_index], "-solver") == 0 )
      {
         arg_index++;
         solver_id = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-rhsisone") == 0 )
      {
         arg_index++;
         build_rhs_type      = 2;
      }
      else if ( strcmp(argv[arg_index], "-rhsrand") == 0 )
      {
         arg_index++;
         build_rhs_type      = 3;
      }    
      else if ( strcmp(argv[arg_index], "-xisone") == 0 )
      {
         arg_index++;
         build_rhs_type      = 4;
      }
      else if ( strcmp(argv[arg_index], "-cljp") == 0 )
      {
         arg_index++;
         coarsen_type      = 0;
      }    
      else if ( strcmp(argv[arg_index], "-pmis") == 0 )
      {
         arg_index++;
         coarsen_type      = 8;
      }    
      else if ( strcmp(argv[arg_index], "-hmis") == 0 )
      {
         arg_index++;
         coarsen_type      = 10;
      }    
      else if ( strcmp(argv[arg_index], "-falgout") == 0 )
      {
         arg_index++;
         coarsen_type      = 6;
      }    
      else if ( strcmp(argv[arg_index], "-rlx") == 0 )
      {
         arg_index++;
         relax_type = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-rlx_coarse") == 0 )
      {
         arg_index++;
         relax_coarse = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-rlx_down") == 0 )
      {
         arg_index++;
         relax_down = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-rlx_up") == 0 )
      {
         arg_index++;
         relax_up = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mxl") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-agg_nl") == 0 )
      {
         arg_index++;
         agg_num_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-npaths") == 0 )
      {
         arg_index++;
         num_paths = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ns") == 0 )
      {
         arg_index++;
         num_sweeps = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ns_coarse") == 0 )
      {
         arg_index++;
         ns_coarse = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ns_down") == 0 )
      {
         arg_index++;
         ns_down = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ns_up") == 0 )
      {
         arg_index++;
         ns_up = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-max_iter") == 0 )
      {
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mg_max_iter") == 0 )
      {
         arg_index++;
         mg_max_iter = atoi(argv[arg_index++]);
      }

      else if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         print_usage = 1;
      }
      else
      {
         arg_index++;
      }
   }


   /* defaults for BoomerAMG */
   if (solver_id == 0 || solver_id == 1 || solver_id == 3 || solver_id == 9)
   {
      strong_threshold = 0.25;
      trunc_factor = 0.;
      cycle_type = 1;
      relax_wt = 1.;
      outer_wt = 1.;
   }

   /* defaults for GMRES */

   k_dim = 10;


   arg_index = 0;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-k") == 0 )
      {
         arg_index++;
         k_dim = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-w") == 0 )
      {
         arg_index++;
         relax_wt = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ow") == 0 )
      {
         arg_index++;
         outer_wt = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-th") == 0 )
      {
         arg_index++;
         strong_threshold  = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-CF") == 0 )
      {
         arg_index++;
         relax_order = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tol") == 0 )
      {
         arg_index++;
         tol  = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mxrs") == 0 )
      {
         arg_index++;
         max_row_sum  = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tr") == 0 )
      {
         arg_index++;
         trunc_factor  = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-Pmx") == 0 )
      {
         arg_index++;
         P_max_elmts  = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-iout") == 0 )
      {
         arg_index++;
         ioutdat  = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mu") == 0 )
      {
         arg_index++;
         cycle_type  = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-interptype") == 0 )
      {
         arg_index++;
         interp_type  = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-agg_Pmx") == 0 )
      {
         arg_index++;
         agg_P_max_elmts  = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-agg_tr") == 0 )
      {
         arg_index++;
         agg_trunc_factor  = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-print") == 0 )
      {
         arg_index++;
         print_system = 1;
      }
      else
      {
         arg_index++;
      }
   }

   /*-----------------------------------------------------------
    * Print usage info
    *-----------------------------------------------------------*/
 
   if ( print_usage )
   {
      if ( myid == 0 )
      {
         hypre_printf("\n");
         hypre_printf("Usage: %s [<options>]\n", argv[0]);
         hypre_printf("\n");
         hypre_printf("Choice of Problem:\n");
         hypre_printf("  -laplacian [<options>] : build 7pt 3D laplacian problem (default) \n");
         hypre_printf("  -9pt [<opts>]          : build 9pt 2D laplacian problem\n");
         hypre_printf("  -27pt [<opts>]         : build 27pt 3D laplacian problem\n");
         hypre_printf("  -difconv [<opts>]      : build convection-diffusion problem\n");
         hypre_printf("    -n <nx> <ny> <nz>    : problem size per process \n");
         hypre_printf("    -P <Px> <Py> <Pz>    : process topology\n");
         hypre_printf("    -c <cx> <cy> <cz>    : diffusion coefficients\n");
         hypre_printf("    -a <ax>              : convection coefficients\n");
         hypre_printf("\n");
         hypre_printf("Choice of right hand side:\n");
         hypre_printf("  -rhsrand               : rhs is random vector\n");
         hypre_printf("  -rhsisone              : rhs is vector with unit components (default)\n");
         hypre_printf("  -xisone                : solution of all ones\n");
         hypre_printf("\n");
         hypre_printf("Choice of solver:\n");
         hypre_printf("   -amg                  : AMG only        \n");
         hypre_printf("   -amgpcg               : AMG-PCG        \n");
         hypre_printf("   -pcg                  : diagonally scaled PCG        \n");
         hypre_printf("   -amggmres             : AMG-GMRES with restart k (default 10)     \n");
         hypre_printf("   -gmres                : diagonally scaled GMRES(k) (default k=10)   \n");
         hypre_printf("   -amgbicgstab          : AMG-BiCGSTAB   \n");
         hypre_printf("   -bicgstab             : diagonally scaled BiCGSTAB     \n");
         hypre_printf("   -k   <val>             : dimension Krylov space for GMRES\n");
         hypre_printf("\n");
         hypre_printf("Choice of coarsening scheme:\n");
         hypre_printf("  -cljp                 : CLJP coarsening \n");
         hypre_printf("  -pmis                 : PMIS coarsening \n");
         hypre_printf("  -hmis                 : HMIS coarsening (default)\n");
         hypre_printf("  -falgout              : local Ruge_Stueben followed by CLJP\n");
         hypre_printf("\n");
         hypre_printf("Choice of interpolation:\n");
         hypre_printf("  -interptype  <val>    : set interpolation type\n");
         hypre_printf("       0=Classical modified interpolation  \n");
         hypre_printf("       4=multipass interpolation  \n");
         hypre_printf("       6=extended classical modified interpolation (default) \n");
         hypre_printf("       7=extended (only if no common C neighbor) interpolation  \n");
         hypre_printf("       8=standard interpolation  \n");
         hypre_printf("\n");
         hypre_printf("Choice of smoother:\n");
         hypre_printf("  -rlx  <val>            : relaxation type\n");
         hypre_printf("       0=Weighted Jacobi  \n");
         hypre_printf("  -w   <val>             : set Jacobi relax weight = val (default:1.0)\n");
         hypre_printf("       3=Hybrid (forward) Gauss-Seidel  \n");
         hypre_printf("       4=Hybrid backward Gauss-Seidel  \n");
         hypre_printf("       6=Hybrid symmetric Gauss-Seidel  \n");
         hypre_printf("       8=symmetric L1-Gauss-Seidel  \n");
         hypre_printf("       13=forward L1-Gauss-Seidel  \n");
         hypre_printf("       14=backward L1-Gauss-Seidel  \n");
         hypre_printf("       18=L1-Jacobi (may be used with -CF) \n");
         hypre_printf("  -rlx_down    <val>       : set relaxation type for down cycle\n");
         hypre_printf("  -rlx_up      <val>       : set relaxation type for up cycle\n");
         hypre_printf("  -CF  : relaxes first C then F-points on the down cycle, and F then C up\n");
         hypre_printf("\n");
         hypre_printf("  -ns <val>              : Use <val> sweeps on each level\n");
         hypre_printf("  -mu   <val>            : set AMG cycles (1=V, 2=W, etc.)\n"); 
         hypre_printf("  -th   <val>            : set AMG threshold Theta = val \n");
         hypre_printf("  -tr   <val>            : set AMG interpolation truncation factor = val \n");
         hypre_printf("  -Pmx  <val>            : set maximal no. of elmts per row for AMG interpolation (default: 4)\n");
         hypre_printf("  -mxrs <val>            : set AMG maximum row sum threshold for dependency weakening \n");

         hypre_printf("\n");


         hypre_printf("  -mxl  <val>            : maximum number of levels (AMG)\n");
         hypre_printf("  -tol  <val>            : set solver convergence tolerance = val\n");
         hypre_printf("  -max_iter  <val>       : set max iterations\n");
         hypre_printf("  -mg_max_iter  <val>    : set max iterations for mg solvers\n");
         hypre_printf("  -agg_nl  <val>         : set number of aggressive coarsening levels (default:0)\n");
         hypre_printf("  -np  <val>             : set number of paths of length 2 for aggr. coarsening\n");
         hypre_printf("\n");
         hypre_printf("  -iout <val>            : set output flag\n");
         hypre_printf("       0=no output    1=matrix stats\n"); 
         hypre_printf("       2=cycle stats  3=matrix & cycle stats\n"); 
         hypre_printf("\n");  
         hypre_printf("  -print                 : print out the system\n");
         hypre_printf("\n");
      }

      goto final;
   }

   /*-----------------------------------------------------------
    * Print driver parameters
    *-----------------------------------------------------------*/
 
   if (myid == 0)
   {
      hypre_printf("Running with these driver parameters:\n");
      hypre_printf("  solver ID    = %d\n\n", solver_id);
   }

   /*-----------------------------------------------------------
    * Set up matrix
    *-----------------------------------------------------------*/

   time_index = hypre_InitializeTiming("Spatial Operator");
   hypre_BeginTiming(time_index);
   if ( build_matrix_type == 2 )
   {
      BuildParLaplacian(argc, argv, build_matrix_arg_index, &parcsr_A);
   }
   else if ( build_matrix_type == 3 )
   {
      BuildParLaplacian9pt(argc, argv, build_matrix_arg_index, &parcsr_A);
   }
   else if ( build_matrix_type == 4 )
   {
      BuildParLaplacian27pt(argc, argv, build_matrix_arg_index, &parcsr_A);
   }
   else if ( build_matrix_type == 5 )
   {
      BuildParDifConv(argc, argv, build_matrix_arg_index, &parcsr_A);
   }
   else if ( build_matrix_type == 6 )
   {
      BuildParVarDifConv(argc, argv, build_matrix_arg_index, &parcsr_A, &b);
      build_rhs_type      = 6;
   }
   else if ( build_matrix_type == 7 )
   {
      BuildParRotate7pt(argc, argv, build_matrix_arg_index, &parcsr_A);
   }

   else
   {
      hypre_printf("You have asked for an unsupported problem with\n");
      hypre_printf("build_matrix_type = %d.\n", build_matrix_type);
      return(-1);
   }

   hypre_EndTiming(time_index);
   hypre_PrintTiming("Generate Matrix", hypre_MPI_COMM_WORLD);
   hypre_FinalizeTiming(time_index);
   hypre_ClearTiming();

   /*-----------------------------------------------------------
    * Set up the RHS and initial guess
    *-----------------------------------------------------------*/

   time_index = hypre_InitializeTiming("RHS and Initial Guess");
   hypre_BeginTiming(time_index);

   HYPRE_ParCSRMatrixGetLocalRange( parcsr_A,
                                          &first_local_row, &last_local_row ,
                                          &first_local_col, &last_local_col );

   local_num_rows = last_local_row - first_local_row + 1;
   local_num_cols = last_local_col - first_local_col + 1;

   if ( build_rhs_type == 2 )
   {
      if (myid == 0)
      {
         hypre_printf("  RHS vector has unit components\n");
         hypre_printf("  Initial guess is 0\n");
      }

      /* RHS */
      HYPRE_IJVectorCreate(hypre_MPI_COMM_WORLD, first_local_row, last_local_row, &ij_b);
      HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_b);

      values = hypre_CTAlloc(HYPRE_Real,  local_num_rows, HYPRE_MEMORY_HOST);
      for (i = 0; i < local_num_rows; i++)
         values[i] = 1.0;
      HYPRE_IJVectorSetValues(ij_b, local_num_rows, NULL, values);
      hypre_TFree(values, HYPRE_MEMORY_HOST);

      HYPRE_IJVectorGetObject( ij_b, &object );
      b = (HYPRE_ParVector) object;

      /* Initial guess */
      HYPRE_IJVectorCreate(hypre_MPI_COMM_WORLD, first_local_col, last_local_col, &ij_x);
      HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_x);

      values = hypre_CTAlloc(HYPRE_Real,  local_num_cols, HYPRE_MEMORY_HOST);
      for (i = 0; i < local_num_cols; i++)
         values[i] = 0.;
      HYPRE_IJVectorSetValues(ij_x, local_num_cols, NULL, values);
      hypre_TFree(values, HYPRE_MEMORY_HOST);

      HYPRE_IJVectorGetObject( ij_x, &object );
      x = (HYPRE_ParVector) object;
   }
   else if ( build_rhs_type == 3 )
   {
      if (myid == 0)
      {
         hypre_printf("  RHS vector has random components and unit 2-norm\n");
         hypre_printf("  Initial guess is 0\n");
      }

      /* RHS */
      HYPRE_IJVectorCreate(hypre_MPI_COMM_WORLD, first_local_row, last_local_row, &ij_b); 
      HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_b);
      HYPRE_IJVectorGetObject( ij_b, &object );
      b = (HYPRE_ParVector) object;

      /* For purposes of this test, HYPRE_ParVector functions are used, but
         these are not necessary.  For a clean use of the interface, the user
         "should" modify components of ij_x by using functions
         HYPRE_IJVectorSetValues or HYPRE_IJVectorAddToValues */

      HYPRE_ParVectorSetRandomValues(b, 22775);
      HYPRE_ParVectorInnerProd(b,b,&norm);
      norm = 1./sqrt(norm);
      HYPRE_ParVectorScale(norm, b);      

      /* Initial guess */
      HYPRE_IJVectorCreate(hypre_MPI_COMM_WORLD, first_local_col, last_local_col, &ij_x);
      HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_x);

      values = hypre_CTAlloc(HYPRE_Real,  local_num_cols, HYPRE_MEMORY_HOST);
      for (i = 0; i < local_num_cols; i++)
         values[i] = 0.;
      HYPRE_IJVectorSetValues(ij_x, local_num_cols, NULL, values);
      hypre_TFree(values, HYPRE_MEMORY_HOST);

      HYPRE_IJVectorGetObject( ij_x, &object );
      x = (HYPRE_ParVector) object;
   }
   else if ( build_rhs_type == 4 )
   {
      if (myid == 0)
      {
         hypre_printf("  RHS vector set for solution with unit components\n");
         hypre_printf("  Initial guess is 0\n");
      }

      /* Temporary use of solution vector */
      HYPRE_IJVectorCreate(hypre_MPI_COMM_WORLD, first_local_col, last_local_col, &ij_x);
      HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_x);

      values = hypre_CTAlloc(HYPRE_Real,  local_num_cols, HYPRE_MEMORY_HOST);
      for (i = 0; i < local_num_cols; i++)
         values[i] = 1.;
      HYPRE_IJVectorSetValues(ij_x, local_num_cols, NULL, values);
      hypre_TFree(values, HYPRE_MEMORY_HOST);

      HYPRE_IJVectorGetObject( ij_x, &object );
      x = (HYPRE_ParVector) object;

      /* RHS */
      HYPRE_IJVectorCreate(hypre_MPI_COMM_WORLD, first_local_row, last_local_row, &ij_b); 
      HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_b);
      HYPRE_IJVectorGetObject( ij_b, &object );
      b = (HYPRE_ParVector) object;

      HYPRE_ParCSRMatrixMatvec(1.,parcsr_A,x,0.,b);

      /* Initial guess */
      values = hypre_CTAlloc(HYPRE_Real,  local_num_cols, HYPRE_MEMORY_HOST);
      for (i = 0; i < local_num_cols; i++)
         values[i] = 0.;
      HYPRE_IJVectorSetValues(ij_x, local_num_cols, NULL, values);
      hypre_TFree(values, HYPRE_MEMORY_HOST);
   }
   else if ( build_rhs_type == 6) 
   {
      ij_b = NULL;
      /* Initial guess */
      HYPRE_IJVectorCreate(hypre_MPI_COMM_WORLD, first_local_col, last_local_col, &ij_x);
      HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_x);

      values = hypre_CTAlloc(HYPRE_Real,  local_num_cols, HYPRE_MEMORY_HOST);
      for (i = 0; i < local_num_cols; i++)
         values[i] = 0.;
      HYPRE_IJVectorSetValues(ij_x, local_num_cols, NULL, values);
      hypre_TFree(values, HYPRE_MEMORY_HOST);

      HYPRE_IJVectorGetObject( ij_x, &object );
      x = (HYPRE_ParVector) object;
   }

   hypre_EndTiming(time_index);
   hypre_PrintTiming("IJ Vector Setup", hypre_MPI_COMM_WORLD);
   hypre_FinalizeTiming(time_index);
   hypre_ClearTiming();
   
   /*-----------------------------------------------------------
    * Print out the system and initial guess
    *-----------------------------------------------------------*/

   if (print_system)
   {
      if (ij_A)
      {
         HYPRE_IJMatrixPrint(ij_A, "IJ.out.A");
      }
      else if (parcsr_A)
      {
         hypre_ParCSRMatrixPrintIJ(parcsr_A, 0, 0, "IJ.out.A");
      }
      if (ij_b)
      {
         HYPRE_IJVectorPrint(ij_b, "IJ.out.b");
      }
      else if (b)
      {
         HYPRE_ParVectorPrint(b, "ParVec.out.b");
      }
      HYPRE_IJVectorPrint(ij_x, "IJ.out.x0");

      /* HYPRE_ParCSRMatrixPrint( parcsr_A, "new_mat.A" );*/
   }

   if (solver_id == 0)
   {
      if (myid == 0) hypre_printf("Solver:  AMG\n");
      time_index = hypre_InitializeTiming("BoomerAMG Setup");
      hypre_BeginTiming(time_index);

      HYPRE_BoomerAMGCreate(&amg_solver); 
     
      HYPRE_BoomerAMGSetInterpType(amg_solver, interp_type);
      HYPRE_BoomerAMGSetCoarsenType(amg_solver, coarsen_type);
      HYPRE_BoomerAMGSetTol(amg_solver, tol);
      HYPRE_BoomerAMGSetStrongThreshold(amg_solver, strong_threshold);
      HYPRE_BoomerAMGSetMaxCoarseSize(amg_solver, coarse_threshold);
      HYPRE_BoomerAMGSetMinCoarseSize(amg_solver, min_coarse_size);
      HYPRE_BoomerAMGSetTruncFactor(amg_solver, trunc_factor);
      HYPRE_BoomerAMGSetPMaxElmts(amg_solver, P_max_elmts);
      HYPRE_BoomerAMGSetPrintLevel(amg_solver, ioutdat);
      HYPRE_BoomerAMGSetCycleType(amg_solver, cycle_type);
      HYPRE_BoomerAMGSetNumSweeps(amg_solver, num_sweeps);
      if (relax_type > -1) HYPRE_BoomerAMGSetRelaxType(amg_solver, relax_type);
      if (relax_down > -1)
         HYPRE_BoomerAMGSetCycleRelaxType(amg_solver, relax_down, 1);
      if (relax_up > -1)
         HYPRE_BoomerAMGSetCycleRelaxType(amg_solver, relax_up, 2);
      if (relax_coarse > -1)
         HYPRE_BoomerAMGSetCycleRelaxType(amg_solver, relax_coarse, 3);
      HYPRE_BoomerAMGSetRelaxOrder(amg_solver, relax_order);
      HYPRE_BoomerAMGSetRelaxWt(amg_solver, relax_wt);
      HYPRE_BoomerAMGSetOuterWt(amg_solver, outer_wt);
      HYPRE_BoomerAMGSetMaxLevels(amg_solver, max_levels);
      if (level_w > -1)
         HYPRE_BoomerAMGSetLevelRelaxWt(amg_solver, relax_wt_level, level_w);
      if (level_ow > -1)
         HYPRE_BoomerAMGSetLevelOuterWt(amg_solver, outer_wt_level, level_ow);
      HYPRE_BoomerAMGSetMaxRowSum(amg_solver, max_row_sum);
      HYPRE_BoomerAMGSetAggNumLevels(amg_solver, agg_num_levels);
      HYPRE_BoomerAMGSetAggTruncFactor(amg_solver, agg_trunc_factor);
      HYPRE_BoomerAMGSetAggPMaxElmts(amg_solver, agg_P_max_elmts);
      HYPRE_BoomerAMGSetNumPaths(amg_solver, num_paths);
      HYPRE_BoomerAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
      if (ns_down > -1)
      {
         HYPRE_BoomerAMGSetCycleNumSweeps(amg_solver, ns_down,   1);
      }
      if (ns_up > -1)
      {
         HYPRE_BoomerAMGSetCycleNumSweeps(amg_solver, ns_up,     2);
      }
      HYPRE_BoomerAMGSetMaxIter(amg_solver, mg_max_iter);

      HYPRE_BoomerAMGSetup(amg_solver, parcsr_A, b, x);

      hypre_EndTiming(time_index);
      hypre_PrintTiming2("Setup phase times", &setup_time, hypre_MPI_COMM_WORLD);
      hypre_FinalizeTiming(time_index);
      hypre_ClearTiming();
 
      time_index = hypre_InitializeTiming("BoomerAMG Solve");
      hypre_BeginTiming(time_index);

      HYPRE_BoomerAMGSolve(amg_solver, parcsr_A, b, x);

      hypre_EndTiming(time_index);
      hypre_PrintTiming2("Solve phase times", &solve_time, hypre_MPI_COMM_WORLD);
      hypre_FinalizeTiming(time_index);
      hypre_ClearTiming();

      HYPRE_BoomerAMGGetNumIterations(amg_solver, &num_iterations);
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(amg_solver, &final_res_norm);

      if (myid == 0)
      {
         hypre_printf("\n");
         hypre_printf("BoomerAMG Iterations = %d\n", num_iterations);
         hypre_printf("Final Relative Residual Norm = %e\n", final_res_norm);
         hypre_printf("Total time = %f\n", (solve_time+setup_time));
         hypre_printf("\n");
      }

      HYPRE_BoomerAMGDestroy(amg_solver);
   }

   /*-----------------------------------------------------------
    * Solve the system using PCG 
    *-----------------------------------------------------------*/
   if (solver_id == 1 || solver_id == 2)
   {
      time_index = hypre_InitializeTiming("PCG Setup");
      hypre_BeginTiming(time_index);
 
      HYPRE_ParCSRPCGCreate(hypre_MPI_COMM_WORLD, &pcg_solver);
      HYPRE_PCGSetMaxIter(pcg_solver, max_iter);
      HYPRE_PCGSetTol(pcg_solver, tol);
      HYPRE_PCGSetTwoNorm(pcg_solver, 1);
      HYPRE_PCGSetRelChange(pcg_solver, rel_change);
      HYPRE_PCGSetPrintLevel(pcg_solver, ioutdat);

      if (solver_id == 1)
      {
         /* use BoomerAMG as preconditioner */
         if (myid == 0) hypre_printf("Solver: AMG-PCG\n");
         HYPRE_BoomerAMGCreate(&pcg_precond); 
         HYPRE_BoomerAMGSetInterpType(pcg_precond, interp_type);
         HYPRE_BoomerAMGSetTol(pcg_precond, pc_tol);
         HYPRE_BoomerAMGSetCoarsenType(pcg_precond, coarsen_type);
         HYPRE_BoomerAMGSetStrongThreshold(pcg_precond, strong_threshold);
         HYPRE_BoomerAMGSetMaxCoarseSize(pcg_precond, coarse_threshold);
         HYPRE_BoomerAMGSetMinCoarseSize(pcg_precond, min_coarse_size);
         HYPRE_BoomerAMGSetTruncFactor(pcg_precond, trunc_factor);
         HYPRE_BoomerAMGSetPMaxElmts(pcg_precond, P_max_elmts);
         HYPRE_BoomerAMGSetPrintLevel(pcg_precond, ioutdat);
         HYPRE_BoomerAMGSetMaxIter(pcg_precond, 1);
         HYPRE_BoomerAMGSetCycleType(pcg_precond, cycle_type);
         HYPRE_BoomerAMGSetNumSweeps(pcg_precond, num_sweeps);
         if (relax_type > -1) HYPRE_BoomerAMGSetRelaxType(pcg_precond, relax_type);
         if (relax_down > -1)
            HYPRE_BoomerAMGSetCycleRelaxType(pcg_precond, relax_down, 1);
         if (relax_up > -1)
            HYPRE_BoomerAMGSetCycleRelaxType(pcg_precond, relax_up, 2);
         if (relax_coarse > -1)
            HYPRE_BoomerAMGSetCycleRelaxType(pcg_precond, relax_coarse, 3);
         HYPRE_BoomerAMGSetRelaxOrder(pcg_precond, relax_order);
         HYPRE_BoomerAMGSetRelaxWt(pcg_precond, relax_wt);
         HYPRE_BoomerAMGSetOuterWt(pcg_precond, outer_wt);
         if (level_w > -1)
            HYPRE_BoomerAMGSetLevelRelaxWt(pcg_precond, relax_wt_level,level_w);
         if (level_ow > -1)
            HYPRE_BoomerAMGSetLevelOuterWt(pcg_precond,outer_wt_level,level_ow);
         HYPRE_BoomerAMGSetMaxLevels(pcg_precond, max_levels);
         HYPRE_BoomerAMGSetMaxRowSum(pcg_precond, max_row_sum);
         HYPRE_BoomerAMGSetAggNumLevels(pcg_precond, agg_num_levels);
         HYPRE_BoomerAMGSetAggTruncFactor(pcg_precond, agg_trunc_factor);
         HYPRE_BoomerAMGSetAggPMaxElmts(pcg_precond, agg_P_max_elmts);
         HYPRE_BoomerAMGSetNumPaths(pcg_precond, num_paths);
         HYPRE_BoomerAMGSetCycleNumSweeps(pcg_precond, ns_coarse, 3);
         HYPRE_PCGSetMaxIter(pcg_solver, mg_max_iter);
         HYPRE_PCGSetPrecond(pcg_solver,
                             (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                             (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                             pcg_precond);
      }
      else if (solver_id == 2)
      {
         
         /* use diagonal scaling as preconditioner */
         if (myid == 0) hypre_printf("Solver: DS-PCG\n");
         pcg_precond = NULL;

         HYPRE_PCGSetPrecond(pcg_solver,
                             (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScale,
                             (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScaleSetup,
                             pcg_precond);
      }
 
      HYPRE_PCGGetPrecond(pcg_solver, &pcg_precond_gotten);
      if (pcg_precond_gotten !=  pcg_precond)
      {
         hypre_printf("HYPRE_ParCSRPCGGetPrecond got bad precond\n");
         return(-1);
      }
      else 
         if (myid == 0)
            hypre_printf("HYPRE_ParCSRPCGGetPrecond got good precond\n");
      HYPRE_PCGSetup(pcg_solver, (HYPRE_Matrix)parcsr_A, 
                     (HYPRE_Vector)b, (HYPRE_Vector)x);
      hypre_EndTiming(time_index);
      hypre_PrintTiming2("Setup phase times", &setup_time, hypre_MPI_COMM_WORLD);
      hypre_FinalizeTiming(time_index);
      hypre_ClearTiming();
   
      time_index = hypre_InitializeTiming("PCG Solve");
      hypre_BeginTiming(time_index);

      HYPRE_PCGSolve(pcg_solver, (HYPRE_Matrix)parcsr_A, 
                     (HYPRE_Vector)b, (HYPRE_Vector)x);

      hypre_EndTiming(time_index);
      hypre_PrintTiming2("Solve phase times", &solve_time, hypre_MPI_COMM_WORLD);
      hypre_FinalizeTiming(time_index);
      hypre_ClearTiming();
 
      HYPRE_PCGGetNumIterations(pcg_solver, &num_iterations);
      HYPRE_PCGGetFinalRelativeResidualNorm(pcg_solver, &final_res_norm);


      HYPRE_ParCSRPCGDestroy(pcg_solver);
 
      if (solver_id == 1)
      {
         HYPRE_BoomerAMGDestroy(pcg_precond);
      }

      if (myid == 0)
      {
         hypre_printf("\n");
         hypre_printf("Iterations = %d\n", num_iterations);
         hypre_printf("Final Relative Residual Norm = %e\n", final_res_norm);
         hypre_printf("Total time = %f\n", (solve_time+setup_time));
         hypre_printf("\n");
      }
 
   }

   /*-----------------------------------------------------------
    * Solve the system using GMRES 
    *-----------------------------------------------------------*/

   if (solver_id == 3 || solver_id == 4 )
   {
      time_index = hypre_InitializeTiming("GMRES Setup");
      hypre_BeginTiming(time_index);
 
      HYPRE_ParCSRGMRESCreate(hypre_MPI_COMM_WORLD, &pcg_solver);
      HYPRE_GMRESSetKDim(pcg_solver, k_dim);
      HYPRE_GMRESSetMaxIter(pcg_solver, max_iter);
      HYPRE_GMRESSetTol(pcg_solver, tol);
      HYPRE_GMRESSetLogging(pcg_solver, 1);
      HYPRE_GMRESSetPrintLevel(pcg_solver, ioutdat);
      HYPRE_GMRESSetRelChange(pcg_solver, rel_change);

      if (solver_id == 3)
      {
         /* use BoomerAMG as preconditioner */
         if (myid == 0) hypre_printf("Solver: AMG-GMRES\n");

         HYPRE_BoomerAMGCreate(&pcg_precond); 

         HYPRE_BoomerAMGSetInterpType(pcg_precond, interp_type);
         HYPRE_BoomerAMGSetTol(pcg_precond, pc_tol);
         HYPRE_BoomerAMGSetCoarsenType(pcg_precond, coarsen_type);
         HYPRE_BoomerAMGSetStrongThreshold(pcg_precond, strong_threshold);
         HYPRE_BoomerAMGSetMaxCoarseSize(pcg_precond, coarse_threshold);
         HYPRE_BoomerAMGSetMinCoarseSize(pcg_precond, min_coarse_size);
         HYPRE_BoomerAMGSetTruncFactor(pcg_precond, trunc_factor);
         HYPRE_BoomerAMGSetPMaxElmts(pcg_precond, P_max_elmts);
         HYPRE_BoomerAMGSetPrintLevel(pcg_precond, ioutdat);
         HYPRE_BoomerAMGSetMaxIter(pcg_precond, 1);
         HYPRE_BoomerAMGSetCycleType(pcg_precond, cycle_type);
         HYPRE_BoomerAMGSetNumSweeps(pcg_precond, num_sweeps);
         if (relax_type > -1) HYPRE_BoomerAMGSetRelaxType(pcg_precond, relax_type);
         if (relax_down > -1)
            HYPRE_BoomerAMGSetCycleRelaxType(pcg_precond, relax_down, 1);
         if (relax_up > -1)
            HYPRE_BoomerAMGSetCycleRelaxType(pcg_precond, relax_up, 2);
         if (relax_coarse > -1)
            HYPRE_BoomerAMGSetCycleRelaxType(pcg_precond, relax_coarse, 3);
         HYPRE_BoomerAMGSetRelaxOrder(pcg_precond, relax_order);
         HYPRE_BoomerAMGSetRelaxWt(pcg_precond, relax_wt);
         HYPRE_BoomerAMGSetOuterWt(pcg_precond, outer_wt);
         if (level_w > -1)
            HYPRE_BoomerAMGSetLevelRelaxWt(pcg_precond, relax_wt_level,level_w);
         if (level_ow > -1)
            HYPRE_BoomerAMGSetLevelOuterWt(pcg_precond,outer_wt_level,level_ow);
         HYPRE_BoomerAMGSetMaxLevels(pcg_precond, max_levels);
         HYPRE_BoomerAMGSetMaxRowSum(pcg_precond, max_row_sum);
         HYPRE_BoomerAMGSetAggNumLevels(pcg_precond, agg_num_levels);
         HYPRE_BoomerAMGSetAggTruncFactor(pcg_precond, agg_trunc_factor);
         HYPRE_BoomerAMGSetAggPMaxElmts(pcg_precond, agg_P_max_elmts);
         HYPRE_BoomerAMGSetNumPaths(pcg_precond, num_paths);
         HYPRE_BoomerAMGSetCycleNumSweeps(pcg_precond, ns_coarse, 3);
         if (ns_down > -1)
         {
            HYPRE_BoomerAMGSetCycleNumSweeps(pcg_precond, ns_down,   1);
         }
         if (ns_up > -1)
         {
            HYPRE_BoomerAMGSetCycleNumSweeps(pcg_precond, ns_up,     2);
         }
         HYPRE_GMRESSetMaxIter(pcg_solver, mg_max_iter);
         HYPRE_GMRESSetPrecond(pcg_solver,
                               (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                               (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                               pcg_precond);
      }
      else if (solver_id == 4)
      {
         /* use diagonal scaling as preconditioner */
         if (myid == 0) hypre_printf("Solver: DS-GMRES\n");
         pcg_precond = NULL;

         HYPRE_GMRESSetPrecond(pcg_solver,
                               (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScale,
                               (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScaleSetup,
                               pcg_precond);
      }
 
      HYPRE_GMRESGetPrecond(pcg_solver, &pcg_precond_gotten);
      if (pcg_precond_gotten != pcg_precond)
      {
         hypre_printf("HYPRE_GMRESGetPrecond got bad precond\n");
         return(-1);
      }
      else
         if (myid == 0)
            hypre_printf("HYPRE_GMRESGetPrecond got good precond\n");
      HYPRE_GMRESSetup
         (pcg_solver, (HYPRE_Matrix)parcsr_A, (HYPRE_Vector)b, (HYPRE_Vector)x);
 
      hypre_EndTiming(time_index);
      hypre_PrintTiming2("Setup phase times", &setup_time, hypre_MPI_COMM_WORLD);
      hypre_FinalizeTiming(time_index);
      hypre_ClearTiming();
   
      time_index = hypre_InitializeTiming("GMRES Solve");
      hypre_BeginTiming(time_index);
 
      HYPRE_GMRESSolve
         (pcg_solver, (HYPRE_Matrix)parcsr_A, (HYPRE_Vector)b, (HYPRE_Vector)x);
 
      hypre_EndTiming(time_index);
      hypre_PrintTiming2("Solve phase times", &solve_time, hypre_MPI_COMM_WORLD);
      hypre_FinalizeTiming(time_index);
      hypre_ClearTiming();
 
      HYPRE_GMRESGetNumIterations(pcg_solver, &num_iterations);
      HYPRE_GMRESGetFinalRelativeResidualNorm(pcg_solver,&final_res_norm);

      HYPRE_ParCSRGMRESDestroy(pcg_solver);
 
      if (solver_id == 3)
      {
         HYPRE_BoomerAMGDestroy(pcg_precond);
      }

      if (myid == 0)
      {
         hypre_printf("\n");
         hypre_printf("GMRES Iterations = %d\n", num_iterations);
         hypre_printf("Final GMRES Relative Residual Norm = %e\n", final_res_norm);
         hypre_printf("Total time = %f\n", (solve_time+setup_time));
         hypre_printf("\n");
      }
   }

   /*-----------------------------------------------------------
    * Solve the system using BiCGSTAB 
    *-----------------------------------------------------------*/

   if (solver_id == 9 || solver_id == 10 )
   {
      time_index = hypre_InitializeTiming("BiCGSTAB Setup");
      hypre_BeginTiming(time_index);
 
      HYPRE_ParCSRBiCGSTABCreate(hypre_MPI_COMM_WORLD, &pcg_solver);
      HYPRE_BiCGSTABSetMaxIter(pcg_solver, max_iter);
      HYPRE_BiCGSTABSetTol(pcg_solver, tol);
      HYPRE_BiCGSTABSetLogging(pcg_solver, ioutdat);
      HYPRE_BiCGSTABSetPrintLevel(pcg_solver, ioutdat);
 
      if (solver_id == 9)
      {
         /* use BoomerAMG as preconditioner */
         if (myid == 0) hypre_printf("Solver: AMG-BiCGSTAB\n");
         HYPRE_BoomerAMGCreate(&pcg_precond); 
         HYPRE_BoomerAMGSetInterpType(pcg_precond, interp_type);
         HYPRE_BoomerAMGSetTol(pcg_precond, pc_tol);
         HYPRE_BoomerAMGSetCoarsenType(pcg_precond, coarsen_type);
         HYPRE_BoomerAMGSetStrongThreshold(pcg_precond, strong_threshold);
         HYPRE_BoomerAMGSetMaxCoarseSize(pcg_precond, coarse_threshold);
         HYPRE_BoomerAMGSetMinCoarseSize(pcg_precond, min_coarse_size);
         HYPRE_BoomerAMGSetTruncFactor(pcg_precond, trunc_factor);
         HYPRE_BoomerAMGSetPMaxElmts(pcg_precond, P_max_elmts);
         HYPRE_BoomerAMGSetPrintLevel(pcg_precond, ioutdat);
         HYPRE_BoomerAMGSetMaxIter(pcg_precond, 1);
         HYPRE_BoomerAMGSetCycleType(pcg_precond, cycle_type);
         HYPRE_BoomerAMGSetNumSweeps(pcg_precond, num_sweeps);
         if (relax_type > -1) HYPRE_BoomerAMGSetRelaxType(pcg_precond, relax_type);
         if (relax_down > -1)
            HYPRE_BoomerAMGSetCycleRelaxType(pcg_precond, relax_down, 1);
         if (relax_up > -1)
            HYPRE_BoomerAMGSetCycleRelaxType(pcg_precond, relax_up, 2);
         if (relax_coarse > -1)
            HYPRE_BoomerAMGSetCycleRelaxType(pcg_precond, relax_coarse, 3);
         HYPRE_BoomerAMGSetRelaxOrder(pcg_precond, relax_order);
         HYPRE_BoomerAMGSetRelaxWt(pcg_precond, relax_wt);
         HYPRE_BoomerAMGSetOuterWt(pcg_precond, outer_wt);
         if (level_w > -1)
            HYPRE_BoomerAMGSetLevelRelaxWt(pcg_precond, relax_wt_level,level_w);
         if (level_ow > -1)
            HYPRE_BoomerAMGSetLevelOuterWt(pcg_precond,outer_wt_level,level_ow);
         HYPRE_BoomerAMGSetMaxLevels(pcg_precond, max_levels);
         HYPRE_BoomerAMGSetMaxRowSum(pcg_precond, max_row_sum);
         HYPRE_BoomerAMGSetAggNumLevels(pcg_precond, agg_num_levels);
         HYPRE_BoomerAMGSetAggTruncFactor(pcg_precond, agg_trunc_factor);
         HYPRE_BoomerAMGSetAggPMaxElmts(pcg_precond, agg_P_max_elmts);
         HYPRE_BoomerAMGSetNumPaths(pcg_precond, num_paths);
         HYPRE_BoomerAMGSetCycleNumSweeps(pcg_precond, ns_coarse, 3);
         HYPRE_BiCGSTABSetMaxIter(pcg_solver, mg_max_iter);
         HYPRE_BiCGSTABSetPrecond(pcg_solver,
                                  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                                  pcg_precond);
      }
      else if (solver_id == 10)
      {
         /* use diagonal scaling as preconditioner */
         if (myid == 0) hypre_printf("Solver: DS-BiCGSTAB\n");
         pcg_precond = NULL;

         HYPRE_BiCGSTABSetPrecond(pcg_solver,
                                  (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScale,
                                  (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScaleSetup,
                                  pcg_precond);
      }
      HYPRE_BiCGSTABSetup(pcg_solver, (HYPRE_Matrix)parcsr_A, 
                          (HYPRE_Vector)b, (HYPRE_Vector)x);
 
      hypre_EndTiming(time_index);
      hypre_PrintTiming2("Setup phase times", &setup_time, hypre_MPI_COMM_WORLD);
      hypre_FinalizeTiming(time_index);
      hypre_ClearTiming();
   
      time_index = hypre_InitializeTiming("BiCGSTAB Solve");
      hypre_BeginTiming(time_index);
 
      HYPRE_BiCGSTABSolve(pcg_solver, (HYPRE_Matrix)parcsr_A, 
                          (HYPRE_Vector)b, (HYPRE_Vector)x);
 
      hypre_EndTiming(time_index);
      hypre_PrintTiming2("Solve phase times", &solve_time, hypre_MPI_COMM_WORLD);
      hypre_FinalizeTiming(time_index);
      hypre_ClearTiming();
 
      HYPRE_BiCGSTABGetNumIterations(pcg_solver, &num_iterations);
      HYPRE_BiCGSTABGetFinalRelativeResidualNorm(pcg_solver,&final_res_norm);

      HYPRE_ParCSRBiCGSTABDestroy(pcg_solver);
 
      if (solver_id == 9)
      {
         HYPRE_BoomerAMGDestroy(pcg_precond);
      }

      if (myid == 0)
      {
         hypre_printf("\n");
         hypre_printf("BiCGSTAB Iterations = %d\n", num_iterations);
         hypre_printf("Final BiCGSTAB Relative Residual Norm = %e\n", final_res_norm);
         hypre_printf("Total time = %f\n", (solve_time+setup_time));
         hypre_printf("\n");
      }
   }

   /*-----------------------------------------------------------
    * Print the solution and other info
    *-----------------------------------------------------------*/

   /* RDF: Why is this here? */
   /*if (!(build_rhs_type ==1 || build_rhs_type ==7))
      HYPRE_IJVectorGetObjectType(ij_b, &j);*/

   if (print_system)
   {
      HYPRE_IJVectorPrint(ij_x, "IJ.out.x");
   }

   /*-----------------------------------------------------------
    * Finalize things
    *-----------------------------------------------------------*/

   HYPRE_ParCSRMatrixDestroy(parcsr_A);

   /* for build_rhs_type = 1 or 7, we did not create ij_b  - just b*/
   if (build_rhs_type ==1 || build_rhs_type ==7 || build_rhs_type==6)
      HYPRE_ParVectorDestroy(b);
   else
      HYPRE_IJVectorDestroy(ij_b);

   HYPRE_IJVectorDestroy(ij_x);

/*
  hypre_FinalizeMemoryDebug();
*/
  final:

   hypre_MPI_Finalize();

   return (0);
}

/*----------------------------------------------------------------------
 * Build standard 7-point laplacian in 3D with grid and anisotropy.
 * Parameters given in command line.
 *----------------------------------------------------------------------*/

HYPRE_Int
BuildParLaplacian( HYPRE_Int                  argc,
                   char                *argv[],
                   HYPRE_Int                  arg_index,
                   HYPRE_ParCSRMatrix  *A_ptr     )
{
   HYPRE_Int                 nx, ny, nz;
   HYPRE_Int                 gnx, gny, gnz;
   HYPRE_Int                 P, Q, R;
   HYPRE_Real          cx, cy, cz;

   HYPRE_ParCSRMatrix  A;

   HYPRE_Int                 num_procs, myid;
   HYPRE_Int                 p, q, r;
   HYPRE_Real         *values;

   
   /*-----------------------------------------------------------
    * Initialize some stuff
    *-----------------------------------------------------------*/

   hypre_MPI_Comm_size(hypre_MPI_COMM_WORLD, &num_procs );
   hypre_MPI_Comm_rank(hypre_MPI_COMM_WORLD, &myid );

   /*-----------------------------------------------------------
    * Set defaults
    *-----------------------------------------------------------*/
 
   nx = 30;
   ny = 30;
   nz = 30;

   P  = num_procs;
   Q  = 1;
   R  = 1;

   cx = 1.;
   cy = 1.;
   cz = 1.;

   /*-----------------------------------------------------------
    * Parse command line
    *-----------------------------------------------------------*/
   arg_index = 0;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-n") == 0 )
      {
         arg_index++;
         nx = atoi(argv[arg_index++]);
         ny = atoi(argv[arg_index++]);
         nz = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-P") == 0 )
      {
         arg_index++;
         P  = atoi(argv[arg_index++]);
         Q  = atoi(argv[arg_index++]);
         R  = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-c") == 0 )
      {
         arg_index++;
         cx = atof(argv[arg_index++]);
         cy = atof(argv[arg_index++]);
         cz = atof(argv[arg_index++]);
      }
      else
      {
         arg_index++;
      }
   }

   /*-----------------------------------------------------------
    * Check a few things
    *-----------------------------------------------------------*/

   if ((P*Q*R) != num_procs)
   {
      hypre_printf("Error: Invalid number of processors or processor topology \n");
      exit(1);
   }

   gnx = P*nx;
   gny = Q*ny;
   gnz = R*nz;

   /*-----------------------------------------------------------
    * Print driver parameters
    *-----------------------------------------------------------*/
 
   if (myid == 0)
   {
      hypre_printf("    (nx, ny, nz) = (%d, %d, %d)\n", nx, ny, nz);
      hypre_printf("    (Px, Py, Pz) = (%d, %d, %d)\n", P,  Q,  R);
      hypre_printf("    (cx, cy, cz) = (%f, %f, %f)\n\n", cx, cy, cz);
      hypre_printf("    Problem size = (%d x %d x %d)\n\n", gnx, gny, gnz);
   }

   /*-----------------------------------------------------------
    * Set up the grid structure
    *-----------------------------------------------------------*/

   /* compute p,q,r from P,Q,R and myid */
   p = myid % P;
   q = (( myid - p)/P) % Q;
   r = ( myid - p - P*q)/( P*Q );

   /*-----------------------------------------------------------
    * Generate the matrix 
    *-----------------------------------------------------------*/
 
   values = hypre_CTAlloc(HYPRE_Real,  4, HYPRE_MEMORY_HOST);

   values[1] = -cx;
   values[2] = -cy;
   values[3] = -cz;

   values[0] = 0.;
   if (nx > 1)
   {
      values[0] += 2.0*cx;
   }
   if (ny > 1)
   {
      values[0] += 2.0*cy;
   }
   if (nz > 1)
   {
      values[0] += 2.0*cz;
   }

   A = (HYPRE_ParCSRMatrix) GenerateLaplacian(hypre_MPI_COMM_WORLD, 
                                                 gnx, gny, gnz, P, Q, R, p, q, r, values);

   hypre_TFree(values, HYPRE_MEMORY_HOST);

   *A_ptr = A;

   return (0);
}

/*----------------------------------------------------------------------
 * returns the sign of a real number
 *  1 : positive
 *  0 : zero
 * -1 : negative
 *----------------------------------------------------------------------*/
static inline HYPRE_Int sign_double(HYPRE_Real a)
{
   return ( (0.0 < a) - (0.0 > a) );
}

/*----------------------------------------------------------------------
 * Build standard 7-point convection-diffusion operator 
 * Parameters given in command line.
 * Operator:
 *
 *  -cx Dxx - cy Dyy - cz Dzz + ax Dx + ay Dy + az Dz = f
 *
 *----------------------------------------------------------------------*/

HYPRE_Int
BuildParDifConv( HYPRE_Int                  argc,
                 char                *argv[],
                 HYPRE_Int                  arg_index,
                 HYPRE_ParCSRMatrix  *A_ptr)
{
   HYPRE_Int           nx, ny, nz;
   HYPRE_Int           gnx, gny, gnz;
   HYPRE_Int           P, Q, R;
   HYPRE_Real          cx, cy, cz;
   HYPRE_Real          ax, ay, az, atype;
   HYPRE_Real          hinx,hiny,hinz;
   HYPRE_Int           sign_prod;

   HYPRE_ParCSRMatrix  A;

   HYPRE_Int           num_procs, myid;
   HYPRE_Int           p, q, r;
   HYPRE_Real         *values;

   /*-----------------------------------------------------------
    * Initialize some stuff
    *-----------------------------------------------------------*/

   hypre_MPI_Comm_size(hypre_MPI_COMM_WORLD, &num_procs );
   hypre_MPI_Comm_rank(hypre_MPI_COMM_WORLD, &myid );

   /*-----------------------------------------------------------
    * Set defaults
    *-----------------------------------------------------------*/
 
   nx = 30;
   ny = 30;
   nz = 30;

   P  = num_procs;
   Q  = 1;
   R  = 1;

   cx = 1.;
   cy = 1.;
   cz = 1.;

   ax = 10.;
   ay = 10.;
   az = 10.;

   atype = 0;

   /*-----------------------------------------------------------
    * Parse command line
    *-----------------------------------------------------------*/
   arg_index = 0;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-n") == 0 )
      {
         arg_index++;
         nx = atoi(argv[arg_index++]);
         ny = atoi(argv[arg_index++]);
         nz = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-P") == 0 )
      {
         arg_index++;
         P  = atoi(argv[arg_index++]);
         Q  = atoi(argv[arg_index++]);
         R  = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-c") == 0 )
      {
         arg_index++;
         cx = atof(argv[arg_index++]);
         cy = atof(argv[arg_index++]);
         cz = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-a") == 0 )
      {
         arg_index++;
         ax = atof(argv[arg_index++]);
         ay = ax;
         az = ax;
      }
      else if ( strcmp(argv[arg_index], "-atype") == 0 )
      {
         arg_index++;
         atype = atoi(argv[arg_index++]);
      }
      else
      {
         arg_index++;
      }
   }

   /*-----------------------------------------------------------
    * Check a few things
    *-----------------------------------------------------------*/

   if ((P*Q*R) != num_procs)
   {
      hypre_printf("Error: Invalid number of processors or processor topology \n");
      exit(1);
   }

   gnx = P*nx;
   gny = Q*ny;
   gnz = R*nz;

   /*-----------------------------------------------------------
    * Print driver parameters
    *-----------------------------------------------------------*/
 
   if (myid == 0)
   {
      hypre_printf("  Convection-Diffusion: \n");
      hypre_printf("    -cx Dxx - cy Dyy - cz Dzz + ax Dx + ay Dy + az Dz = f\n");  
      hypre_printf("    (nx, ny, nz) = (%d, %d, %d)\n", nx, ny, nz);
      hypre_printf("    (Px, Py, Pz) = (%d, %d, %d)\n", P,  Q,  R);
      hypre_printf("    (cx, cy, cz) = (%f, %f, %f)\n", cx, cy, cz);
      hypre_printf("    (ax, ay, az) = (%f, %f, %f)\n\n", ax, ay, az);
      hypre_printf("    Problem size = (%f x %f x %f)\n\n", gnx, gny, gnz);
   }

   /*-----------------------------------------------------------
    * Set up the grid structure
    *-----------------------------------------------------------*/

   /* compute p,q,r from P,Q,R and myid */
   p = myid % P;
   q = (( myid - p)/P) % Q;
   r = ( myid - p - P*q)/( P*Q );

   hinx = 1./(gnx+1);
   hiny = 1./(gny+1);
   hinz = 1./(gnz+1);

   /*-----------------------------------------------------------
    * Generate the matrix 
    *-----------------------------------------------------------*/
   /* values[7]:
    *    [0]: center
    *    [1]: X-
    *    [2]: Y-
    *    [3]: Z-
    *    [4]: X+
    *    [5]: Y+
    *    [6]: Z+
    */    
   values = hypre_CTAlloc(HYPRE_Real,  7, HYPRE_MEMORY_HOST);

   values[0] = 0.;

   if (0 == atype) /* forward scheme for conv */
   {
      values[1] = -cx/(hinx*hinx);
      values[2] = -cy/(hiny*hiny);
      values[3] = -cz/(hinz*hinz);
      values[4] = -cx/(hinx*hinx) + ax/hinx;
      values[5] = -cy/(hiny*hiny) + ay/hiny;
      values[6] = -cz/(hinz*hinz) + az/hinz;

      if (nx > 1)
      {
         values[0] += 2.0*cx/(hinx*hinx) - 1.*ax/hinx;
      }
      if (ny > 1)
      {
         values[0] += 2.0*cy/(hiny*hiny) - 1.*ay/hiny;
      }
      if (nz > 1)
      {
         values[0] += 2.0*cz/(hinz*hinz) - 1.*az/hinz;
      }
   } 
   else if (1 == atype) /* backward scheme for conv */
   {
      values[1] = -cx/(hinx*hinx) - ax/hinx;
      values[2] = -cy/(hiny*hiny) - ay/hiny;
      values[3] = -cz/(hinz*hinz) - az/hinz;
      values[4] = -cx/(hinx*hinx);
      values[5] = -cy/(hiny*hiny);
      values[6] = -cz/(hinz*hinz);

      if (nx > 1)
      {
         values[0] += 2.0*cx/(hinx*hinx) + 1.*ax/hinx;
      }
      if (ny > 1)
      {
         values[0] += 2.0*cy/(hiny*hiny) + 1.*ay/hiny;
      }
      if (nz > 1)
      {
         values[0] += 2.0*cz/(hinz*hinz) + 1.*az/hinz;
      }
   }
   else if (3 == atype) /* upwind scheme */
   {
      sign_prod = sign_double(cx) * sign_double(ax);
      if (sign_prod == 1) /* same sign use back scheme */
      {
         values[1] = -cx/(hinx*hinx) - ax/hinx;
         values[4] = -cx/(hinx*hinx);
         if (nx > 1)
         {
            values[0] += 2.0*cx/(hinx*hinx) + 1.*ax/hinx;
         }
      }
      else /* diff sign use forward scheme */
      {
         values[1] = -cx/(hinx*hinx);
         values[4] = -cx/(hinx*hinx) + ax/hinx;
         if (nx > 1)
         {
            values[0] += 2.0*cx/(hinx*hinx) - 1.*ax/hinx;
         }
      }

      sign_prod = sign_double(cy) * sign_double(ay);
      if (sign_prod == 1) /* same sign use back scheme */
      {
         values[2] = -cy/(hiny*hiny) - ay/hiny;
         values[5] = -cy/(hiny*hiny);
         if (ny > 1)
         {
            values[0] += 2.0*cy/(hiny*hiny) + 1.*ay/hiny;
         }
      }
      else /* diff sign use forward scheme */
      {
         values[2] = -cy/(hiny*hiny);
         values[5] = -cy/(hiny*hiny) + ay/hiny;
         if (ny > 1)
         {
            values[0] += 2.0*cy/(hiny*hiny) - 1.*ay/hiny;
         }
      }

      sign_prod = sign_double(cz) * sign_double(az);
      if (sign_prod == 1) /* same sign use back scheme */
      {
         values[3] = -cz/(hinz*hinz) - az/hinz;
         values[6] = -cz/(hinz*hinz);
         if (nz > 1)
         {
            values[0] += 2.0*cz/(hinz*hinz) + 1.*az/hinz;
         }
      }
      else /* diff sign use forward scheme */
      {
         values[3] = -cz/(hinz*hinz);
         values[6] = -cz/(hinz*hinz) + az/hinz;
         if (nz > 1)
         {
            values[0] += 2.0*cz/(hinz*hinz) - 1.*az/hinz;
         }
      }
   }
   else /* centered difference scheme */
   {
      values[1] = -cx/(hinx*hinx) - ax/(2.*hinx);
      values[2] = -cy/(hiny*hiny) - ay/(2.*hiny);
      values[3] = -cz/(hinz*hinz) - az/(2.*hinz);
      values[4] = -cx/(hinx*hinx) + ax/(2.*hinx);
      values[5] = -cy/(hiny*hiny) + ay/(2.*hiny);
      values[6] = -cz/(hinz*hinz) + az/(2.*hinz);

      if (nx > 1)
      {
         values[0] += 2.0*cx/(hinx*hinx);
      }
      if (ny > 1)
      {
         values[0] += 2.0*cy/(hiny*hiny);
      }
      if (nz > 1)
      {
         values[0] += 2.0*cz/(hinz*hinz);
      }
   }

   A = (HYPRE_ParCSRMatrix) GenerateDifConv(hypre_MPI_COMM_WORLD,
                                            gnx, gny, gnz, P, Q, R, p, q, r, values);

   hypre_TFree(values, HYPRE_MEMORY_HOST);

   *A_ptr = A;

   return (0);
}

/*----------------------------------------------------------------------
 * Build matrix from one file on Proc. 0. Expects matrix to be in
 * CSR format. Distributes matrix across processors giving each about
 * the same number of rows.
 * Parameters given in command line.
 *----------------------------------------------------------------------*/

HYPRE_Int
BuildParFromOneFile( HYPRE_Int                  argc,
                     char                *argv[],
                     HYPRE_Int                  arg_index,
                     HYPRE_Int                  num_functions,
                     HYPRE_ParCSRMatrix  *A_ptr     )
{
   char               *filename;

   HYPRE_ParCSRMatrix  A;
   HYPRE_CSRMatrix  A_CSR = NULL;

   HYPRE_Int                 myid, numprocs;
   HYPRE_Int                 i, rest, size, num_nodes, num_dofs;
   HYPRE_Int		      *row_part;
   HYPRE_Int		      *col_part;

   /*-----------------------------------------------------------
    * Initialize some stuff
    *-----------------------------------------------------------*/

   hypre_MPI_Comm_rank(hypre_MPI_COMM_WORLD, &myid );
   hypre_MPI_Comm_size(hypre_MPI_COMM_WORLD, &numprocs );

   /*-----------------------------------------------------------
    * Parse command line
    *-----------------------------------------------------------*/

   if (arg_index < argc)
   {
      filename = argv[arg_index];
   }
   else
   {
      hypre_printf("Error: No filename specified \n");
      exit(1);
   }

   /*-----------------------------------------------------------
    * Print driver parameters
    *-----------------------------------------------------------*/
 
   if (myid == 0)
   {
      hypre_printf("  FromFile: %s\n", filename);

      /*-----------------------------------------------------------
       * Generate the matrix 
       *-----------------------------------------------------------*/
 
      A_CSR = HYPRE_CSRMatrixRead(filename);
   }

   row_part = NULL;
   col_part = NULL;
   if (myid == 0 && num_functions > 1)
   {
      HYPRE_CSRMatrixGetNumRows(A_CSR, &num_dofs);
      num_nodes = num_dofs/num_functions;
      if (num_dofs != num_functions*num_nodes)
      {
	 row_part = NULL;
	 col_part = NULL;
      }
      else
      {
         row_part = hypre_CTAlloc(HYPRE_Int,  numprocs+1, HYPRE_MEMORY_HOST);
	 row_part[0] = 0;
	 size = num_nodes/numprocs;
	 rest = num_nodes-size*numprocs;
	 for (i=0; i < numprocs; i++)
	 {
	    row_part[i+1] = row_part[i]+size*num_functions;
	    if (i < rest) row_part[i+1] += num_functions;
         }
         col_part = row_part;
      }
   }

   HYPRE_CSRMatrixToParCSRMatrix(hypre_MPI_COMM_WORLD, A_CSR, row_part, col_part, &A);

   *A_ptr = A;

   if (myid == 0) HYPRE_CSRMatrixDestroy(A_CSR);

   return (0);
}


/*----------------------------------------------------------------------
 * Build Rhs from one file on Proc. 0. Distributes vector across processors 
 * giving each about using the distribution of the matrix A.
 *----------------------------------------------------------------------*/

HYPRE_Int
BuildRhsParFromOneFile( HYPRE_Int                  argc,
                        char                *argv[],
                        HYPRE_Int                  arg_index,
                        HYPRE_Int                 *partitioning,
                        HYPRE_ParVector     *b_ptr     )
{
   char           *filename;

   HYPRE_ParVector b;
   HYPRE_Vector    b_CSR=NULL;

   HYPRE_Int             myid;

   /*-----------------------------------------------------------
    * Initialize some stuff
    *-----------------------------------------------------------*/

   hypre_MPI_Comm_rank(hypre_MPI_COMM_WORLD, &myid );

   /*-----------------------------------------------------------
    * Parse command line
    *-----------------------------------------------------------*/

   if (arg_index < argc)
   {
      filename = argv[arg_index];
   }
   else
   {
      hypre_printf("Error: No filename specified \n");
      exit(1);
   }

   /*-----------------------------------------------------------
    * Print driver parameters
    *-----------------------------------------------------------*/
 
   if (myid == 0)
   {
      hypre_printf("  Rhs FromFile: %s\n", filename);

      /*-----------------------------------------------------------
       * Generate the matrix 
       *-----------------------------------------------------------*/
 
      b_CSR = HYPRE_VectorRead(filename);
   }
   HYPRE_VectorToParVector(hypre_MPI_COMM_WORLD, b_CSR, partitioning,&b); 

   *b_ptr = b;

   HYPRE_VectorDestroy(b_CSR);

   return (0);
}

/*----------------------------------------------------------------------
 * Build standard 9-point laplacian in 2D with grid and anisotropy.
 * Parameters given in command line.
 *----------------------------------------------------------------------*/

HYPRE_Int
BuildParLaplacian9pt( HYPRE_Int                  argc,
                      char                *argv[],
                      HYPRE_Int                  arg_index,
                      HYPRE_ParCSRMatrix  *A_ptr     )
{
   HYPRE_Int                 nx, ny;
   HYPRE_Int                 gnx, gny;
   HYPRE_Int                 P, Q;

   HYPRE_ParCSRMatrix  A;

   HYPRE_Int                 num_procs, myid;
   HYPRE_Int                 p, q;
   HYPRE_Real         *values;

   /*-----------------------------------------------------------
    * Initialize some stuff
    *-----------------------------------------------------------*/

   hypre_MPI_Comm_size(hypre_MPI_COMM_WORLD, &num_procs );
   hypre_MPI_Comm_rank(hypre_MPI_COMM_WORLD, &myid );

   /*-----------------------------------------------------------
    * Set defaults
    *-----------------------------------------------------------*/
 
   nx = 10;
   ny = 10;

   P  = num_procs;
   Q  = 1;

   /*-----------------------------------------------------------
    * Parse command line
    *-----------------------------------------------------------*/
   arg_index = 0;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-n") == 0 )
      {
         arg_index++;
         nx = atoi(argv[arg_index++]);
         ny = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-P") == 0 )
      {
         arg_index++;
         P  = atoi(argv[arg_index++]);
         Q  = atoi(argv[arg_index++]);
      }
      else
      {
         arg_index++;
      }
   }

   /*-----------------------------------------------------------
    * Check a few things
    *-----------------------------------------------------------*/

   if ((P*Q) != num_procs)
   {
      hypre_printf("Error: Invalid number of processors or processor topology \n");
      exit(1);
   }

   gnx = P*nx;
   gny = Q*ny;

   /*-----------------------------------------------------------
    * Print driver parameters
    *-----------------------------------------------------------*/
 
   if (myid == 0)
   {
      hypre_printf("  Laplacian 9pt:\n");
      hypre_printf("    (nx, ny) = (%d, %d)\n", nx, ny);
      hypre_printf("    (Px, Py) = (%d, %d)\n\n", P,  Q);
      hypre_printf("problem size = (%d x %d)\n", gnx, gny);
   }

   /*-----------------------------------------------------------
    * Set up the grid structure
    *-----------------------------------------------------------*/

   /* compute p,q from P,Q and myid */
   p = myid % P;
   q = ( myid - p)/P;

   /*-----------------------------------------------------------
    * Generate the matrix 
    *-----------------------------------------------------------*/
 
   values = hypre_CTAlloc(HYPRE_Real,  2, HYPRE_MEMORY_HOST);

   values[1] = -1.;

   values[0] = 0.;
   if (nx > 1)
   {
      values[0] += 2.0;
   }
   if (ny > 1)
   {
      values[0] += 2.0;
   }
   if (nx > 1 && ny > 1)
   {
      values[0] += 4.0;
   }

   A = (HYPRE_ParCSRMatrix) GenerateLaplacian9pt(hypre_MPI_COMM_WORLD,
                                                 gnx, gny, P, Q, p, q, values);

   hypre_TFree(values, HYPRE_MEMORY_HOST);

   *A_ptr = A;

   return (0);
}
/*----------------------------------------------------------------------
 * Build 27-point laplacian in 3D,
 * Parameters given in command line.
 *----------------------------------------------------------------------*/

HYPRE_Int
BuildParLaplacian27pt( HYPRE_Int                  argc,
                       char                *argv[],
                       HYPRE_Int                  arg_index,
                       HYPRE_ParCSRMatrix  *A_ptr     )
{
   HYPRE_Int                 nx, ny, nz;
   HYPRE_Int                 gnx, gny, gnz;
   HYPRE_Int                 P, Q, R;

   HYPRE_ParCSRMatrix  A;

   HYPRE_Int                 num_procs, myid;
   HYPRE_Int                 p, q, r;
   HYPRE_Real         *values;

   /*-----------------------------------------------------------
    * Initialize some stuff
    *-----------------------------------------------------------*/

   hypre_MPI_Comm_size(hypre_MPI_COMM_WORLD, &num_procs );
   hypre_MPI_Comm_rank(hypre_MPI_COMM_WORLD, &myid );

   /*-----------------------------------------------------------
    * Set defaults
    *-----------------------------------------------------------*/
 
   nx = 10;
   ny = 10;
   nz = 10;

   P  = num_procs;;
   Q  = 1;
   R  = 1;

   /*-----------------------------------------------------------
    * Parse command line
    *-----------------------------------------------------------*/
   arg_index = 0;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-n") == 0 )
      {
         arg_index++;
         nx = atoi(argv[arg_index++]);
         ny = atoi(argv[arg_index++]);
         nz = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-P") == 0 )
      {
         arg_index++;
         P  = atoi(argv[arg_index++]);
         Q  = atoi(argv[arg_index++]);
         R  = atoi(argv[arg_index++]);
      }
      else
      {
         arg_index++;
      }
   }

   gnx = P*nx;
   gny = Q*ny;
   gnz = R*nz;

   /*-----------------------------------------------------------
    * Check a few things
    *-----------------------------------------------------------*/

   if ((P*Q*R) != num_procs)
   {
      hypre_printf("Error: Invalid number of processors or processor topology \n");
      exit(1);
   }

   /*-----------------------------------------------------------
    * Print driver parameters
    *-----------------------------------------------------------*/
 
   if (myid == 0)
   {
      hypre_printf("  Laplacian_27pt:\n");
      hypre_printf("    (nx, ny, nz) = (%d, %d, %d)\n", nx, ny, nz);
      hypre_printf("    (Px, Py, Pz) = (%d, %d, %d)\n\n", P,  Q,  R);
      hypre_printf("    problem size = (%d x %d x %d)\n", gnx, gny, gnz);
   }

   /*-----------------------------------------------------------
    * Set up the grid structure
    *-----------------------------------------------------------*/

   /* compute p,q,r from P,Q,R and myid */
   p = myid % P;
   q = (( myid - p)/P) % Q;
   r = ( myid - p - P*q)/( P*Q );

   /*-----------------------------------------------------------
    * Generate the matrix 
    *-----------------------------------------------------------*/
 
   values = hypre_CTAlloc(HYPRE_Real,  2, HYPRE_MEMORY_HOST);

   values[0] = 26.0;
   if (nx == 1 || ny == 1 || nz == 1)
      values[0] = 8.0;
   if (nx*ny == 1 || nx*nz == 1 || ny*nz == 1)
      values[0] = 2.0;
   values[1] = -1.;

   A = (HYPRE_ParCSRMatrix) GenerateLaplacian27pt(hypre_MPI_COMM_WORLD,
                                                  gnx, gny, gnz, P, Q, R, p, q, r, values);

   hypre_TFree(values, HYPRE_MEMORY_HOST);

   *A_ptr = A;

   return (0);
}


/*----------------------------------------------------------------------
 * Build 7-point in 2D 
 * Parameters given in command line.
 *----------------------------------------------------------------------*/

HYPRE_Int
BuildParRotate7pt( HYPRE_Int                  argc,
                   char                *argv[],
                   HYPRE_Int                  arg_index,
                   HYPRE_ParCSRMatrix  *A_ptr     )
{
   HYPRE_Int                 nx, ny;
   HYPRE_Int                 gnx, gny;
   HYPRE_Int                 P, Q;

   HYPRE_ParCSRMatrix  A;

   HYPRE_Int                 num_procs, myid;
   HYPRE_Int                 p, q;
   HYPRE_Real          eps, alpha;

   /*-----------------------------------------------------------
    * Initialize some stuff
    *-----------------------------------------------------------*/

   hypre_MPI_Comm_size(hypre_MPI_COMM_WORLD, &num_procs );
   hypre_MPI_Comm_rank(hypre_MPI_COMM_WORLD, &myid );

   /*-----------------------------------------------------------
    * Set defaults
    *-----------------------------------------------------------*/

   nx = 10;
   ny = 10;

   P  = num_procs;
   Q  = 1;

   /*-----------------------------------------------------------
    * Parse command line
    *-----------------------------------------------------------*/
   arg_index = 0;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-n") == 0 )
      {
         arg_index++;
         nx = atoi(argv[arg_index++]);
         ny = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-P") == 0 )
      {
         arg_index++;
         P  = atoi(argv[arg_index++]);
         Q  = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-alpha") == 0 )
      {
         arg_index++;
         alpha  = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-eps") == 0 )
      {
         arg_index++;
         eps  = atof(argv[arg_index++]);
      }
      else
      {
         arg_index++;
      }
   }

   /*-----------------------------------------------------------
    * Check a few things
    *-----------------------------------------------------------*/

   if ((P*Q) != num_procs)
   {
      hypre_printf("Error: Invalid number of processors or processor topology \n");
      exit(1);
   }

   gnx = P*nx;
   gny = Q*ny;

   /*-----------------------------------------------------------
    * Print driver parameters
    *-----------------------------------------------------------*/

   if (myid == 0)
   {
      hypre_printf("  Rotate 7pt:\n");
      hypre_printf("    alpha = %f, eps = %f\n", alpha,eps);
      hypre_printf("    (nx, ny) = (%d, %d)\n", nx, ny);
      hypre_printf("    (Px, Py) = (%d, %d)\n", P,  Q);
      hypre_printf("problem size = (%d x %d)\n", gnx, gny);
   }

   /*-----------------------------------------------------------
    * Set up the grid structure
    *-----------------------------------------------------------*/

   /* compute p,q from P,Q and myid */
   p = myid % P;
   q = ( myid - p)/P;

   /*-----------------------------------------------------------
    * Generate the matrix 
    *-----------------------------------------------------------*/

   A = (HYPRE_ParCSRMatrix) GenerateRotate7pt(hypre_MPI_COMM_WORLD,
                                              gnx, gny, P, Q, p, q, alpha, eps);

   *A_ptr = A;

   return (0);
}

/*----------------------------------------------------------------------
 * Build standard 7-point difference operator using centered differences
 *
 *  eps*(a(x,y,z) ux)x + (b(x,y,z) uy)y + (c(x,y,z) uz)z 
 *  d(x,y,z) ux + e(x,y,z) uy + f(x,y,z) uz + g(x,y,z) u
 *
 *  functions a,b,c,d,e,f,g need to be defined inside par_vardifconv.c
 *
 *----------------------------------------------------------------------*/

HYPRE_Int
BuildParVarDifConv( HYPRE_Int                  argc,
                    char                *argv[],
                    HYPRE_Int                  arg_index,
                    HYPRE_ParCSRMatrix  *A_ptr    ,
                    HYPRE_ParVector  *rhs_ptr     )
{
   HYPRE_Int                 nx, ny, nz;
   HYPRE_Int                 gnx, gny, gnz;
   HYPRE_Int                 P, Q, R;

   HYPRE_ParCSRMatrix  A;
   HYPRE_ParVector  rhs;

   HYPRE_Int           num_procs, myid;
   HYPRE_Int           p, q, r;
   HYPRE_Int           type;
   HYPRE_Real          eps;

   /*-----------------------------------------------------------
    * Initialize some stuff
    *-----------------------------------------------------------*/

   hypre_MPI_Comm_size(hypre_MPI_COMM_WORLD, &num_procs );
   hypre_MPI_Comm_rank(hypre_MPI_COMM_WORLD, &myid );

   /*-----------------------------------------------------------
    * Set defaults
    *-----------------------------------------------------------*/

   nx = 10;
   ny = 10;
   nz = 10;
   P  = num_procs;
   Q  = 1;
   R  = 1;
   eps = 1.0;

   /* type: 0   : default FD;
    *       1-3 : FD and examples 1-3 in Ruge-Stuben paper */
   type = 0;

   /*-----------------------------------------------------------
    * Parse command line
    *-----------------------------------------------------------*/
   arg_index = 0;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-n") == 0 )
      {
         arg_index++;
         nx = atoi(argv[arg_index++]);
         ny = atoi(argv[arg_index++]);
         nz = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-P") == 0 )
      {
         arg_index++;
         P  = atoi(argv[arg_index++]);
         Q  = atoi(argv[arg_index++]);
         R  = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-eps") == 0 )
      {
         arg_index++;
         eps  = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-vardifconvRS") == 0 )
      {
         arg_index++;
         type = atoi(argv[arg_index++]);
      }
      else
      {
         arg_index++;
      }
   }

   /*-----------------------------------------------------------
    * Check a few things
    *-----------------------------------------------------------*/

   if ((P*Q*R) != num_procs)
   {
      hypre_printf("Error: Invalid number of processors or processor topology \n");
      exit(1);
   }

   gnx = P*nx;
   gny = Q*ny;
   gnz = R*nz;

   /*-----------------------------------------------------------
    * Print driver parameters
    *-----------------------------------------------------------*/

   if (myid == 0)
   {
      hypre_printf("  ell PDE: eps = %f\n", eps);
      hypre_printf("    Dx(aDxu) + Dy(bDyu) + Dz(cDzu) + d Dxu + e Dyu + f Dzu  + g u= f\n");
      hypre_printf("    (nx, ny, nz) = (%d, %d, %d)\n", nx, ny, nz);
      hypre_printf("    (Px, Py, Pz) = (%d, %d, %d)\n", P,  Q,  R);
      hypre_printf("    Problem size = (%d x %d x %d)\n", gnx, gny, gnz);
   }
   /*-----------------------------------------------------------
    * Set up the grid structure
    *-----------------------------------------------------------*/

   /* compute p,q,r from P,Q,R and myid */
   p = myid % P;
   q = (( myid - p)/P) % Q;
   r = ( myid - p - P*q)/( P*Q );

   /*-----------------------------------------------------------
    * Generate the matrix
    *-----------------------------------------------------------*/

   if (0 == type)
   {
      A = (HYPRE_ParCSRMatrix) GenerateVarDifConv(hypre_MPI_COMM_WORLD,
                                                  gnx, gny, gnz, P, Q, R, p, q, r, eps, &rhs);
   }
   else
   {
      A = (HYPRE_ParCSRMatrix) GenerateRSVarDifConv(hypre_MPI_COMM_WORLD,
                                                    gnx, gny, gnz, P, Q, R, p, q, r, eps, &rhs,
                                                    type);
   }

   *A_ptr = A;
   *rhs_ptr = rhs;

   return (0);
}

HYPRE_Int
hypre_PrintTiming2( const char     *heading,
                   HYPRE_Real     *wall_time_ptr,
                   MPI_Comm        comm  )
{
   HYPRE_Int  ierr = 0;

   HYPRE_Real  local_wall_time;
   HYPRE_Real  local_cpu_time;
   HYPRE_Real  wall_time;
   HYPRE_Real  cpu_time;
   HYPRE_Real  wall_mflops;
   HYPRE_Real  cpu_mflops;

   HYPRE_Int     i;
   HYPRE_Int     myrank;

   if (hypre_global_timing == NULL)
      return ierr;

   hypre_MPI_Comm_rank(comm, &myrank );

   /* print heading */
   if (myrank == 0)
   {
      hypre_printf("=============================================\n");
      hypre_printf("%s:\n", heading);
      hypre_printf("=============================================\n");
   }

   for (i = 0; i < (hypre_global_timing -> size); i++)
   {
      if (hypre_TimingNumRegs(i) > 0)
      {
         local_wall_time = hypre_TimingWallTime(i);
         local_cpu_time  = hypre_TimingCPUTime(i);
         hypre_MPI_Allreduce(&local_wall_time, &wall_time, 1,
                       hypre_MPI_REAL, hypre_MPI_MAX, comm);
         hypre_MPI_Allreduce(&local_cpu_time, &cpu_time, 1,
                       hypre_MPI_REAL, hypre_MPI_MAX, comm);

         if (myrank == 0)
         {
            hypre_printf("%s:\n", hypre_TimingName(i));
            *wall_time_ptr = wall_time;

            /* print wall clock info */
            hypre_printf("  wall clock time = %f seconds\n", wall_time);
            if (wall_time)
               wall_mflops = hypre_TimingFLOPS(i) / wall_time / 1.0E6;
            else
               wall_mflops = 0.0;
            hypre_printf("  wall MFLOPS     = %f\n", wall_mflops);

            /* print CPU clock info */
            hypre_printf("  cpu clock time  = %f seconds\n", cpu_time);
            if (cpu_time)
               cpu_mflops = hypre_TimingFLOPS(i) / cpu_time / 1.0E6;
            else
               cpu_mflops = 0.0;
            hypre_printf("  cpu MFLOPS      = %f\n\n", cpu_mflops);
         }
      }
   }

   return ierr;
}



