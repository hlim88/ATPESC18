//                     MFEM DG Advection Example with PETSc
//
// Compile with: make advection-ode
//
// Sample runs:
//    mpirun -np 4 advection-ode -m ../../data/periodic-hexagon.mesh
//    mpirun -np 4 advection-ode -m ../../data/periodic-hexagon.mesh -implicit
//
// Description:  This example code solves the time-dependent advection equation
//               du/dt + v.grad(u) = 0, where v is a given fluid velocity, and
//               u0(x)=u(0,x) is a given initial condition.
//
//               The example demonstrates the use of Discontinuous Galerkin (DG)
//               bilinear forms in MFEM (face integrators), the use of explicit
//               ODE time integrators, the definition of periodic boundary
//               conditions through periodic meshes, as well as the use of GLVis
//               for persistent visualization of a time-evolving solution. The
//               saving of time-dependent data files for external visualization
//               with VisIt (visit.llnl.gov) is also illustrated.
//
//               The example also demonstrates how to use PETSc ODE solvers and
//               customize them by command line (see rc_ex9p_expl and
//               rc_ex9p_impl). The split in left-hand side and right-hand side
//               of the TimeDependentOperator is amenable for IMEX methods.
//               When using fully implicit methods, just the left-hand side of
//               the operator should be provided for efficiency reasons when
//               assembling the Jacobians. Here, we provide two Jacobian
//               routines just to illustrate the capabilities of the
//               PetscODESolver class.  We also show how to monitor the time
//               dependent solution inside a call to PetscODESolver:Mult.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v);

// Initial condition
double u0_function(const Vector &x);

// Inflow boundary condition
double inflow_function(const Vector &x);

// Mesh bounding box
Vector bb_min, bb_max;


/** A time-dependent operator for the ODE as F(u,du/dt,t) = G(u,t)
    The DG weak form of du/dt = -v.grad(u) is M du/dt = K u + b, where M and K are the mass
    and advection matrices, and b describes the flow on the boundary. This can
    be also written as a general ODE with the right-hand side only as
    du/dt = M^{-1} (K u + b).
    This class is used to evaluate the right-hand side and the left-hand side. */
class FE_Evolution : public TimeDependentOperator
{
private:
   HypreParMatrix &M, &K;
   const Vector &b;
   HypreSmoother M_prec;
   CGSolver M_solver;

   mutable Vector z;
   mutable PetscParMatrix* iJacobian;
   mutable PetscParMatrix* rJacobian;

public:
   FE_Evolution(HypreParMatrix &_M, HypreParMatrix &_K, const Vector &_b,
                bool M_in_lhs);

   virtual void ExplicitMult(const Vector &x, Vector &y) const;
   virtual void ImplicitMult(const Vector &x, const Vector &xp, Vector &y) const;
   virtual void Mult(const Vector &x, Vector &y) const;
   virtual Operator& GetExplicitGradient(const Vector &x) const;
   virtual Operator& GetImplicitGradient(const Vector &x, const Vector &xp,
                                         double shift) const;
   virtual ~FE_Evolution() { delete iJacobian; delete rJacobian; }
};


// Monitor the solution at time step "step", explicitly in the time loop
class UserMonitor : public PetscSolverMonitor
{
private:
   socketstream&    sout;
   ParMesh*         pmesh;
   ParGridFunction* u;
   int              vt;
   bool             pause;

public:
   UserMonitor(socketstream& _s, ParMesh* _m, ParGridFunction* _u, int _vt) :
      PetscSolverMonitor(true,false), sout(_s), pmesh(_m), u(_u), vt(_vt),
      pause(true) {}

   void MonitorSolution(PetscInt step, PetscReal norm, const Vector &X)
   {
      if (step % vt == 0)
      {
         int  num_procs, myid;

         *u = X;
         MPI_Comm_size(pmesh->GetComm(),&num_procs);
         MPI_Comm_rank(pmesh->GetComm(),&myid);
         sout << "parallel " << num_procs << " " << myid << "\n";
         sout << "solution\n" << *pmesh << *u;
         if (pause) { sout << "pause\n"; }
         sout << flush;
         if (pause)
         {
            pause = false;
            if (myid == 0)
            {
               cout << "GLVis visualization paused."
                    << " Press space (in the GLVis window) to resume it.\n";
            }
         }
      }
   }
};


int main(int argc, char *argv[])
{
   // 1. Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   // 2. Parse command-line options.
   const char *mesh_file = "periodic-hexagon.mesh";
   int ser_ref_levels = 3;
   int par_ref_levels = 0;
   int order = 2;
   double t_final = 3.0;
   double dt = 0.01;
   bool visualization = true;
   bool visit = true;
   int vis_steps = 10;
   bool implicit = false;
   bool use_step = false;

   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                  "--no-visit-datafiles",
                  "Save data files for VisIt (visit.llnl.gov) visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.AddOption(&use_step, "-usestep", "--usestep", "-no-step",
                  "--no-step",
                  "Use the Step() or Run() method to solve the ODE system.");
   args.AddOption(&implicit, "-implicit", "--implicit", "-no-implicit",
                  "--no-implicit",
                  "Use or not an implicit method in PETSc to solve the ODE system.");
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(cout);
      }
      MPI_Finalize();
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }
   const char *petscrc_expl = "rc-advection-expl";
   const char *petscrc_impl = "rc-advection-impl";
   const char *petscrc_file = (implicit) ? petscrc_impl : petscrc_expl;

   // 3. Read the serial mesh from the given mesh file on all processors. We can
   //    handle geometrically periodic meshes in this code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 4. Define the ODE solver used for time integration. Several explicit
   //    Runge-Kutta methods are available.
   ODESolver *ode_solver = NULL;
   PetscODESolver *pode_solver = NULL;
   UserMonitor *pmon = NULL;
   PetscInitialize(NULL, NULL, petscrc_file, NULL);
   ode_solver = pode_solver = new PetscODESolver(MPI_COMM_WORLD);


   // 5. Refine the mesh in serial to increase the resolution. In this example
   //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
   //    a command-line parameter. If the mesh is of NURBS type, we convert it
   //    to a (piecewise-polynomial) high-order mesh.
   for (int lev = 0; lev < ser_ref_levels; lev++)
   {
      mesh->UniformRefinement();
   }
   mesh->GetBoundingBox(bb_min, bb_max, max(order, 1));

   // 6. Define the parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   for (int lev = 0; lev < par_ref_levels; lev++)
   {
      pmesh->UniformRefinement();
   }

   // 7. Define the parallel discontinuous DG finite element space on the
   //    parallel refined mesh of the given polynomial order.
   DG_FECollection fec(order, dim);
   ParFiniteElementSpace *fes = new ParFiniteElementSpace(pmesh, &fec);

   HYPRE_Int global_vSize = fes->GlobalTrueVSize();
   if (myid == 0)
   {
      cout << "Number of unknowns: " << global_vSize << endl;
   }

   // 8. Set up and assemble the parallel bilinear and linear forms (and the
   //    parallel hypre matrices) corresponding to the DG discretization. The
   //    DGTraceIntegrator involves integrals over mesh interior faces.
   VectorFunctionCoefficient velocity(dim, velocity_function);
   FunctionCoefficient inflow(inflow_function);
   FunctionCoefficient u0(u0_function);

   ParBilinearForm *m = new ParBilinearForm(fes);
   m->AddDomainIntegrator(new MassIntegrator);
   ParBilinearForm *k = new ParBilinearForm(fes);
   k->AddDomainIntegrator(new ConvectionIntegrator(velocity, -1.0));
   k->AddInteriorFaceIntegrator(
      new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));
   k->AddBdrFaceIntegrator(
      new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));

   ParLinearForm *b = new ParLinearForm(fes);
   b->AddBdrFaceIntegrator(
      new BoundaryFlowIntegrator(inflow, velocity, -1.0, -0.5));

   tic();
   m->Assemble();
   m->Finalize();
   int skip_zeros = 0;
   k->Assemble(skip_zeros);
   k->Finalize(skip_zeros);
   b->Assemble();

   HypreParMatrix *M = m->ParallelAssemble();
   HypreParMatrix *K = k->ParallelAssemble();
   HypreParVector *B = b->ParallelAssemble();
   cout << "Time for Matrix and vector assembly (s):  " << toc() << endl;

   // 9. Define the initial conditions, save the corresponding grid function to
   //    a file and (optionally) save data in the VisIt format and initialize
   //    GLVis visualization.
   ParGridFunction *u = new ParGridFunction(fes);
   u->ProjectCoefficient(u0);
   HypreParVector *U = u->GetTrueDofs();
   {
      ostringstream mesh_name, sol_name;
      mesh_name << "adv-mesh." << setfill('0') << setw(6) << myid;
      sol_name << "adv-init." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
      ofstream osol(sol_name.str().c_str());
      osol.precision(precision);
      u->Save(osol);
   }

   // Create data collection for solution output VisItDataCollection
   DataCollection *dc = NULL;
   if (visit)
   {
      dc = new VisItDataCollection("adv", pmesh);
      dc->SetPrecision(precision);
      dc->RegisterField("solution", u);
      dc->SetCycle(0);
      dc->SetTime(0.0);
      dc->Save();
   }

   socketstream sout;
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      sout.open(vishost, visport);

      // Set the monitoring routine for the PetscODESolver.
      sout.precision(precision);
      pmon = new UserMonitor(sout,pmesh,u,vis_steps);
      pode_solver->SetMonitor(pmon);
   }

   // 10. Define the time-dependent evolution operator describing the ODE
   FE_Evolution *adv = new FE_Evolution(*M, *K, *B, implicit);

   double t = 0.0;
   adv->SetTime(t);
   pode_solver->Init(*adv,PetscODESolver::ODE_SOLVER_LINEAR);


   // Explicitly perform time-integration (looping over the time iterations, ti,
   // with a time-step dt), or use the Run method of the ODE solver class.
   tic();
   if (use_step)
   {
      bool done = false;
      for (int ti = 0; !done; )
      {
         double dt_real = min(dt, t_final - t);
         ode_solver->Step(*U, t, dt_real);
         ti++;

         done = (t >= t_final - 1e-8*dt);

         if (done || ti % vis_steps == 0)
         {
            if (myid == 0)
            {
               cout << "time step: " << ti << ", time: " << t << endl;
            }
            // 11. Extract the parallel grid function corresponding to the finite
            //     element approximation U (the local solution on each processor).
            *u = *U;

            if (visualization)
            {
               sout << "parallel " << num_procs << " " << myid << "\n";
               sout << "solution\n" << *pmesh << *u << flush;
            }

            if (visit)
            {
               dc->SetCycle(ti);
               dc->SetTime(t);
               dc->Save();
            }
         }
      }
   }
   else { ode_solver->Run(*U, t, dt, t_final); }
   cout << "Time for all ODE steps (s):  " << toc() << endl;

   // 12. Save the final solution in parallel. This output can be viewed later
   {
      *u = *U;
      ostringstream sol_name;
      sol_name << "adv-final." << setfill('0') << setw(6) << myid;
      ofstream osol(sol_name.str().c_str());
      osol.precision(precision);
      u->Save(osol);
   }

   FunctionCoefficient u_exact_coef(u0_function);
   double l2_err = u->ComputeL2Error(u_exact_coef);
   cout << "L2 Error after trip around hexagon:  " << l2_err << endl;

   // 13. Free the used memory.
   delete U;
   delete u;
   delete B;
   delete b;
   delete K;
   delete k;
   delete M;
   delete m;
   delete fes;
   delete pmesh;
   delete ode_solver;
   delete dc;
   delete adv;

   delete pmon;

   // We finalize PETSc
   PetscFinalize();

   MPI_Finalize();
   return 0;
}


// Implementation of class FE_Evolution
FE_Evolution::FE_Evolution(HypreParMatrix &_M, HypreParMatrix &_K,
                           const Vector &_b,bool M_in_lhs)
   : TimeDependentOperator(_M.Height(), 0.0,
                           M_in_lhs ? TimeDependentOperator::IMPLICIT
                           : TimeDependentOperator::EXPLICIT),
     M(_M), K(_K), b(_b), M_solver(M.GetComm()), z(_M.Height()),
     iJacobian(NULL), rJacobian(NULL)
{
   if (isExplicit())
   {
      M_prec.SetType(HypreSmoother::Jacobi);
      M_solver.SetPreconditioner(M_prec);
      M_solver.SetOperator(M);

      M_solver.iterative_mode = false;
      M_solver.SetRelTol(1e-9);
      M_solver.SetAbsTol(0.0);
      M_solver.SetMaxIter(100);
      M_solver.SetPrintLevel(0);
   }
}

// RHS evaluation
void FE_Evolution::ExplicitMult(const Vector &x, Vector &y) const
{
   if (isExplicit())
   {
      // y = M^{-1} (K x + b)
      K.Mult(x, z);
      z += b;
      M_solver.Mult(z, y);
   }
   else
   {
      // y = K x + b
      K.Mult(x, y);
      y += b;
   }
}

// LHS evaluation
void FE_Evolution::ImplicitMult(const Vector &x, const Vector &xp,
                                Vector &y) const
{
   if (isImplicit())
   {
      M.Mult(xp, y);
   }
   else
   {
      y = xp;
   }
}

void FE_Evolution::Mult(const Vector &x, Vector &y) const
{
   // y = M^{-1} (K x + b)
   K.Mult(x, z);
   z += b;
   M_solver.Mult(z, y);
}

// RHS Jacobian
Operator& FE_Evolution::GetExplicitGradient(const Vector &x) const
{
   delete rJacobian;
   if (isImplicit())
   {
      rJacobian = new PetscParMatrix(&K, Operator::PETSC_MATAIJ);
   }
   else
   {
      mfem_error("FE_Evolution::GetExplicitGradient(x): Capability not coded!");
   }
   return *rJacobian;
}

// LHS Jacobian, evaluated as shift*F_du/dt + F_u
Operator& FE_Evolution::GetImplicitGradient(const Vector &x, const Vector &xp,
                                            double shift) const
{
   delete iJacobian;
   if (isImplicit())
   {
      iJacobian = new PetscParMatrix(&M, Operator::PETSC_MATAIJ);
      *iJacobian *= shift;
   }
   else
   {
      mfem_error("FE_Evolution::GetImplicitGradient(x,xp,shift):"
                 " Capability not coded!");
   }
   return *iJacobian;
}

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v)
{
   int dim = x.Size();

   // Translations in 1D, 2D, and 3D
   switch (dim)
   {
      case 1: v(0) = 1.0; break;
      case 2: v(0) = 1.0; v(1) = 0.0; break;
      case 3: v(0) = 1.0; v(1) = 0.0; v(2) = 0.0; break;
   }
}

// Initial condition
double u0_function(const Vector &x)
{
   int dim = x.Size();

   // map to the reference [-1,1] domain
   Vector X(dim);
   for (int i = 0; i < dim; i++)
   {
      double center = (bb_min[i] + bb_max[i]) * 0.5;
      X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
   }

   switch (dim)
   {
      case 1:
         return exp(-40.*pow(X(0)-0.5,2));
      case 2:
      case 3:
      {
         double rx = 0.45, ry = 0.25, cx = 0., cy = -0.2, w = 10.;
         if (dim == 3)
         {
            const double s = (1. + 0.25*cos(2*M_PI*X(2)));
            rx *= s;
            ry *= s;
         }
         return ( erfc(w*(X(0)-cx-rx))*erfc(-w*(X(0)-cx+rx)) *
                  erfc(w*(X(1)-cy-ry))*erfc(-w*(X(1)-cy+ry)) )/16;
      }
   }

   return 0.0;
}

// Inflow boundary condition (zero for the problems considered in this example)
double inflow_function(const Vector &x)
{
   return 0.0;
}
