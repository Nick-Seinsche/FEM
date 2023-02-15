#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   // Path to the mesh
   const char *mesh_file = "../Libraries/mfem-4.5/data/ref-square.mesh";
   //const char *mesh_file = "unit-square.mesh";
   int order = 1;

   // 2. Read the mesh from the given mesh file, and refine uniformly.
   Mesh mesh(mesh_file);
   mesh.UniformRefinement();
   mesh.UniformRefinement();
   mesh.UniformRefinement();
   mesh.UniformRefinement();
   mesh.UniformRefinement();
   mesh.UniformRefinement();
   mesh.UniformRefinement();

   // 3. Define a finite element space on the mesh. Here we use P1 continuous
   H1_FECollection fec(order, mesh.Dimension());
   FiniteElementSpace fespace(&mesh, &fec);

   cout << "Number of unknowns: " << fespace.GetTrueVSize() << endl;

   // 4. Extract the list of all the boundary DOFs. These will be marked as
   //    Dirichlet in order to enforce zero boundary conditions.
   Array<int> boundary_dofs;
   fespace.GetBoundaryTrueDofs(boundary_dofs);

   // 5. Define the solution x as a finite element grid function in fespace. Set
   //    the initial guess to zero, which also sets the boundary conditions.
   GridFunction x_eps(&fespace);
   x_eps = 0.0;

   GridFunction x_hom(&fespace);
   x_hom = 0.0;

   // 6. Set up the linear form b(.) corresponding to the right-hand side.
   std::function<double (const mfem::Vector &)> f = [](const Vector& x){return sin(4 * M_PI * x[0]) * cos(4 * M_PI * x[1]);};
   FunctionCoefficient fc(f);

   ConstantCoefficient one(1.0);
   LinearForm b(&fespace);
   //b.AddDomainIntegrator(new DomainLFIntegrator(one));
   b.AddDomainIntegrator(new DomainLFIntegrator(fc));
   b.Assemble();

   // 7. Set up the bilinear form a(.,.) corresponding to left hand side of the (weak formulation of the) pde
   double epsilon = 0.01;
   double low = 0.25, high = 4;

   std::function<void(const Vector &, DenseMatrix &)> mf = [epsilon, low, high](const Vector& v, DenseMatrix& A) {
      if ((int) floor(2 * v[0] / epsilon) % 2 == 0) {
         A(0,0) = low;
         A(1,0) = 0;
         A(0,1) = 0;
         A(1,1) = low;
      } else {
         A(0,0) = high;
         A(1,0) = 0;
         A(0,1) = 0;
         A(1,1) = high;
      }
   };
   int A_dim = 2;
   MatrixFunctionCoefficient mfc(A_dim, mf);

   // Define the homogenized matrix
   double homogenzied_array[2][2] = {{0.5 * high * low / (high + low), 0}, {0, 0.5 * high + 0.5 * low}};
   DenseMatrix hom_matrix(homogenzied_array);
   MatrixConstantCoefficient mc(hom_matrix);

   BilinearForm a_eps(&fespace);
   a_eps.AddDomainIntegrator(new DiffusionIntegrator(mfc));
   a_eps.Assemble();

   BilinearForm a_hom(&fespace);
   a_hom.AddDomainIntegrator(new DiffusionIntegrator(mc));
   a_hom.Assemble();

   // 8. Form the linear system A X = B. This includes eliminating boundary
   //    conditions, applying AMR constraints, and other transformations.
   SparseMatrix A_eps, A_hom;
   Vector B_eps, B_hom, X_eps, X_hom;
   a_eps.FormLinearSystem(boundary_dofs, x_eps, b, A_eps, X_eps, B_eps);

   a_hom.FormLinearSystem(boundary_dofs, x_hom, b, A_hom, X_hom, B_hom);

   // 9. Solve the system using PCG with symmetric Gauss-Seidel preconditioner.
   GSSmoother M_eps(A_eps);
   PCG(A_eps, M_eps, B_eps, X_eps, 1, 200, 1e-12, 0.0);

   GSSmoother M_hom(A_hom);
   PCG(A_hom, M_hom, B_hom, X_hom, 1, 200, 1e-12, 0.0);

   // 10. Recover the solution x as a grid function and save to file. The output
   //     can be viewed using GLVis as follows: "glvis -m mesh.mesh -g sol.gf"
   a_eps.RecoverFEMSolution(X_eps, b, x_eps);

   a_hom.RecoverFEMSolution(X_hom, b, x_hom);

   x_eps.Save("sol_eps.gf");

   x_hom.Save("sol_hom.gf");

   mesh.Save("mesh.mesh");

   return 0;
}