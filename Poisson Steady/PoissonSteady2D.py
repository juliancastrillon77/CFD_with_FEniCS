# Julian Castrillon P.
# Poisson Steady Equation 2D

import os
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

##

mesh = fe.Mesh('Domain2D.xml')
V = fe.FunctionSpace(mesh, 'Lagrange', 1)

## Boundaries
DomainBoundaries = fe.MeshFunction('size_t', mesh, 'Domain2D_facet_region.xml') # size_t ?

Entry = 35
Exit1 = 37
Exit2 = 36
Walls = 38

Expres = fe.Expression('x[1]*(1-x[1])', degree = 2)
Const = fe.Constant(0)

EntryBC = fe.DirichletBC(V, Const, DomainBoundaries,Entry)
Exit1BC = fe.DirichletBC(V, Const, DomainBoundaries,Exit1)
Exit2BC = fe.DirichletBC(V, Const, DomainBoundaries,Exit2)
WallsBC = fe.DirichletBC(V, Expres, DomainBoundaries,Walls)
BCset = [EntryBC,Exit1BC,Exit2BC,WallsBC]

## Define Weak Formulation
u = fe.TrialFunction(V)
w = fe.TestFunction(V)

a = fe.dot(fe.grad(w), fe.grad(u))*fe.dx
f = fe.Expression('x[0]*x[1]', degree = 2)
b = w*f*fe.dx

Solution = fe.Function(V)
#fe.solve(a==-b, Solution, BCset)

Problem = fe.LinearVariationalProblem(a, -b, Solution, BCset)
Solver  = fe.LinearVariationalSolver(Problem)
Solver.parameters['linear_solver'] = 'cg'
Solver.parameters['preconditioner']  = 'sor'
Solver.parameters['krylov_solver']['monitor_convergence'] = True
Solver.parameters['krylov_solver']['relative_tolerance'] = fe.Constant(1e-6)
Solver.parameters['krylov_solver']['absolute_tolerance'] = fe.Constant(1e-8)
Solver.solve()

vtkFile = fe.File('Results/PoissonSteady2D.pvd')
Solution.rename('u','')

vtkFile << Solution









