# Julian Castrillon
# CFD - 2D Steady Poisson Equation

import os
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

Mesh = fe.Mesh('../Mesh/2D Mesh/2DMesh.xml')

FS = fe.FunctionSpace(Mesh, 'Lagrange', 1)

u = fe.TrialFunction(FS)
w = fe.TestFunction(FS)

TF = fe.Function(FS)

DomainBoundaries = fe.MeshFunction('size_t', Mesh, '../Mesh/2D Mesh/2DMesh_facet_region.xml')

Entry         = 8
BottomWall    = 9
Exit          = 10
TopWall       = 11
Circle        = 12
TriangleLeft  = 13
TriangleRight = 14

Expres = fe.Expression('x[1]*(1-x[1])', degree = 2)
Const  = fe.Constant(0)

EntryBC         = fe.DirichletBC(FS, Const,  DomainBoundaries, Entry)
BottomWallBC    = fe.DirichletBC(FS, Const,  DomainBoundaries, BottomWall)
ExitBC          = fe.DirichletBC(FS, Const,  DomainBoundaries, Exit)
TopWallBC       = fe.DirichletBC(FS, Const,  DomainBoundaries, TopWall)
CircleBC        = fe.DirichletBC(FS, Expres, DomainBoundaries, Circle)
TriangleLeftBC  = fe.DirichletBC(FS, Expres, DomainBoundaries, TriangleLeft)
TriangleRightBC = fe.DirichletBC(FS, Expres, DomainBoundaries, TriangleRight)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, CircleBC, TriangleLeftBC, TriangleRightBC]

a = fe.dot(fe.grad(w), fe.grad(u))*fe.dx
f = fe.Expression('x[0]*x[1]', degree = 2)
b = w*f*fe.dx

Problem = fe.LinearVariationalProblem(a, -b, TF, BCs)
Solver  = fe.LinearVariationalSolver(Problem)

Solver.parameters['linear_solver']  = 'cg'
Solver.parameters['preconditioner'] = 'ilu'
Solver.parameters['krylov_solver']['monitor_convergence'] = True
Solver.parameters['krylov_solver']['relative_tolerance'] = fe.Constant(1e-6)
Solver.parameters['krylov_solver']['absolute_tolerance'] = fe.Constant(1e-8)

Solver.solve()

TF.rename('u','u')
FileU = fe.File('Results/u.pvd')
FileU << TF











