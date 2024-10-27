# Julian Castrillon
# 2D Heat Conduction

import os
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

Mesh = fe.Mesh('Mesh/2D Mesh/2DMesh.xml')

Temp1 = fe.Constant(25)
Temp2 = fe.Constant(7)

FS = fe.FunctionSpace(Mesh, 'Lagrange', 1)

T = fe.TrialFunction(FS)
w = fe.TestFunction(FS)

TF = fe.Function(FS)

WeakForm = fe.dot(fe.grad(w), fe.grad(T))*fe.dx

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DMesh_facet_region.xml')

Entry         = 8
BottomWall    = 9
Exit          = 10
TopWall       = 11
Circle        = 12
TriangleLeft  = 13
TriangleRight = 14

EntryBC         = fe.DirichletBC(FS, Temp1, DomainBoundaries, Entry)
BottomWallBC    = fe.DirichletBC(FS, Temp1, DomainBoundaries, BottomWall)
ExitBC          = fe.DirichletBC(FS, Temp1, DomainBoundaries, Exit)
TopWallBC       = fe.DirichletBC(FS, Temp1, DomainBoundaries, TopWall)
CircleBC        = fe.DirichletBC(FS, Temp2, DomainBoundaries, Circle)
TriangleLeftBC  = fe.DirichletBC(FS, Temp2, DomainBoundaries, TriangleLeft)
TriangleRightBC = fe.DirichletBC(FS, Temp2, DomainBoundaries, TriangleRight)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, CircleBC, TriangleLeftBC, TriangleRightBC]

a = fe.lhs(WeakForm)
b = fe.rhs(WeakForm)

Problem = fe.LinearVariationalProblem(a, b, TF, BCs)
Solver  = fe.LinearVariationalSolver(Problem)

Parameters = Solver.parameters
Solver.parameters['linear_solver']  = 'cg'
Solver.parameters['preconditioner'] = 'ilu'
Solver.parameters['krylov_solver']['monitor_convergence'] = True
Solver.parameters['krylov_solver']['relative_tolerance'] = fe.Constant(1e-6)
Solver.parameters['krylov_solver']['absolute_tolerance'] = fe.Constant(1e-8)

Solver.solve()

TF.rename('Temperature','Temperature')

FileTemp = fe.File('Results/Temperature.pvd')
FileTemp << TF

