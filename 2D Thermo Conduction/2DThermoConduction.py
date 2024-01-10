# Julian Castrillon
# CFD - 2D Thermo Conduction

import os
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

#################################################################################################
## Parameters ###################################################################################
MeshFile  = 'Mesh.xml'
FacetFile = 'Mesh_facet_region.xml'
OutFileT  = 'Results/Temperature.pvd'

Mesh = fe.Mesh(MeshFile)

Te = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
FS = fe.FunctionSpace(Mesh, Te)

T = fe.TrialFunction(FS)
w = fe.TestFunction(FS)

TF = fe.Function(FS)
#################################################################################################
## Weak formulation #############################################################################

WeakForm = fe.dot(fe.grad(w), fe.grad(T))*fe.dx

#################################################################################################
## Boundary Conditions ##########################################################################
DomainBoundaries = fe.MeshFunction('size_t', Mesh, FacetFile)
ds = fe.ds(subdomain_data = DomainBoundaries)

Entry         = 8
BottomWall    = 9
Exit          = 10
TopWall       = 11
Circle        = 12
TriangleLeft  = 13
TriangleRight = 14

Temp1 = fe.Constant(25)
Temp2 = fe.Constant(7)

EntryBC         = fe.DirichletBC(FS, Temp1, DomainBoundaries, Entry)
BottomWallBC    = fe.DirichletBC(FS, Temp1, DomainBoundaries, BottomWall)
ExitBC          = fe.DirichletBC(FS, Temp1, DomainBoundaries, Exit)
TopWallBC       = fe.DirichletBC(FS, Temp1, DomainBoundaries, TopWall)
CircleBC        = fe.DirichletBC(FS, Temp2, DomainBoundaries, Circle)
TriangleLeftBC  = fe.DirichletBC(FS, Temp2, DomainBoundaries, TriangleLeft)
TriangleRightBC = fe.DirichletBC(FS, Temp1, DomainBoundaries, TriangleRight)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, CircleBC, TriangleLeftBC, TriangleRightBC]

## Boundary check
#import sys
#fe.File('BoundaryCheck.pvd') << DomainBoundaries
#sys.exit()

#################################################################################################
## Solver #######################################################################################
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

TempFile = fe.File(OutFileT)
TF.rename('Temperature','Temperature')
TempFile << TF

#################################################################################################
## END ##########################################################################################
#################################################################################################