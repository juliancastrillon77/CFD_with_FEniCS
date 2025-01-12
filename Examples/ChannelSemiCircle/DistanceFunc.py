# Julian Castrillon
# Minimum distance from wall function
import os
import fenics as fe
import numpy as np
import matplotlib.pyplot as plt

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

Mesh = fe.Mesh('Mesh/2D Mesh/2DChannelSemiCircle.xml')

FS = fe.FunctionSpace(Mesh, 'Lagrange', 1)

d = fe.Function(FS)
w = fe.TestFunction(FS)

Weakform = fe.dot(fe.grad(w),fe.grad(d))*fe.dx \
         - fe.Constant(1)*w*fe.dx

J = fe.derivative(Weakform, d)  

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DChannelSemiCircle_facet_region.xml')

Entry      = 1
BottomWall = 2
Exit       = 3
TopWall    = 4

Const  = fe.Constant(0)

EntryBC      = fe.DirichletBC(FS, Const,  DomainBoundaries, Entry)
BottomWallBC = fe.DirichletBC(FS, Const,  DomainBoundaries, BottomWall)
ExitBC       = fe.DirichletBC(FS, Const,  DomainBoundaries, Exit)
TopWallBC    = fe.DirichletBC(FS, Const,  DomainBoundaries, TopWall)

BCs = [BottomWallBC, TopWallBC]

InitialVel  = fe.interpolate(fe.Constant(1), FS) 
fe.assign(d,  InitialVel)

solutions = []

for i, BC in enumerate(BCs):
   
    Problem = fe.NonlinearVariationalProblem(Weakform, d, BC, J)
    Solver  = fe.NonlinearVariationalSolver(Problem)

    Parameters = Solver.parameters
    Parameters['newton_solver']['absolute_tolerance']   = 1e-8
    Parameters['newton_solver']['relative_tolerance']   = 1e-7
    Parameters['newton_solver']['maximum_iterations']   = 20
    Parameters['newton_solver']['relaxation_parameter'] = 1.0

    Solver.solve()
    solutions.append(d.copy(deepcopy=True))

MinDist = fe.Function(FS)
FinalVec = solutions[0].vector()[:]

for sol in solutions[1:]:
    FinalVec = np.minimum(FinalVec, sol.vector()[:])

MinDist.vector().set_local(FinalVec)

# d = fe.sqrt(fe.assemble(fe.dot(fe.grad(MinDist),fe.grad(MinDist))*fe.dx)+2*MinDist) \
#   - fe.sqrt(fe.assemble(fe.dot(MinDist,MinDist)*fe.dx))
# E = fe.project(d, FS)

MinDist.rename('Distance','Distance')
fe.File(f'Results/Distance.pvd') << MinDist
fe.File(f'Results/Distance.xml') << MinDist

