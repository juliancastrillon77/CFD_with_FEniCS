# Julian Castrillon
# 3D Unsteady Heat Conduction

import os
import numpy as np
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

Mesh = fe.Mesh('Mesh/3D Mesh/3DMesh.xml')

Temp1 = fe.Constant(25)
Temp2 = fe.Constant(7)
alpha = fe.Constant(0.07) # Thermal diffusivity

dt = 0.01
t_end = 5
theta = fe.Constant(0.5)

FS = fe.FunctionSpace(Mesh, 'Lagrange', 1)

T = fe.TrialFunction(FS)
w = fe.TestFunction(FS)

TF = fe.Function(FS)

WeakFormB = alpha*fe.dot(fe.grad(w), fe.grad(T))*fe.dx

T0 = fe.Function(FS)

WeakFormF = alpha*fe.dot(fe.grad(w), fe.grad(T0))*fe.dx

WeakForm = fe.inner((T-T0)/dt,w)*fe.dx + (theta)*WeakFormB + (1-theta)*WeakFormF

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/3D Mesh/3DMesh_facet_region.xml')

LeftWall            = 46
RightWall           = 47
Entry               = 48
Exit                = 49
BottomWall          = 50
TopWall             = 51
SphereForwardRight  = 52
SphereForwardLeft   = 53
SphereBackwardLeft  = 54
SphereBackwardRight = 55
ConeForwardRight    = 56
ConeForwardLeft     = 57
ConeBackwardLeft    = 58
ConeBackwardRight   = 59

LeftWallBC            = fe.DirichletBC(FS, Temp1, DomainBoundaries, LeftWall)
RightWallBC           = fe.DirichletBC(FS, Temp1, DomainBoundaries, RightWall)
EntryBC               = fe.DirichletBC(FS, Temp1, DomainBoundaries, Entry)
ExitBC                = fe.DirichletBC(FS, Temp1, DomainBoundaries, Exit)
BottomWallBC          = fe.DirichletBC(FS, Temp1, DomainBoundaries, BottomWall)
TopWallBC             = fe.DirichletBC(FS, Temp1, DomainBoundaries, TopWall)
SphereForwardRightBC  = fe.DirichletBC(FS, Temp2, DomainBoundaries, SphereForwardRight)
SphereForwardLeftBC   = fe.DirichletBC(FS, Temp2, DomainBoundaries, SphereForwardLeft)
SphereBackwardLeftBC  = fe.DirichletBC(FS, Temp2, DomainBoundaries, SphereBackwardLeft)
SphereBackwardRightBC = fe.DirichletBC(FS, Temp2, DomainBoundaries, SphereBackwardRight)
ConeForwardRightBC    = fe.DirichletBC(FS, Temp2, DomainBoundaries, ConeForwardRight)
ConeForwardLeftBC     = fe.DirichletBC(FS, Temp2, DomainBoundaries, ConeForwardLeft)
ConeBackwardLeftBC    = fe.DirichletBC(FS, Temp2, DomainBoundaries, ConeBackwardLeft)
ConeBackwardRightBC   = fe.DirichletBC(FS, Temp2, DomainBoundaries, ConeBackwardRight)

BCs = [LeftWallBC, RightWallBC, EntryBC, ExitBC, BottomWallBC, TopWallBC, SphereForwardRightBC, \
       SphereForwardLeftBC, SphereBackwardLeftBC, SphereBackwardRightBC, ConeForwardRightBC,    \
       ConeForwardLeftBC, ConeBackwardLeftBC, ConeBackwardRightBC]

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

FileTemp = fe.File('Results/Temperature.pvd')

t  = dt
tn = 0

while t < t_end:

    print("t = ", np.around(t,3))
    print("Solving ....")
    
    Solver.solve()
    
    if tn%1 == 0:   # Save to file only if the time step increases by a step of 1 
        TF.rename('Temperature','Temperature')
        FileTemp  << TF
        print('Written temperature data')
    
    T0.assign(TF)
    t   += dt
    tn  += 1
    print('\n')

