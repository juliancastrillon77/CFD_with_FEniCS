# Julian Castrillon
# CFD - 2D Unsteady Stokes Flow

import os
import numpy as np
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

fe.set_log_level(fe.LogLevel.INFO)
fe.parameters['form_compiler']['representation']     = 'uflacs'
fe.parameters['form_compiler']['optimize']           = True
fe.parameters['form_compiler']['cpp_optimize']       = True
fe.parameters["form_compiler"]["cpp_optimize_flags"] = '-O2 -funroll-loops'

Mesh = fe.Mesh('../Mesh/2D Mesh/2DMesh.xml')

dt    = 0.01
t_end = 0.03

# Fluid Properties
mu  = fe.Constant(1)      # Dynamic viscosity  [kg/ms]
rho = fe.Constant(1)      # Density            [kg/m3]
b   = fe.Constant((0, 0)) # Body accelerations [m/s2]

theta = fe.Constant(0.5)      
n     = fe.FacetNormal(Mesh)   
h     = fe.CellDiameter(Mesh) 

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres])
FS   = fe.FunctionSpace(Mesh, M)

TF     = fe.TrialFunction(FS)
(w, q) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
(v, p) = fe.split(TFsol)

WeakFormB =   mu*fe.inner(fe.grad(w),fe.grad(v))*fe.dx \
           - fe.div(w)*p*fe.dx                        \
           - rho*fe.dot(w,b)*fe.dx
 
WeakFormBB = - q*fe.div(v)*fe.dx

TFsol0   = fe.Function(FS)
(v0, p0) = fe.split(TFsol0)

WeakFormF =   mu*fe.inner(fe.grad(w),fe.grad(v0))*fe.dx \
           - fe.div(w)*p*fe.dx                        \
           - rho*fe.dot(w,b)*fe.dx

WeakFormFF = - q*fe.div(v0)*fe.dx
                

F1 = fe.inner((v-v0)/dt,w)*fe.dx + (theta)*WeakFormB + (1-theta)*WeakFormF
F2 = fe.inner((v-v0)/dt,q)*fe.dx + (theta)*WeakFormBB + (1-theta)*WeakFormFF
WeakForm = F1+F2

#WeakForm = fe.inner((v-v0)/dt,w)*fe.dx + (theta)*WeakFormB + (1-theta)*WeakFormF


### Petrov Galerkin Pressure Stabilzation (PSPG) stabilization for pressure field // The Ladyzhenskaya-Babuska-Brezzi condition not met
#h    = fe.CellDiameter(Mesh)
#tau  = (1/3)*(h*h)/(4.0*mu)
#R    = fe.grad(p)
#PSPG = -tau*fe.inner(fe.grad(q),R)*fe.dx(metadata={'quadrature_degree':4})
#WeakForm += PSPG

DomainBoundaries = fe.MeshFunction('size_t', Mesh, '../Mesh/2D Mesh/2DMesh_facet_region.xml')

Entry         = 8
BottomWall    = 9
Exit          = 10
TopWall       = 11
Circle        = 12
TriangleLeft  = 13
TriangleRight = 14

NoSlip    = fe.Constant((0, 0))
POut      = fe.Constant(0)
InletFlow = fe.Expression(['6*x[1]*(1-x[1])','0'], degree=2)

EntryBC         = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
BottomWallBC    = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
ExitBC          = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
TopWallBC       = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)
CircleBC        = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, Circle)
TriangleLeftBC  = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleLeft)
TriangleRightBC = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleRight)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, CircleBC, TriangleLeftBC, TriangleRightBC]

InitialVel  = fe.interpolate(fe.Constant((fe.Constant(1), fe.Constant(0))), FS.sub(0).collapse()) 
InitialPres = fe.interpolate(fe.Constant(10), FS.sub(1).collapse())

fe.assign(TFsol.sub(0),  InitialVel)
fe.assign(TFsol.sub(1),  InitialPres)

a = fe.lhs(WeakForm)
b = fe.rhs(WeakForm)

Problem = fe.LinearVariationalProblem(a, b, TFsol, BCs)
Solver  = fe.LinearVariationalSolver(Problem)

Solver.parameters['linear_solver']  = 'petsc'
Solver.parameters['preconditioner'] = 'ilu'
Solver.parameters['krylov_solver']['monitor_convergence'] = True
Solver.parameters['krylov_solver']['relative_tolerance'] = fe.Constant(1e-6)
Solver.parameters['krylov_solver']['absolute_tolerance'] = fe.Constant(1e-8)

FileVel  = fe.File('Results/Velocity.pvd')
FilePres = fe.File('Results/Pressure.pvd')

t  = dt
tn = 0

while t < t_end:

    print("t = ", np.around(t,3))
    print("Solving ....")
    
    Solver.solve()
    
    (Velocity, Pressure) = TFsol.split(deepcopy = True)
    
    if tn%1 == 0:   # Save to file only if the time step increases by a step of 1 
        Velocity.rename('Velocity', 'Velocity')
        Pressure.rename('Pressure', 'Pressure')
        FileVel  << Velocity
        FilePres << Pressure
        print('Written velocity & pressure data')
    
    TFsol0.assign(TFsol)
    t   += dt
    tn  += 1
    print('\n')