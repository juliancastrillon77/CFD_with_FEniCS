# Julian Castrillon
# 2D Steady Stokes Flow

import os
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

fe.set_log_level(fe.LogLevel.PROGRESS)
fe.parameters['form_compiler']['representation']     = 'uflacs'
fe.parameters['form_compiler']['optimize']           = True
fe.parameters['form_compiler']['cpp_optimize']       = True
fe.parameters["form_compiler"]["cpp_optimize_flags"] = '-O2 -funroll-loops'

Mesh = fe.Mesh('Mesh/2D Mesh/2DMesh.xml')

# Fluid Properties
mu  = fe.Constant(1)      # Dynamic viscosity  [kg/ms]
rho = fe.Constant(1)      # Density            [kg/m3]
b   = fe.Constant((0, 0)) # Body accelerations [m/s2]

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres])
FS   = fe.FunctionSpace(Mesh, M)

(u, p) = fe.TrialFunctions(FS)
(w, q) = fe.TestFunctions(FS)

TF  = fe.Function(FS)

WeakForm =   mu*fe.inner(fe.grad(w),fe.grad(u))*fe.dx \
           - fe.div(w)*p*fe.dx                        \
           - rho*fe.dot(w,b)*fe.dx                    \
           - q*fe.div(u)*fe.dx

### Petrov Galerkin Pressure Stabilzation (PSPG) stabilization for pressure field // The Ladyzhenskaya-Babuska-Brezzi condition not met
#h    = fe.CellDiameter(Mesh)
#tau  = (1/3)*(h*h)/(4.0*mu)
#R    = fe.grad(p)
#PSPG = -tau*fe.inner(fe.grad(q),R)*fe.dx(metadata={'quadrature_degree':4})
#WeakForm += PSPG

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DMesh_facet_region.xml')

Entry         = 8
BottomWall    = 9
Exit          = 10
TopWall       = 11
Circle        = 12
TriangleLeft  = 13
TriangleRight = 14

NoSlip    = fe.Constant((0,0))
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

fe.assign(TF.sub(0), InitialVel)
fe.assign(TF.sub(1), InitialPres)

a = fe.lhs(WeakForm)
b = fe.rhs(WeakForm)

Problem = fe.LinearVariationalProblem(a, b, TF, BCs)
Solver  = fe.LinearVariationalSolver(Problem)

Solver.parameters['linear_solver']  = 'petsc'
Solver.parameters['preconditioner'] = 'ilu'
Solver.parameters['krylov_solver']['monitor_convergence'] = True
Solver.parameters['krylov_solver']['relative_tolerance'] = fe.Constant(1e-6)
Solver.parameters['krylov_solver']['absolute_tolerance'] = fe.Constant(1e-8)

Solver.solve()

(Velocity,Pressure) = TF.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')

FileVel  = fe.File('Results/Velocity.pvd')
FilePres = fe.File('Results/Pressure.pvd')
FileVel << Velocity
FilePres << Pressure
