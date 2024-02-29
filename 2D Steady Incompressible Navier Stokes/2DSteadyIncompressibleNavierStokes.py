# Julian Castrillon
# CFD - 2D Steady Incompressible Navier Stokes Equations

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

Mesh = fe.Mesh('../Mesh/2D Mesh/2DMesh.xml')

vi  = 1                       # Inlet velocity [m/s]
mu  = fe.Constant(0.001002)   # Dynamic viscosity  [kg/ms]
rho = fe.Constant(1000)       # Density            [kg/m3]
b   = fe.Constant((0, -9.81)) # Body accelerations [m/s2]

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres])
FS   = fe.FunctionSpace(Mesh, M)

TF     = fe.TrialFunction(FS)
(w, q) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
(v, p) = fe.split(TFsol)

WeakForm =   rho*fe.inner(w,fe.grad(v)*v)*fe.dx       \
           - fe.div(w)*p*fe.dx                        \
           - rho*fe.dot(w,b)*fe.dx                    \
           + mu*fe.inner(fe.grad(w),fe.grad(v))*fe.dx \
           + q*fe.div(v)*fe.dx

J = fe.derivative(WeakForm, TFsol, TF)  

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
InletFlow = fe.Constant((vi, 0))

EntryBC         = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
BottomWallBC    = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
ExitBC          = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
TopWallBC       = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)
CircleBC        = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, Circle)
TriangleLeftBC  = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleLeft)
TriangleRightBC = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleRight)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, CircleBC, TriangleLeftBC, TriangleRightBC]

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['maximum_iterations']   = 20
Parameters['newton_solver']['relaxation_parameter'] = 1.0

Solver.solve()

(Velocity, Pressure) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')

FileVel  = fe.File('Results/Velocity.pvd')
FilePres = fe.File('Results/Pressure.pvd')
FileVel  << Velocity
FilePres << Pressure