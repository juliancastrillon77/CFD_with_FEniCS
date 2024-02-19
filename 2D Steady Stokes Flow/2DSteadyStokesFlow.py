# Julian Castrillon
# CFD - 2D Steady Stokes Flow

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

h = fe.CellDiameter(Mesh)

vi  = 1
mu  = fe.Constant(1000)
rho = fe.Constant(1)
b   = fe.Constant((0, 0))

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres])
FS   = fe.FunctionSpace(Mesh, M)

(v, p) = fe.TrialFunctions(FS)
(w, q) = fe.TestFunctions(FS)

TF  = fe.Function(FS)

WeakForm = mu*fe.inner(fe.grad(w),fe.grad(v))*fe.dx \
           - p*fe.div(w)*fe.dx                      \
           - q*fe.div(v)*fe.dx

DomainBoundaries = fe.MeshFunction('size_t', Mesh, '../Mesh/2D Mesh/2DMesh_facet_region.xml')

Entry         = 8
BottomWall    = 9
Exit          = 10
TopWall       = 11
Circle        = 12
TriangleLeft  = 13
TriangleRight = 14

NoSlip    = fe.Constant((0,0))
POut      = fe.Constant(0)
InletFlow = fe.Constant((vi,0))

EntryBC         = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
BottomWallBC    = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
ExitBC          = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
TopWallBC       = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)
CircleBC        = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, Circle)
TriangleLeftBC  = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleLeft)
TriangleRightBC = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleRight)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, CircleBC, TriangleLeftBC, TriangleRightBC]

InitialVel  = fe.interpolate(fe.Constant((fe.Constant(vi), fe.Constant(0))), FS.sub(0).collapse()) 
InitialPres = fe.interpolate(fe.Constant(1e-3), FS.sub(1).collapse())

fe.assign(TF.sub(0), InitialVel)
fe.assign(TF.sub(1), InitialPres)

a = fe.lhs(WeakForm) # Bilinear form
b = fe.rhs(WeakForm) # Linear form

Problem = fe.LinearVariationalProblem(a, b, TF, BCs)
Solver  = fe.LinearVariationalSolver(Problem)

Solver.parameters['linear_solver']  = 'cg'
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
