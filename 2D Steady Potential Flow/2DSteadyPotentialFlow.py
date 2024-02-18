# Julian Castrillon
# CFD - 2d Steady Potential Flow / Irrotational Flow

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

Mesh = fe.Mesh('Mesh.xml')

h = fe.CellDiameter(Mesh)

vi  = 0.1
mu  = fe.Constant(1)
rho = fe.Constant(1)
b   = fe.Constant((0, 0))

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres])
FS   = fe.FunctionSpace(Mesh, M)

TF     = fe.TrialFunction(FS)
(w, q) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
(v, p) = fe.split(TFsol)

WeakForm = fe.inner(w,fe.grad(v)*v)*fe.dx \
           -fe.div(w)*p*fe.dx             \
           -fe.dot(w,b)*fe.dx             \
           -q*fe.div(v)*fe.dx

vnorm = fe.sqrt(fe.dot(v,v))
tau = ((2*vnorm/h)**2+9*(4/h**2)**2)**(-0.5)

#R = fe.grad(v)*v + fe.grad(p)
#
#R =   fe.grad(u)*u           \
#    - mu*fe.div(fe.grad(u))  \
#    + fe.grad(p)             \
#    - b) 

#WeakForm += tau*fe.inner(fe.grad(w)*v,R)*fe.dx(metadata={'quadrature_degree':8})

J = fe.derivative(WeakForm, TFsol, TF)

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh_facet_region.xml')

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

fe.assign(TFsol.sub(0), InitialVel)
fe.assign(TFsol.sub(1), InitialPres)

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['absolute_tolerance'] = 1e-4
Parameters['newton_solver']['relative_tolerance'] = 1e-4
Parameters['newton_solver']['maximum_iterations'] = 4
Parameters['newton_solver']['relaxation_parameter'] = 1.0

FileVel  = fe.File('Results/Velocity.pvd')
FilePres = fe.File('Results/Pressure.pvd')

Solver.solve()

(Velocity,Pressure) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')

FileVel << Velocity
FilePres << Pressure
