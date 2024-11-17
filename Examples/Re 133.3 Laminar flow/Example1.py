# Julian Castrillon
# 2D Steady Incompressible

import os
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

fe.parameters['form_compiler']['representation']     = 'uflacs'
fe.parameters['form_compiler']['optimize']           = True
fe.parameters['form_compiler']['cpp_optimize']       = True
fe.parameters["form_compiler"]["cpp_optimize_flags"] = '-O2 -funroll-loops'

Mesh = fe.Mesh('Mesh/2D Mesh/2DChannelSemiCircle.xml')

mu  = fe.Constant(0.001)  # Dynamic viscosity  [kg/ms]
rho = fe.Constant(1)      # Density            [kg/m3]
b   = fe.Constant((0, 0)) # Body accelerations [m/s2]

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres])
FS   = fe.FunctionSpace(Mesh, M)

TF     = fe.TrialFunction(FS)
(w, q) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
(u, p) = fe.split(TFsol)

S = 0.5*(fe.grad(u)+fe.grad(u).T)      # Strain rate tensor
I = fe.Identity(Mesh.topology().dim()) # Identity tensor
Cst = -(p*I)+(2*mu*S)                  # Cauchy stress tensor

dx = fe.dx(metadata={"quadrature_degree":4})

Mnt = fe.inner(w,fe.grad(rho*u)*u)*dx \
    - fe.dot(w,(rho*b))*dx            \
    + fe.inner(fe.grad(w),Cst)*dx

Cnt = q*fe.div(rho*u)*dx

WeakForm = Mnt + Cnt

J = fe.derivative(WeakForm, TFsol, TF)  

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DChannelSemiCircle_facet_region.xml')

Entry         = 1
BottomWall    = 2
Exit          = 3
TopWall       = 4

NoSlip    = fe.Constant((0, 0))
POut      = fe.Constant(0)
InletFlow = fe.Expression(['2*((4*x[1]*(0.1-x[1]))/(0.1*0.1))','0'], degree=2)

EntryBC         = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
BottomWallBC    = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
ExitBC          = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
TopWallBC       = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC]

InitialVel  = fe.interpolate(fe.Constant((fe.Constant(1), fe.Constant(0.01))), FS.sub(0).collapse()) 
InitialPres = fe.interpolate(fe.Constant(10), FS.sub(1).collapse())

fe.assign(TFsol.sub(0),  InitialVel)
fe.assign(TFsol.sub(1),  InitialPres)

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['maximum_iterations']   = 20
Parameters['newton_solver']['relaxation_parameter'] = 1.0

Solver.solve()

FS2 = fe.FunctionSpace(Mesh, 'Lagrange', 1)
L   = fe.Constant(0.1)
uMag = fe.sqrt(fe.dot(u,u))
Re  = (rho*uMag*L)/mu
Reynolds = fe.project(Re, FS2)

SSFS = fe.TensorFunctionSpace(Mesh, 'P', 1)
ShearStressV = (2*mu*S) 
ShearStress  = fe.project(ShearStressV, SSFS)

(Velocity, Pressure) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')
Reynolds.rename('Reynolds','Reynolds')
ShearStress.rename('ShearStress','ShearStress')
DomainBoundaries.rename('DomainBoundaries','DomainBoundaries')

fe.File('Results/Velocity.xml')    << Velocity
fe.File('Results/Pressure.xml')    << Pressure
fe.File('Results/Reynolds.pvd')    << Reynolds
fe.File('Results/ShearStress.xml') << ShearStress
fe.File('Results/Boundaries.pvd')  << DomainBoundaries






