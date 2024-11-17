# Julian Castrillon
# 2D Steady Incompressible Heat Transfer

import os
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

fe.parameters['form_compiler']['representation']     = 'uflacs'
fe.parameters['form_compiler']['optimize']           = True
fe.parameters['form_compiler']['cpp_optimize']       = True
fe.parameters["form_compiler"]["cpp_optimize_flags"] = '-O2 -funroll-loops'

Mesh = fe.Mesh('Mesh/2D Mesh/2DMesh.xml')

## Air Properties at 25ºC 1atm
#mu  = fe.Constant(0.00001825) # Dynamic viscosity              [kg/ms]
mu  = fe.Constant(0.01825)    # Dynamic viscosity               [kg/ms]
rho = fe.Constant(1.225)      # Density                         [kg/m3]
k   = fe.Constant(0.02551)    # Thermal conductivity            [W/mK]
cp  = fe.Constant(1007)       # Heat capacity constant pressure [J/kgK]
cv  = fe.Constant(718)        # Heat capacity constant volume   [J/kgK]
RR  = fe.Constant(287)        # Gas constant                    [J/kgK]
b   = fe.Constant((0, 0))     # Body accelerations              [m/s2]
           
Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
Temp = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres, Temp])
FS   = fe.FunctionSpace(Mesh, M)

TF        = fe.TrialFunction(FS)
(w, q, s) = fe.TestFunctions(FS)

TFsol     = fe.Function(FS)
(u, p, T) = fe.split(TFsol)

S = 0.5*(fe.grad(u)+fe.grad(u).T)      # Strain rate tensor
I = fe.Identity(Mesh.topology().dim()) # Identity tensor
Cst = -(p*I)+(2*mu*S)                  # Cauchy stress tensor

dx = fe.dx(metadata={"quadrature_degree":4})

Mnt = fe.inner(w,fe.grad(rho*u)*u)*dx \
    - fe.dot(w,(rho*b))*dx            \
    + fe.inner(fe.grad(w),Cst)*dx

Cnt = q*fe.div(rho*u)*dx

Egy = s*rho*cp*fe.dot(u,fe.grad(T))*dx    \
    + fe.dot(fe.grad(s), k*fe.grad(T))*dx 
    
WeakForm = Mnt + Cnt + Egy

J = fe.derivative(WeakForm, TFsol, TF)  

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DMesh_facet_region.xml')

Entry         = 8
BottomWall    = 9
Exit          = 10
TopWall       = 11
Circle        = 12
TriangleLeft  = 13
TriangleRight = 14

NoSlip    = fe.Constant((0, 0))
POut      = fe.Constant(0)
InletFlow = fe.Constant((0.2, 0))
Temp1     = fe.Constant(400)
Temp2     = fe.Constant(273)

EntryBC          = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
BottomWallBC     = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
ExitBC           = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
TopWallBC        = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)
CircleBC         = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, Circle)
TriangleLeftBC   = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleLeft)
TriangleRightBC  = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleRight)

EntryBCT         = fe.DirichletBC(FS.sub(2), Temp2,     DomainBoundaries, Entry)
BottomWallBCT    = fe.DirichletBC(FS.sub(2), Temp2,     DomainBoundaries, BottomWall)
TopWallBCT       = fe.DirichletBC(FS.sub(2), Temp2,     DomainBoundaries, TopWall)
CircleBCT        = fe.DirichletBC(FS.sub(2), Temp1,     DomainBoundaries, Circle)
TriangleLeftBCT  = fe.DirichletBC(FS.sub(2), Temp1,     DomainBoundaries, TriangleLeft)
TriangleRightBCT = fe.DirichletBC(FS.sub(2), Temp1,     DomainBoundaries, TriangleRight)

BCs = [EntryBC,  BottomWallBC,  ExitBC,  TopWallBC,  CircleBC,  TriangleLeftBC,  TriangleRightBC, \
       EntryBCT, BottomWallBCT, TopWallBCT, CircleBCT, TriangleLeftBCT, TriangleRightBCT]

InitialVel  = fe.interpolate(fe.Constant((fe.Constant(0.1), fe.Constant(0.1))), FS.sub(0).collapse()) 
InitialPres = fe.interpolate(fe.Constant(0.1), FS.sub(1).collapse())
InitialTemp = fe.interpolate(fe.Constant(273),  FS.sub(1).collapse())

fe.assign(TFsol.sub(0), InitialVel)
fe.assign(TFsol.sub(1), InitialPres)
fe.assign(TFsol.sub(2), InitialTemp)

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['linear_solver'] = 'petsc'
Parameters['newton_solver']['absolute_tolerance']   = 1e-7
Parameters['newton_solver']['relative_tolerance']   = 1e-6
Parameters['newton_solver']['maximum_iterations']   = 5
Parameters['newton_solver']['relaxation_parameter'] = 0.9
Parameters['newton_solver']['error_on_nonconvergence'] = False

Solver.solve()

FS2 = fe.FunctionSpace(Mesh, 'Lagrange', 1)
L   = fe.Constant(1.0)
vel = fe.sqrt(fe.dot(u,u))
Re  = (rho*vel*L)/mu
Reynolds = fe.project(Re, FS2)

(Velocity, Pressure, Temperature) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')
Temperature.rename('Temperature','Temperature')
Reynolds.rename('Reynolds','Reynolds')

FileVel  = fe.File('Results/Velocity.pvd')
FilePres = fe.File('Results/Pressure.pvd')
FileTemp = fe.File('Results/Temperature.pvd')
FileRey = fe.File('Results/Reynolds.pvd')
FileVel  << Velocity
FilePres << Pressure
FileTemp << Temperature
FileRey  << Reynolds
