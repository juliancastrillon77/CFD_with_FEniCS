# Julian Castrillon
# CFD - 2D Steady Compressible Navier Stokes Equations

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

vi = 1 # Inlet velocity [m/s]
## Air Properties at 25ºC 1atm
mu = fe.Constant(0.1849)     # Dynamic viscosity                [kg/ms]
k  = fe.Constant(0.02551)    # Thermal conductivity             [W/mK]
cp = fe.Constant(1007)       # Heat capacity constant pressure  [J/kgK]
cv = fe.Constant(718)        # Heat capacity constant volume    [J/kgK]
RR = fe.Constant(287)        # Gas constant                     [J/kgK]
B  = fe.Constant(1)          # Coefficient of thermal expansion [-]
b  = fe.Constant((0, -9.81)) # Body accelerations               [m/s2]
           
Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
Temp = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres, Temp])
FS   = fe.FunctionSpace(Mesh, M)

TF        = fe.TrialFunction(FS)
(w, q, s) = fe.TestFunctions(FS)

TFsol     = fe.Function(FS)
(v, p, T) = fe.split(TFsol)

WeakForm =   fe.dot(w*p/(RR*T),fe.grad(v)*v)*fe.dx     \
           - fe.div(w)*p*fe.dx                        \
           - fe.dot(w*p/(RR*T),b)*fe.dx                \
           + mu*fe.inner(fe.grad(w),fe.grad(v))*fe.dx \
               \
           + q*fe.div(v)*fe.dx \
               \
           - s*fe.dot(v,fe.grad(p))*fe.dx                            \
           + (cp/(k))*(p/(RR*T))*fe.dot(s,fe.dot(v,fe.grad(T)))*fe.dx \
           + fe.dot(fe.grad(s),fe.grad(T))*fe.dx

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
Temp1     = fe.Constant(293)
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
ExitBCT          = fe.DirichletBC(FS.sub(2), Temp2,     DomainBoundaries, Exit)
TopWallBCT       = fe.DirichletBC(FS.sub(2), Temp2,     DomainBoundaries, TopWall)
CircleBCT        = fe.DirichletBC(FS.sub(2), Temp1,     DomainBoundaries, Circle)
TriangleLeftBCT  = fe.DirichletBC(FS.sub(2), Temp1,     DomainBoundaries, TriangleLeft)
TriangleRightBCT = fe.DirichletBC(FS.sub(2), Temp1,     DomainBoundaries, TriangleRight)

BCs = [EntryBC,  BottomWallBC,  ExitBC,  TopWallBC,  CircleBC,  TriangleLeftBC,  TriangleRightBC, \
       EntryBCT, BottomWallBCT, ExitBCT, TopWallBCT, CircleBCT, TriangleLeftBCT, TriangleRightBCT]

InitialVel  = fe.interpolate(fe.Constant((fe.Constant(1.5), fe.Constant(0.2))), FS.sub(0).collapse()) 
InitialPres = fe.interpolate(fe.Constant(1e-3), FS.sub(1).collapse())
InitialTemp = fe.interpolate(fe.Constant(273),  FS.sub(1).collapse())

fe.assign(TFsol.sub(0), InitialVel)
fe.assign(TFsol.sub(1), InitialPres)
fe.assign(TFsol.sub(2), InitialTemp)

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['maximum_iterations']   = 5
Parameters['newton_solver']['relaxation_parameter'] = 1.0

Solver.solve()

(Velocity, Pressure, Temperature) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')
Temperature.rename('Temperature','Temperature')

FileVel  = fe.File('Results/Velocity.pvd')
FilePres = fe.File('Results/Pressure.pvd')
FileTemp = fe.File('Results/Temperature.pvd')
FileVel  << Velocity
FilePres << Pressure
FileTemp << Temperature
