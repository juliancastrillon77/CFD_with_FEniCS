# Julian Castrillon
# CFD - 2D Steady Compressible Navier Stokes Equation

import os
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

#################################################################################################
## Definition of some global FEniCS optimization parameters #####################################
fe.set_log_level(fe.LogLevel.PROGRESS)
fe.parameters['form_compiler']['representation']     = 'uflacs'
fe.parameters['form_compiler']['optimize']           = True
fe.parameters['form_compiler']['cpp_optimize']       = True
fe.parameters["form_compiler"]["cpp_optimize_flags"] = '-O2 -funroll-loops'

#################################################################################################
## Parameters ###################################################################################
MeshFile  = 'Mesh.xml'
FacetFile = 'Mesh_facet_region.xml'
OutFileV  = 'Results/Velocity.pvd'
OutFileP  = 'Results/Pressure.pvd'
OutFileT  = 'Results/Temperature.pvd'

Mesh = fe.Mesh(MeshFile)

## Air Properties at 25ºC 1atm
mu = fe.Constant(0.1849)     # Dynamic viscosity [kg/ms]
k  = fe.Constant(0.02551)    # Thermal conductivity [W/mK]
cp = fe.Constant(1007)       # Heat capacity constant pressure [J/kgK]
cv = fe.Constant(718)        # Heat capacity constant volume[J/kgK]
R  = fe.Constant(287)        # Gas constant [J/kgK]
B  = fe.Constant(1)          # Coefficient of thermal expansion
b  = fe.Constant((0, -9.81)) # Body forces []

U0  = 1

Velo  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres  = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
Temp  = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M     = fe.MixedElement([Velo, Pres, Temp])
FS    = fe.FunctionSpace(Mesh, M)

TF        = fe.TrialFunction(FS)
(w, q, s) = fe.TestFunctions(FS)

TFsol     = fe.Function(FS)
(u, P, T) = fe.split(TFsol)

#################################################################################################
## Weak formulation #############################################################################
WeakForm =   fe.dot(w*P/(R*T),fe.grad(u)*u)*fe.dx \
           - fe.div(w)*P*fe.dx                                        \
           - fe.dot(w*P/(R*T),b)*fe.dx              \
           + mu*fe.inner(fe.grad(w),fe.grad(u))*fe.dx                 \
               \
               \
           + q*fe.div(u)*fe.dx    \
               \
               \
           - s*fe.dot(u,fe.grad(P))*fe.dx \
           + (cp/(k))*(P/(R*T))*fe.dot(s,fe.dot(u,fe.grad(T)))*fe.dx \
           + fe.dot(fe.grad(s),fe.grad(T))*fe.dx

#################################################################################################
## Boundary Conditions ##########################################################################
DomainBoundaries = fe.MeshFunction('size_t', Mesh, FacetFile)
ds = fe.ds(subdomain_data = DomainBoundaries)

Entry         = 8
BottomWall    = 9
Exit          = 10
TopWall       = 11
Circle        = 12
TriangleLeft  = 13
TriangleRight = 14

NoSlip     = fe.Constant((0, 0))
POut       = fe.Constant(0)
InletFlow  = fe.Constant((U0, 0))
Temp1      = fe.Constant(293)
Temp2      = fe.Constant(273)

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


BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, CircleBC, TriangleLeftBC, TriangleRightBC, \
       EntryBCT, BottomWallBCT, ExitBCT, TopWallBCT, CircleBCT, TriangleLeftBCT, TriangleRightBCT]

## Boundary check
#import sys
#fe.File('BoundaryCheck.pvd') << DomainBoundaries
#sys.exit()

## Initial conditions
Intu = fe.interpolate(fe.Constant((fe.Constant(1.5), fe.Constant(0.2))), FS.sub(0).collapse()) 
IntP = fe.interpolate(fe.Constant(1e-3), FS.sub(1).collapse())
IntT = fe.interpolate(fe.Constant(273),  FS.sub(1).collapse())

fe.assign(TFsol.sub(0), Intu)
fe.assign(TFsol.sub(1), IntP)
fe.assign(TFsol.sub(2), IntT)

## Jacobian
J = fe.derivative(WeakForm, TFsol, TF)  

#################################################################################################
## Solver #######################################################################################
Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['maximum_iterations']   = 5
Parameters['newton_solver']['relaxation_parameter'] = 1.0

VelFile = fe.File(OutFileV)
PreFile = fe.File(OutFileP)
TemFile = fe.File(OutFileT)

Solver.solve()

(Velocity, Pressure, Temperature) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')
Temperature.rename('Temperature','Temperature')
VelFile << Velocity
PreFile << Pressure
TemFile << Temperature

#################################################################################################
## END ##########################################################################################
#################################################################################################