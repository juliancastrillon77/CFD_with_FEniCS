# Julian Castrillon
# CFD - 2D Steady Incompressible Navier Stokes Equation

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

Mesh = fe.Mesh(MeshFile)

U0  = 1
mu  = fe.Constant(1)
rho = fe.Constant(10)
b   = fe.Constant((0, -9.810))

V  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
P  = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M  = fe.MixedElement([V, P])
FS = fe.FunctionSpace(Mesh, M)

TF     = fe.TrialFunction(FS) # Note(1)
(w, q) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS) # Note(1)
(u, p) = fe.split(TFsol)

#################################################################################################
## Weak formulation #############################################################################
WeakForm =   rho*fe.inner(w,fe.grad(u)*u)*fe.dx       \
           - fe.div(w)*p*fe.dx                        \
           - rho*fe.dot(w,b)*fe.dx                    \
           + mu*fe.inner(fe.grad(w),fe.grad(u))*fe.dx \
           + q*fe.div(u)*fe.dx

J = fe.derivative(WeakForm, TFsol, TF)  

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

EntryBC         = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
BottomWallBC    = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
ExitBC          = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
TopWallBC       = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)
CircleBC        = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, Circle)
TriangleLeftBC  = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleLeft)
TriangleRightBC = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleRight)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, CircleBC, TriangleLeftBC, TriangleRightBC]

## Boundary check
#import sys
#fe.File('BoundaryCheck.pvd') << DomainBoundaries
#sys.exit()

#################################################################################################
## Solver #######################################################################################
Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['maximum_iterations']   = 7
Parameters['newton_solver']['relaxation_parameter'] = 1.0

VelFile = fe.File(OutFileV)
PreFile = fe.File(OutFileP)

Solver.solve()

(Velocity, Pressure) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')
VelFile << Velocity
PreFile << Pressure


# Note(1)
# - TrialFunction objects must always be used for the unknowns in the problem specification
#   while Function objects must be used for quantities that are computed.

#################################################################################################
## END ##########################################################################################
#################################################################################################