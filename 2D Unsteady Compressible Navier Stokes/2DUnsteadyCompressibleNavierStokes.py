# Julian Castrillon
# CFD - 2D Compressible Unsteady Navier Stokes Convection

import os
import numpy as np
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

#################################################################################################
## Definition of some global FEniCS optimization parameters #####################################
fe.set_log_level(fe.LogLevel.INFO)
fe.parameters['form_compiler']['representation']     = 'uflacs'
fe.parameters['form_compiler']['optimize']           = True
fe.parameters['form_compiler']['cpp_optimize']       = True
fe.parameters["form_compiler"]["cpp_optimize_flags"] = '-O2 -funroll-loops'

#################################################################################################
## Parameters ###################################################################################
MeshFile  = 'Mesh.xml'                # Input mesh
FacetFile = 'Mesh_facet_region.xml'   # Input mesh boundaries
OutFileV  = 'Results/Velocity.pvd'    # Output file velocity
OutFileP  = 'Results/Pressure.pvd'    # Output file pressure
OutFileT  = 'Results/Temperature.pvd' # Output file temperature

Mesh = fe.Mesh(MeshFile)     # Create the mesh and import mesh into the solver

# Air Properties at 25ºC 1atm
mu = fe.Constant(1/500)     # Dynamic viscosity [kg/ms]
k  = fe.Constant(0.1)    # Thermal conductivity [W/mK]
cp = fe.Constant(100)       # Heat capacity constant pressure [J/kgK]
cv = fe.Constant(718)        # Heat capacity constant volume[J/kgK]
RRR  = fe.Constant(287)        # Gas constant [J/kgK]
B  = fe.Constant(1)          # Coefficient of thermal expansion
b  = fe.Constant((0, 0)) # Body forces []
rho = fe.Constant(1)

dt    = 0.01 # Timestep
t_end = 10  # Length of simulation

U0    = 0.5    # Inlet Velocity
theta = fe.Constant(0.5)      # Crank-Nicholson time-stepping scheme
n     = fe.FacetNormal(Mesh)  # Normal vector  
h     = fe.CellDiameter(Mesh) # Mesh size

# Define the mixed vector function space operating on this meshed domain
Velo  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres  = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
Temp  = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M     = fe.MixedElement([Velo, Pres, Temp])
FS    = fe.FunctionSpace(Mesh, M) # Function space

TF        = fe.TrialFunction(FS) # Trial function
(w, q, s) = fe.TestFunctions(FS) # Test  functions

#################################################################################################
## Weak formulation #############################################################################
TFsol     = fe.Function(FS) # Timestep n+1
(u, P, T) = fe.split(TFsol)

WeakFormB =   fe.dot(w*rho,fe.grad(u)*u)*fe.dx         \
            - fe.div(w)*P*fe.dx                        \
            - fe.dot(w*rho,b)*fe.dx                    \
            + mu*fe.inner(fe.grad(w),fe.grad(u))*fe.dx \
            \
            \
            + q*fe.div(u)*fe.dx \
            \
            \
            - s*fe.dot(u,fe.grad(P))*fe.dx                      \
            + (cp/(k))*rho*fe.dot(s,fe.dot(u,fe.grad(T)))*fe.dx \
            + fe.dot(fe.grad(s),fe.grad(T))*fe.dx               \

TFsol0       = fe.Function(FS) # Timestep n
(u0, P0, T0) = fe.split(TFsol0)

WeakFormF =   fe.dot(w*rho,fe.grad(u0)*u0)*fe.dx        \
            - fe.div(w)*P*fe.dx                         \
            - fe.dot(w*rho,b)*fe.dx                     \
            + mu*fe.inner(fe.grad(w),fe.grad(u0))*fe.dx \
            \
            \
            + q*fe.div(u0)*fe.dx \
            \
            \
            - s*fe.dot(u0,fe.grad(P))*fe.dx                        \
            + (cp/(k))*(rho)*fe.dot(s,fe.dot(u0,fe.grad(T)))*fe.dx \
            + fe.dot(fe.grad(s),fe.grad(T))*fe.dx
                
WeakForm = fe.inner((u-u0)/dt,w)*fe.dx + (theta)*WeakFormB + (1-theta)*WeakFormF

### Streamwise Upwind Petrov Galerkin (SUPG) stabilization for convection #######################
vnorm = fe.sqrt(fe.dot(u0,u0))
tau   = ((1/(theta*dt))**2+(2*vnorm/h)**2+9*(4*mu/h**2)**2)**(-0.5)

# Residual of the strong form of Navier-Stokes and continuity
#R =   idt*(u-u0)             \
#    + theta *                \
#    ( fe.grad(u)*u           \
#    - mu*fe.div(fe.grad(u))  \
#    + fe.grad(P)             \
#    - b)                     \
#            

R = (u-u0)/dt                \
    + theta *                \
    (fe.grad(u)*u            \
    + fe.grad(P)             \
    - rho*b                  \
    - mu*fe.div(fe.grad(u))) \
    \
    + (1-theta) *            \
    (fe.grad(u0)*u0            \
    + fe.grad(P)             \
    - rho*b                  \
    - mu*fe.div(fe.grad(u0)))# \
#\
#\
#+ q*fe.div(u)*fe.dx \
#\
#\
#- s*fe.dot(u,fe.grad(P))*fe.dx                            \
#+ (cp/(k))*rho*fe.dot(s,fe.dot(u,fe.grad(T)))*fe.dx \
#+ fe.dot(fe.grad(s),fe.grad(T))*fe.dx                     \

#    + (1-theta) *            \
#    ( fe.grad(u0)*u0         \
#    - mu*fe.div(fe.grad(u0)) \
#    + fe.grad(P)             \
#    - b)

WeakForm += tau*fe.inner(fe.grad(w)*u0,R)*fe.dx(metadata={'quadrature_degree':8})

# Petrov Galerkin Pressure Stabilzation (PSPG) stabilization for pressure field // The Ladyzhenskaya-Babuska-Brezzi condition not met
WeakForm += -tau*fe.inner(fe.grad(s),R)*fe.dx(metadata={'quadrature_degree':8})

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
Temp1      = fe.Constant(340)
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
# fe.File('New.pvd') << domainBoundaries
# sys.exit()

## Initial conditions
Intu = fe.interpolate(fe.Constant((fe.Constant(1.5), fe.Constant(0.2))), FS.sub(0).collapse()) 
IntP = fe.interpolate(fe.Constant(1e-3), FS.sub(1).collapse())
IntT = fe.interpolate(fe.Constant(273),  FS.sub(1).collapse())

fe.assign(TFsol.sub(0), Intu)
fe.assign(TFsol.sub(1), IntP)
fe.assign(TFsol.sub(2), IntT)
fe.assign(TFsol0.sub(0), Intu)
fe.assign(TFsol0.sub(1), IntP)
fe.assign(TFsol0.sub(2), IntT)

# Jacobian
J = fe.derivative(WeakForm, TFsol, TF) 

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

# Create files for storing solution
VelFile = fe.File(OutFileV)
PreFile = fe.File(OutFileP)
TemFile = fe.File(OutFileT)

# Loops #####################################
# Inner loop - Write data so many time steps
# Outer loop - Time steps

t  = dt
tn = 0
   
while t < t_end:

    print("t = ", np.around(t,3))
    print("Solving ....")
    
    Solver.solve()
      
    (Velocity, Pressure, Temperature) = TFsol.split(deepcopy = True)
    
    if tn%1 == 0:   # Save to file only if the time step increases by a step of 1 
               
        Velocity.rename('Velocity','Velocity')
        Pressure.rename('Pressure','Pressure')
        Temperature.rename('Temperature','Temperature')
        VelFile << Velocity
        PreFile << Pressure
        TemFile << Temperature
    
    TFsol0.assign(TFsol)
    t   += dt
    tn  += 1
    print('\n')

#################################################################################################
## END ##########################################################################################
#################################################################################################