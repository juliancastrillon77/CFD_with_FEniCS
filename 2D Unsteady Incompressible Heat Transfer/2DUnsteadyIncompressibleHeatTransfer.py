# Julian Castrillon
# 2D Unsteady Incompressible Heat Transfer

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
Mesh = fe.Mesh('Mesh/2D Mesh/2DMesh.xml') # Create the mesh and import mesh into the solver

dt    = 0.01 # Timestep             [s]
t_end = 10   # Length of simulation [s]

vi = 0.5 # Inlet Velocity [m/s]
T1 = 340 # Temperature 1  [K]
T2 = 273 # Temperature 2  [K] 

## Air Properties at 25ºC 1atm
#mu  = fe.Constant(0.00001825) # Dynamic viscosity              [kg/ms]
mu  = fe.Constant(0.01825)    # Dynamic viscosity               [kg/ms]
rho = fe.Constant(1.225)      # Density                         [kg/m3]
k   = fe.Constant(0.02551)    # Thermal conductivity            [W/mK]
cp  = fe.Constant(1007)       # Heat capacity constant pressure [J/kgK]
cv  = fe.Constant(718)        # Heat capacity constant volume   [J/kgK]
RR  = fe.Constant(287)        # Gas constant                    [J/kgK]
b   = fe.Constant((0, 0))     # Body accelerations              [m/s2]

theta = fe.Constant(0.5)      # Crank-Nicholson time-stepping scheme
n     = fe.FacetNormal(Mesh)  # Normal vector  
h     = fe.CellDiameter(Mesh) # Mesh size

# Define the mixed vector function space operating on this meshed domain
Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
Temp = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres, Temp])
FS   = fe.FunctionSpace(Mesh, M) # Function space

TF        = fe.TrialFunction(FS) # Trial function
(w, q, s) = fe.TestFunctions(FS) # Test  functions

TFsol     = fe.Function(FS) # Timestep n+1
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

WeakFormB = Mnt + Cnt + Egy

TFsol0       = fe.Function(FS) # Timestep n
(u0, p0, T0) = fe.split(TFsol0)

S0 = 0.5*(fe.grad(u0)+fe.grad(u0).T)      # Strain rate tensor
Cst0 = -(p*I)+(2*mu*S0)                   # Cauchy stress tensor

Mnt0 = fe.inner(w,fe.grad(rho*u0)*u0)*dx \
     - fe.dot(w,(rho*b))*dx            \
     + fe.inner(fe.grad(w),Cst0)*dx

Cnt0 = q*fe.div(rho*u0)*dx

Egy0 = s*rho*cp*fe.dot(u0,fe.grad(T))*dx    \
     + fe.dot(fe.grad(s), k*fe.grad(T))*dx 

WeakFormF = Mnt0 + Cnt0 + Egy0
          
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
            
R = (u-u0)/dt                \
    + theta *                \
    (fe.grad(u)*u            \
    + fe.grad(p)             \
    - rho*b                  \
    - mu*fe.div(fe.grad(u))) \
    \
    + (1-theta) *            \
    (fe.grad(u0)*u0          \
    + fe.grad(p)             \
    - rho*b                  \
    - mu*fe.div(fe.grad(u0)))# \
#\
#\
#+ q*fe.div(u)*fe.dx \
#\
#\
#- s*fe.dot(u,fe.grad(P))*fe.dx                      \
#+ (cp/(k))*rho*fe.dot(s,fe.dot(u,fe.grad(T)))*fe.dx \
#+ fe.dot(fe.grad(s),fe.grad(T))*fe.dx               \

#    + (1-theta) *            \
#    ( fe.grad(u0)*u0         \
#    - mu*fe.div(fe.grad(u0)) \
#    + fe.grad(P)             \
#    - b)

#WeakForm += tau*fe.inner(fe.grad(w)*u0,R)*fe.dx(metadata={'quadrature_degree':8})

# Petrov Galerkin Pressure Stabilzation (PSPG) stabilization for pressure field // The Ladyzhenskaya-Babuska-Brezzi condition not met
# WeakForm += -tau*fe.inner(fe.grad(s),R)*fe.dx(metadata={'quadrature_degree':8})

J = fe.derivative(WeakForm, TFsol, TF) # Jacobian

#################################################################################################
## Boundary Conditions ##########################################################################
DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DMesh_facet_region.xml')

# Identification of all correct boundary markers needed for the domain
Entry         = 8
BottomWall    = 9
Exit          = 10
TopWall       = 11
Circle        = 12
TriangleLeft  = 13
TriangleRight = 14

NoSlip    = fe.Constant((0, 0))
POut      = fe.Constant(0)
InletFlow = fe.Expression(['6*x[1]*(1-x[1])','0'], degree=2)
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

fe.assign(TFsol.sub(0),  InitialVel)
fe.assign(TFsol.sub(1),  InitialPres)
fe.assign(TFsol.sub(2),  InitialTemp)
fe.assign(TFsol0.sub(0), InitialVel)
fe.assign(TFsol0.sub(1), InitialPres)
fe.assign(TFsol0.sub(2), InitialTemp)

#################################################################################################
## Solver #######################################################################################
Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['linear_solver'] = 'mumps'
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['maximum_iterations']   = 7
Parameters['newton_solver']['relaxation_parameter'] = 1.0

# Create files for storing solution
FileVel  = fe.File('Results/Velocity.pvd')
FilePres = fe.File('Results/Pressure.pvd')
FileTemp = fe.File('Results/Temperature.pvd')

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
        FileVel  << Velocity
        FilePres << Pressure
        FileTemp << Temperature
        print('Written velocity, pressure & temperature data')
        
    TFsol0.assign(TFsol)
    t   += dt
    tn  += 1
    print('\n')

#################################################################################################
## END ##########################################################################################
#################################################################################################