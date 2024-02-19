# Julian Castrillon
# CFD - 2D Unsteady Incompressible Navier Stokes Equations

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
Mesh = fe.Mesh('../Mesh/2D Mesh/2DMesh.xml') # Create the mesh and import mesh into the solver

dt    = 0.01 # Timestep             [s]
t_end = 0.5  # Length of simulation [s]

vi = 1 # Inlet Velocity [m/s]
# Fluid Properties
mu = fe.Constant(1/500)  # Dynamic viscosity  [kg/ms]
b  = fe.Constant((0, 0)) # Body accelerations [m/s2]

theta = fe.Constant(0.5)      # Crank-Nicholson time-stepping scheme
n     = fe.FacetNormal(Mesh)  # Normal vector  
h     = fe.CellDiameter(Mesh) # Mesh size

# Define the mixed vector function space operating on this meshed domain
Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres])
FS   = fe.FunctionSpace(Mesh, M) # Function space

TF     = fe.TrialFunction(FS) # Trial function
(w, q) = fe.TestFunctions(FS) # Test  functions

#################################################################################################
## Weak formulation #############################################################################
TFsol  = fe.Function(FS) # Timestep n+1
(v, p) = fe.split(TFsol)

WeakFormB =  fe.inner(w,fe.grad(v)*v)*fe.dx             \
            + mu*fe.inner(fe.grad(w), fe.grad(v))*fe.dx \
            - p*fe.div(w)*fe.dx                         \
            - q*fe.div(v)*fe.dx                         \
            - fe.dot(w,b)*fe.dx                         \

TFsol0   = fe.Function(FS) # Timestep n
(v0, p0) = fe.split(TFsol0)

WeakFormF =   fe.inner(w,fe.grad(v0)*v0)*fe.dx          \
            + mu*fe.inner(fe.grad(w),fe.grad(v0))*fe.dx \
            - p*fe.div(w)*fe.dx                         \
            - q*fe.div(v0)*fe.dx                        \
            - fe.dot(w,b)*fe.dx                         \
                
WeakForm = fe.inner((v-v0)/dt,w)*fe.dx + (theta)*WeakFormB + (1-theta)*WeakFormF

### Streamwise Upwind Petrov Galerkin (SUPG) stabilization for convection #######################
# vnorm = fe.sqrt(fe.dot(u0,u0))
# tau   = ((1/(theta*dt))**2+(2*vnorm/h)**2+9*(4*mu/h**2)**2)**(-0.5)

# Residual of the strong form of Navier-Stokes and continuity
#R =   idt*(u-u0)             \
#    + theta *                \
#    ( fe.grad(u)*u           \
#    - mu*fe.div(fe.grad(u))  \
#    + fe.grad(p)             \
#    - b)                     \
#                             \
#    + (1-theta) *            \
#    ( fe.grad(u0)*u0         \
#    - mu*fe.div(fe.grad(u0)) \
#    + fe.grad(p)             \
#    - b)

# WeakForm += tau*fe.inner(fe.grad(w)*u0,R)*fe.dx(metadata={'quadrature_degree':4})

# Petrov Galerkin Pressure Stabilzation (PSPG) stabilization for pressure field // The Ladyzhenskaya-Babuska-Brezzi condition not met
# WeakForm += -tau*fe.inner(fe.grad(q),R)*fe.dx(metadata={'quadrature_degree':4})

J = fe.derivative(WeakForm, TFsol, TF) # Jacobian

#################################################################################################
## Boundary Conditions ##########################################################################
DomainBoundaries = fe.MeshFunction('size_t', Mesh, '../Mesh/2D Mesh/2DMesh_facet_region.xml')

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
InletFlow = fe.Constant((vi, 0))

EntryBC         = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
BottomWallBC    = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
ExitBC          = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
TopWallBC       = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)
CircleBC        = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, Circle)
TriangleLeftBC  = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleLeft)
TriangleRightBC = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleRight)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, CircleBC, TriangleLeftBC, TriangleRightBC]

## Boundary check
# import sys
# fe.File('New.pvd') << DomainBoundaries
# sys.exit()

#################################################################################################
## Solver #######################################################################################
Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['maximum_iterations']   = 4
Parameters['newton_solver']['relaxation_parameter'] = 1.0

# Create files for storing solution
FileVel  = fe.File('Results/Velocity.pvd')
FilePres = fe.File('Results/Pressure.pvd')

# Loops #####################################
# Inner loop - Write data so many time steps
# Outer loop - Time steps

t  = dt
tn = 0

while t < t_end:

    print("t = ", np.around(t,3))
    print("Solving ....")
    
    Solver.solve()
    
    (Velocity, Pressure) = TFsol.split(deepcopy = True)
    
    if tn%1 == 0:   # Save to file only if the time step increases by a step of 1 
        Velocity.rename('Velocity', 'Velocity')
        Pressure.rename('Pressure', 'Pressure')
        FileVel  << Velocity
        FilePres << Pressure
        print('Written velocity & pressure data')
    
    TFsol0.assign(TFsol)
    t   += dt
    tn  += 1
    print('\n')

#################################################################################################
## END ##########################################################################################
#################################################################################################