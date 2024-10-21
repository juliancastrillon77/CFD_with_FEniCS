# Julian Castrillon
# CFD - Turbulent model K-E 

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

# Fluid Properties
mu  = fe.Constant(1.0/5.0)  # Dynamic viscosity  [kg/ms]
rho = fe.Constant(1.0)        # Density            [kg/m3]
b   = fe.Constant((0.0, 0.0)) # Body accelerations [m/s2]

# Turbulent constants by Launder & Sharma in 1974
Cu     = fe.Constant(0.09)
Ce1    = fe.Constant(1.44)
Ce2    = fe.Constant(1.92)
sigmak = fe.Constant(1.0)
sigmae = fe.Constant(1.3)

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
K    = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
Epsi = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)

M    = fe.MixedElement([Vel, Pres, K, Epsi])
FS   = fe.FunctionSpace(Mesh, M)

TF   = fe.TrialFunction(FS)
(w, q, y, z) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
(v, p, k, e) = fe.split(TFsol)

StrTnsr = fe.Constant(0.5)*(fe.grad(v)+fe.grad(v).T)
mut     = Cu*(k**2)/e
P       = fe.Constant(2.0)*mut*fe.inner(StrTnsr,StrTnsr)

RANSeqn =   rho*fe.inner(w,fe.grad(v)*v)*fe.dx                              \
           - fe.div(w)*p*fe.dx                                              \
           - rho*fe.dot(w,b)*fe.dx                                          \
           + mu*fe.inner(fe.grad(w),fe.grad(v))*fe.dx                       \
           + mu*fe.inner(fe.grad(w),fe.grad(v).T)*fe.dx                     \
           + mut*fe.inner(fe.grad(w),fe.grad(v))*fe.dx                      \
           + mut*fe.inner(fe.grad(w),fe.grad(v).T)*fe.dx                    \
           + fe.inner(w,fe.grad(k))*fe.Constant(2.0)/fe.Constant(3.0)*fe.dx     \
           - q*fe.div(v)*fe.dx

#            + w*fe.Constant(2.0)/fe.Constant(3.0)*k*fe.Identity(2)*fe.dx \

keqn =  y*fe.dot(fe.grad(k),v)*fe.dx                            \
      - fe.inner(fe.grad(y),((mu+mut/sigmak)*fe.grad(k)))*fe.dx \
      - y*P*fe.dx                                               \
      + y*e*fe.dx
      
eeqn =  z*fe.dot(fe.grad(e),v)*fe.dx                            \
      - fe.inner(fe.grad(z),((mu+mut/sigmae)*fe.grad(e)))*fe.dx \
      - z*Ce1*P*e/k*fe.dx                                       \
      + z*Ce2*(e**2)/k*fe.dx      

WeakForm = RANSeqn + keqn + eeqn 
# Order for both equations:
# Advection
# Diffusion
# Production
# Dissipation

#vnorm = fe.sqrt(fe.dot(v,v))
#h = fe.CellDiameter(Mesh)
#tau = ((2.0*vnorm/h)**2 + ((4.0*mu)/h**2)**2)**(-0.5)
#R = fe.grad(v)*v+fe.grad(p)-mu*fe.div(fe.grad(v))-b
#SUPG = tau*fe.inner(fe.grad(w)*v,R)*fe.dx(metadata={'quadrature_degree':4})
#WeakForm += SUPG

### Petrov Galerkin Pressure Stabilzation (PSPG) stabilization for pressure field // The Ladyzhenskaya-Babuska-Brezzi condition not met
#PSPG = -tau*fe.inner(fe.grad(q),R)*fe.dx(metadata={'quadrature_degree':4})
#WeakForm += PSPG

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
InletFlow = fe.Constant((1, 0))
KIn       = fe.Constant(0.9)
EIn       = fe.Constant(0.0017)      

EntryBC          = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
BottomWallBC     = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
ExitBC           = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
TopWallBC        = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)
CircleBC         = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, Circle)
TriangleLeftBC   = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleLeft)
TriangleRightBC  = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleRight)

EntryBCK         = fe.DirichletBC(FS.sub(2), KIn,       DomainBoundaries, Entry)
TopWallBCK       = fe.DirichletBC(FS.sub(2), POut,      DomainBoundaries, TopWall)
BottomWallBCK    = fe.DirichletBC(FS.sub(2), POut,      DomainBoundaries, BottomWall)
CircleBCK        = fe.DirichletBC(FS.sub(2), POut,      DomainBoundaries, Circle)
TriangleLeftBCK  = fe.DirichletBC(FS.sub(2), POut,      DomainBoundaries, TriangleLeft)
TriangleRightBCK = fe.DirichletBC(FS.sub(2), POut,      DomainBoundaries, TriangleRight)

EntryBCE         = fe.DirichletBC(FS.sub(3), EIn,       DomainBoundaries, Entry)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, CircleBC, TriangleLeftBC, TriangleRightBC, EntryBCK, EntryBCE, \
       BottomWallBCK, TopWallBCK, CircleBCK, TriangleLeftBCK, TriangleRightBCK]

InitialVel  = fe.interpolate(fe.Constant((fe.Constant(1), fe.Constant(0))), FS.sub(0).collapse()) 
InitialPres = fe.interpolate(fe.Constant(10), FS.sub(1).collapse())
InitialK = fe.interpolate(KIn, FS.sub(2).collapse())
InitialE = fe.interpolate(EIn, FS.sub(3).collapse())

fe.assign(TFsol.sub(0),  InitialVel)
fe.assign(TFsol.sub(1),  InitialPres)
fe.assign(TFsol.sub(2),  InitialK)
fe.assign(TFsol.sub(3),  InitialE)

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['linear_solver'] = 'petsc'
Parameters['newton_solver']['absolute_tolerance']   = 1e-7
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['maximum_iterations']   = 5
Parameters['newton_solver']['relaxation_parameter'] = 1.0

Solver.solve()

(Velocity, Pressure, KineticEnergy, Dissipation) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')
KineticEnergy.rename('KineticEnergy','KineticEnergy')
Dissipation.rename('Dissipation','Dissipation')

FileVel  = fe.File('Results/Velocity.pvd')
FilePres = fe.File('Results/Pressure.pvd')
FileK    = fe.File('Results/K.pvd')
FileEpsi = fe.File('Results/Epsilon.pvd')
FileVel  << Velocity
FilePres << Pressure
FileK    << KineticEnergy
FileEpsi << Dissipation