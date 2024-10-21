# Julian Castrillon
# CFD - Turbulent model K-E 

import os
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

#fe.set_log_level(fe.LogLevel.PROGRESS)
fe.parameters['form_compiler']['representation']     = 'uflacs'
fe.parameters['form_compiler']['optimize']           = True
fe.parameters['form_compiler']['cpp_optimize']       = True
fe.parameters["form_compiler"]["cpp_optimize_flags"] = '-O2 -funroll-loops'

Mesh = fe.Mesh('Mesh/2D Mesh/2DMeshJet.xml')

# Fluid Properties
mu  = fe.Constant(1.0/10.0)  # Dynamic viscosity  [kg/ms]
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
           - mut*fe.inner(fe.grad(w),fe.grad(v))*fe.dx                      \
           - mut*fe.inner(fe.grad(w),fe.grad(v).T)*fe.dx                    \
           - fe.dot(w,fe.grad(k))*fe.Constant(2.0)/fe.Constant(3.0)*fe.dx
           
Contieqn =  q*fe.div(v)*fe.dx

keqn =   y*fe.dot(v,fe.grad(k))*fe.dx                            \
       + fe.inner(fe.grad(y),((mu+mut/sigmak)*fe.grad(k)))*fe.dx \
       - y*P*fe.dx                                               \
       + y*e*fe.dx
      
eeqn =   z*fe.dot(v,fe.grad(e))*fe.dx                            \
       + fe.inner(fe.grad(z),((mu+mut/sigmae)*fe.grad(e)))*fe.dx \
       - z*Ce1*P*e/k*fe.dx                                       \
       + z*Ce2*(e**2)/k*fe.dx          

WeakForm = RANSeqn + Contieqn + keqn + eeqn 
# Order for both equations:
# Advection
# Diffusion
# Production
# Dissipation

vnorm = fe.sqrt(fe.dot(v,v))
h = fe.CellDiameter(Mesh)
tau = ((2.0*vnorm/h)**2+((4.0*mu)/h**2)**2)**(-0.5)

R = fe.grad(v)*v+fe.grad(p)-mu*fe.div(fe.grad(v))-b+fe.div(-1*((fe.grad(v)+fe.grad(v).T)))

SUPG = tau*fe.inner(fe.grad(w)*v,R)*fe.dx
#WeakForm += SUPG

J = fe.derivative(WeakForm, TFsol, TF)  

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DMeshJet_facet_region.xml')

EntryTop      = 8
Nozzle        = 12
EntryBottom   = 13
BottomWall    = 9
Exit          = 10
TopWall       = 11

I = 0.005 # Turbulent intensity
Rratio = 10 # mut/mu

KInnum = 3/2*(2*I)**2
EInnum = Cu*KInnum**2/(Rratio*mu)

POut             = fe.Constant(0)
InletFlowTop     = fe.Constant((0, 0))
InletFlowNozzle  = fe.Constant((1, 0))
InletFlowBottom  = fe.Constant((0, 0))
KIn              = fe.Constant(KInnum)
EIn              = fe.Constant(EInnum)  

EntryTopBC       = fe.DirichletBC(FS.sub(0), InletFlowTop,    DomainBoundaries, EntryTop)
NozzleBC         = fe.DirichletBC(FS.sub(0), InletFlowNozzle, DomainBoundaries, Nozzle)
EntryBottomBC    = fe.DirichletBC(FS.sub(0), InletFlowBottom, DomainBoundaries, EntryBottom)
ExitBC           = fe.DirichletBC(FS.sub(1), POut,            DomainBoundaries, Exit)
Exit2BC          = fe.DirichletBC(FS.sub(1), POut,           DomainBoundaries, TopWall)
Exit3BC          = fe.DirichletBC(FS.sub(1), POut,           DomainBoundaries, BottomWall)

EntryTopBCK      = fe.DirichletBC(FS.sub(2), KIn,             DomainBoundaries, EntryTop)
NozzleBCK        = fe.DirichletBC(FS.sub(2), KIn,             DomainBoundaries, Nozzle)
EntryBottomBCK   = fe.DirichletBC(FS.sub(2), KIn,             DomainBoundaries, EntryBottom)
EntryTopBCE      = fe.DirichletBC(FS.sub(3), EIn,             DomainBoundaries, EntryTop)
NozzleBCE        = fe.DirichletBC(FS.sub(3), EIn,             DomainBoundaries, Nozzle)
EntryBottomBCE   = fe.DirichletBC(FS.sub(3), EIn,             DomainBoundaries, EntryBottom)

BCs = [EntryTopBC, NozzleBC, EntryBottomBC, ExitBC, Exit2BC, Exit3BC, EntryTopBCK, NozzleBCK, EntryBottomBCK,  EntryTopBCE, \
       NozzleBCE, EntryBottomBCE]

InitialVel  = fe.interpolate(fe.Constant((fe.Constant(11), fe.Constant(0.1))), FS.sub(0).collapse()) 
InitialPres = fe.interpolate(fe.Constant(100), FS.sub(1).collapse())
InitialK    = fe.interpolate(KIn, FS.sub(2).collapse())
InitialE    = fe.interpolate(EIn, FS.sub(3).collapse())

fe.assign(TFsol.sub(0),  InitialVel)
fe.assign(TFsol.sub(1),  InitialPres)
fe.assign(TFsol.sub(2),  InitialK)
fe.assign(TFsol.sub(3),  InitialE)

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['linear_solver'] = 'petsc'
Parameters['newton_solver']['absolute_tolerance']   = 1e-6
Parameters['newton_solver']['relative_tolerance']   = 1e-6
Parameters['newton_solver']['maximum_iterations']   = 4
Parameters['newton_solver']['relaxation_parameter'] = 1.0
Parameters['newton_solver']['error_on_nonconvergence'] = False

Total_iterations = 1
for i in range(Total_iterations):
      print(f"Iteration {i+1}:")
      Solver.solve()
      if i == Total_iterations - 1:
            print("Reached maximum limit of iterations")
            break

(Velocity, Pressure, KineticEnergy, Dissipation) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')
KineticEnergy.rename('KineticEnergy','KineticEnergy')
Dissipation.rename('Dissipation','Dissipation')

FileVel  = fe.File('ResultsJet3/Velocity.pvd')
FilePres = fe.File('ResultsJet3/Pressure.pvd')
FileK    = fe.File('ResultsJet3/K.pvd')
FileEpsi = fe.File('ResultsJet3/Epsilon.pvd')
FileVel  << Velocity
FilePres << Pressure
FileK    << KineticEnergy
FileEpsi << Dissipation