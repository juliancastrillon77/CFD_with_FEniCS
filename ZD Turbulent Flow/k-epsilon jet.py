# Julian Castrillon
# Turbulent k-epsilon 

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
mu  = fe.Constant(1.0/1.0)    # Dynamic viscosity  [kg/ms]
rho = fe.Constant(1.0)        # Density            [kg/m3]
b   = fe.Constant((0.0, 0.0)) # Body accelerations [m/s2]

Vel   = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres  = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
KiEng = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
Epsi  = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)

M    = fe.MixedElement([Vel, Pres, KiEng, Epsi])
FS   = fe.FunctionSpace(Mesh, M)

TF   = fe.TrialFunction(FS)
(w, q, y, z) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
(u, p, k, e) = fe.split(TFsol)

# Turbulent constants by Launder & Sharma in 1974
Cu     = fe.Constant(0.09)
Ce1    = fe.Constant(1.44)
Ce2    = fe.Constant(1.92)
sigmak = fe.Constant(1.0)
sigmae = fe.Constant(1.3)

StrTnsr = 0.5*(fe.grad(u)+fe.grad(u).T)     # Strain rate tensor
mut     = Cu*(k**2)/e                       # Turbulent viscosity
P       = 2.0*mut*fe.inner(StrTnsr,StrTnsr) # Production

dx = fe.dx(metadata={"quadrature_degree":4})

RANSeqn =   rho*fe.inner(w,fe.grad(u)*u)*dx           \
           - fe.div(w)*p*dx                           \
           - rho*fe.dot(w,b)*fe.dx                    \
           + mu*fe.inner(fe.grad(w),fe.grad(u))*dx    \
           - mut*fe.inner(fe.grad(w),fe.grad(u))*dx   \
           - mut*fe.inner(fe.grad(w),fe.grad(u).T)*dx \
           - fe.dot(w,fe.grad(k))*(2/3)*dx

Contieqn =  q*fe.div(u)*dx

keqn =   y*fe.dot(u,fe.grad(k))*dx                            \
       + fe.inner(fe.grad(y),((mu+mut/sigmak)*fe.grad(k)))*dx \
       - y*P*dx                                               \
       + y*e*dx
      
eeqn =   z*fe.dot(u,fe.grad(e))*fe.dx                         \
       + fe.inner(fe.grad(z),((mu+mut/sigmae)*fe.grad(e)))*dx \
       - z*Ce1*P*e/k*dx                                       \
       + z*Ce2*(e**2)/k*dx          

WeakForm = RANSeqn + Contieqn + keqn + eeqn 

J = fe.derivative(WeakForm, TFsol, TF)  

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DMeshJet_facet_region.xml')

EntryTop      = 8
Nozzle        = 12
EntryBottom   = 13
BottomWall    = 9
Exit          = 10
TopWall       = 11

I = 0.05 # Turbulent intensity
KInnum = 3/2*(10*I)**2
EInnum = Cu**(3/4)*KInnum**(3/2)/(0.07*1)

POut             = fe.Constant(0)
InletFlowTop     = fe.Constant((0, 0))

InletFlowNozzle  = fe.Constant((15, 0.5))
InletFlowBottom  = fe.Constant((0, 0))
KIn              = fe.Constant(KInnum)
EIn              = fe.Constant(EInnum)  

EntryTopBC       = fe.DirichletBC(FS.sub(0), InletFlowTop,    DomainBoundaries, EntryTop)
NozzleBC         = fe.DirichletBC(FS.sub(0), InletFlowNozzle, DomainBoundaries, Nozzle)
EntryBottomBC    = fe.DirichletBC(FS.sub(0), InletFlowBottom, DomainBoundaries, EntryBottom)
ExitBC           = fe.DirichletBC(FS.sub(1), POut,            DomainBoundaries, Exit)

NozzleBCK        = fe.DirichletBC(FS.sub(2), KIn,             DomainBoundaries, Nozzle)
NozzleBCE        = fe.DirichletBC(FS.sub(3), EIn,             DomainBoundaries, Nozzle)

BCs = [EntryTopBC, NozzleBC, EntryBottomBC, ExitBC, NozzleBCK, NozzleBCE]

InitialVel  = fe.interpolate(fe.Constant((fe.Constant(10), fe.Constant(1))), FS.sub(0).collapse()) 
InitialPres = fe.interpolate(fe.Constant(1), FS.sub(1).collapse())
InitialK    = fe.interpolate(fe.Constant(0.01), FS.sub(2).collapse())
InitialE    = fe.interpolate(fe.Constant(0.5), FS.sub(3).collapse())

fe.assign(TFsol.sub(0), InitialVel)
fe.assign(TFsol.sub(1), InitialPres)
fe.assign(TFsol.sub(2), InitialK)
fe.assign(TFsol.sub(3), InitialE)

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['absolute_tolerance']   = 1e-6
Parameters['newton_solver']['relative_tolerance']   = 1e-6
Parameters['newton_solver']['maximum_iterations']   = 13
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

FileVel  = fe.File('ResultsJetKepsilon/Velocity.pvd')
FilePres = fe.File('ResultsJetKepsilon/Pressure.pvd')
FileK    = fe.File('ResultsJetKepsilon//K.pvd')
FileEpsi = fe.File('ResultsJetKepsilon//Epsilon.pvd')
FileVel  << Velocity
FilePres << Pressure
FileK    << KineticEnergy
FileEpsi << Dissipation