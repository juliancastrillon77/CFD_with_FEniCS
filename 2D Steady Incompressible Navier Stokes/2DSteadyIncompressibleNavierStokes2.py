# Julian Castrillon
# CFD - 2D Steady Incompressible Navier Stokes Equations

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

Mesh = fe.Mesh('Mesh/2D Mesh/2DMeshJet.xml')

# Fluid Properties
mu  = fe.Constant(1/500)   # Dynamic viscosity  [kg/ms]
rho = fe.Constant(1)      # Density            [kg/m3]
b   = fe.Constant((0, 0)) # Body accelerations [m/s2]

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres])
FS   = fe.FunctionSpace(Mesh, M)

TF     = fe.TrialFunction(FS)
(w, q) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
(v, p) = fe.split(TFsol)

WeakForm =   rho*fe.inner(w,fe.grad(v)*v)*fe.dx       \
           - fe.div(w)*p*fe.dx                        \
           - rho*fe.dot(w,b)*fe.dx                    \
           + mu*fe.inner(fe.grad(w),fe.grad(v))*fe.dx \
           - q*fe.div(v)*fe.dx

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

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DMeshJet_facet_region.xml')

EntryTop      = 8
Nozzle        = 12
EntryBottom   = 13
BottomWall    = 9
Exit          = 10
TopWall       = 11

POut             = fe.Constant(0)
NoSlip           = fe.Constant((0, 0))
InletFlowTop     = fe.Constant((0, 0))
InletFlowNozzle  = fe.Constant((10, 0))
InletFlowBottom  = fe.Constant((0, 0))
nuhatOut         = fe.Constant(3*mu)

EntryTopBC       = fe.DirichletBC(FS.sub(0), InletFlowTop,    DomainBoundaries, EntryTop)
NozzleBC         = fe.DirichletBC(FS.sub(0), InletFlowNozzle, DomainBoundaries, Nozzle)
EntryBottomBC    = fe.DirichletBC(FS.sub(0), InletFlowBottom, DomainBoundaries, EntryBottom)
ExitBC           = fe.DirichletBC(FS.sub(1), POut,            DomainBoundaries, Exit)

TopWallBC        = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)
BottomWallBC     = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)

#TopWallBCK       = fe.DirichletBC(FS.sub(2), nuhatIn,         DomainBoundaries, TopWall)

#BCs = [EntryTopBC, NozzleBC, EntryBottomBC, ExitBC, Exit2BC, Exit3BC, EntryTopBCK, NozzleBCK, EntryBottomBCK,  BottomWallBCK, \
#       ExitBCK, TopWallBCK]

BCs = [EntryTopBC, NozzleBC, EntryBottomBC, ExitBC, TopWallBC, BottomWallBC]

InitialVel  = fe.interpolate(fe.Constant((fe.Constant(1), fe.Constant(0))), FS.sub(0).collapse()) 
InitialPres = fe.interpolate(fe.Constant(10), FS.sub(1).collapse())

fe.assign(TFsol.sub(0),  InitialVel)
fe.assign(TFsol.sub(1),  InitialPres)

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
#Parameters['newton_solver']['linear_solver'] = 'cg'
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['maximum_iterations']   = 7
Parameters['newton_solver']['relaxation_parameter'] = 1.0
Parameters['newton_solver']['error_on_nonconvergence'] = False

Solver.solve()

(Velocity, Pressure) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')

FileVel  = fe.File('Results3/Velocity.pvd')
FilePres = fe.File('Results3/Pressure.pvd')
FileVel  << Velocity
FilePres << Pressure