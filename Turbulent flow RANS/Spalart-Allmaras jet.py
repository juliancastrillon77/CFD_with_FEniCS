# Julian Castrillon
# CFD - Turbulent Spalart Allmaras

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
VisHt = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)

M    = fe.MixedElement([Vel, Pres, VisHt])
FS   = fe.FunctionSpace(Mesh, M)

TF   = fe.TrialFunction(FS)
(w, q, y) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
(v, p, nuhat) = fe.split(TFsol)

# Model constants 
cb1   = fe.Constant(0.1355)
cb2   = fe.Constant(0.622)
cw2   = fe.Constant(0.3)
cw3   = fe.Constant(2.0)
cv1   = fe.Constant(7.1)
ct3   = fe.Constant(1.2)
ct4   = fe.Constant(0.5)
sigma = fe.Constant(2/3)
k     = fe.Constant(0.41)
cw1   = (cb1/(k**2))+((1.0+cb2)/(sigma))

#class distance_to_wall(fe.UserExpression):

#    def eval(self, value, x):
#        if x[1] > 1/2:
#            value[0] = 1 - x[1]
#        else:
#            value[0] = x[1]
#        return(value)
#d = distance_to_wall()
#W = FS.sub(2)
#d = fe.interpolate(d, W.collapse())

# Model equations
d       = fe.Constant(10)
X       = nuhat/mu
fv1     = (X**3)/((X**3)+(cv1**3))
nutr    = fv1*nuhat
fv2     = 1-(X/(1+(X*fv1)))
RotTnsr = 0.5*(fe.grad(v)-fe.grad(v).T)
Omega   = fe.sqrt(2*fe.inner(RotTnsr,RotTnsr))
S       = Omega+fv2*nuhat/((k**2)*(d**2))
S       = fe.conditional(fe.le(S,0), 0.3*Omega, S)
r       = fe.conditional(fe.le(nuhat/(S*(k**2)*(d**2)),10),nuhat/(S*(k**2)*(d**2)), 10)
ft2     = ct3*fe.exp(-ct4*(X**2))
g       = r+cw2*((r**6)-r)
fw      = g*((1+(cw3**6))/((g**6)+(cw3**6)))**(1.0/6.0)

dx = fe.dx(metadata={"quadrature_degree":4})
StrTnsr = fe.Constant(0.5)*(fe.grad(v)+fe.grad(v).T)

RANSeqn =   rho*fe.inner(w,fe.grad(v)*v)*dx               \
           - fe.div(w)*p*dx                               \
           + (mu+nutr)*fe.inner(fe.grad(w),fe.grad(v))*dx \
           - rho*fe.dot(w,b)*dx

Contieqn =  q*fe.div(v)*dx

MVisceqn =    y*fe.dot(v,fe.grad(nuhat))*dx                             \
            - y*(cb1*(1-ft2)*S*nuhat)*dx                                \
            + y*(((cw1*fw)-(cb1/(k**2)))*(nuhat/d)**2)*dx               \
            + (1/sigma)*fe.dot(fe.grad(y),(nuhat+mu)*fe.grad(nuhat))*dx \
            - y*(cb2/sigma)*fe.dot(fe.grad(nuhat),fe.grad(nuhat))*dx

WeakForm = RANSeqn + Contieqn + MVisceqn

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
InletFlowNozzle  = fe.Constant((200, 10))
InletFlowBottom  = fe.Constant((0, 0))
nuhatOut         = fe.Constant(3*mu)

EntryTopBC        = fe.DirichletBC(FS.sub(0), InletFlowTop,    DomainBoundaries, EntryTop)
NozzleBC          = fe.DirichletBC(FS.sub(0), InletFlowNozzle, DomainBoundaries, Nozzle)
EntryBottomBC     = fe.DirichletBC(FS.sub(0), InletFlowBottom, DomainBoundaries, EntryBottom)
ExitBC            = fe.DirichletBC(FS.sub(1), POut,            DomainBoundaries, Exit)

TopWallBC         = fe.DirichletBC(FS.sub(0), NoSlip,          DomainBoundaries, TopWall)
BottomWallBC      = fe.DirichletBC(FS.sub(0), NoSlip,          DomainBoundaries, BottomWall)
TopWallBCnuhat    = fe.DirichletBC(FS.sub(2), POut,            DomainBoundaries, TopWall)
BottomWallBCnuhat = fe.DirichletBC(FS.sub(2), POut,            DomainBoundaries, BottomWall)

ExitBCnuhat       = fe.DirichletBC(FS.sub(2), nuhatOut,        DomainBoundaries, Exit)
NozzleBCnuhat     = fe.DirichletBC(FS.sub(2), nuhatOut,        DomainBoundaries, Nozzle)

BCs = [EntryTopBC, NozzleBC, EntryBottomBC, ExitBC, ExitBCnuhat, NozzleBCnuhat, TopWallBC, BottomWallBC, TopWallBCnuhat, BottomWallBCnuhat]

InitialVel   = fe.interpolate(fe.Constant((fe.Constant(0.1), fe.Constant(0.1))), FS.sub(0).collapse()) 
InitialPres  = fe.interpolate(fe.Constant(0.1), FS.sub(1).collapse())
Initialnuhat = fe.interpolate(fe.Constant(0.1), FS.sub(2).collapse())

fe.assign(TFsol.sub(0),  InitialVel)
fe.assign(TFsol.sub(1),  InitialPres)
fe.assign(TFsol.sub(2),  Initialnuhat)

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
#Parameters['newton_solver']['linear_solver'] = 'petsc'
Parameters['newton_solver']['absolute_tolerance']   = 1e-7
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['maximum_iterations']   = 10
Parameters['newton_solver']['relaxation_parameter'] = 1.0
Parameters['newton_solver']['error_on_nonconvergence'] = False

Solver.solve()

(Velocity, Pressure, ViscosityHat) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')
ViscosityHat.rename('KineticEnergy','KineticEnergy')

FileVel   = fe.File('ResultsSA/Velocity.pvd')
FilePres  = fe.File('ResultsSA/Pressure.pvd')
FileVisHt = fe.File('ResultsSA/ViscosityHat.pvd')
FileVel   << Velocity
FilePres  << Pressure
FileVisHt << ViscosityHat