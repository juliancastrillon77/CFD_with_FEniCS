# Julian Castrillon
# Turbulence Spalart Allmaras

import os
import numpy as np
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

fe.parameters['form_compiler']['representation']     = 'uflacs'
fe.parameters['form_compiler']['optimize']           = True
fe.parameters['form_compiler']['cpp_optimize']       = True
fe.parameters["form_compiler"]["cpp_optimize_flags"] = '-O2 -funroll-loops'

Mesh = fe.Mesh('Mesh/2D Mesh/2DMeshJet.xml')

# Fluid Properties
mu  = fe.Constant(1.0/1.0)      # Dynamic viscosity  [kg/ms]
rho = fe.Constant(1.0)          # Density            [kg/m3]
b   = fe.Constant((0.0, 0.0))   # Body accelerations [m/s2]

Vel   = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres  = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
VisHt = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)

M    = fe.MixedElement([Vel, Pres, VisHt])
FS   = fe.FunctionSpace(Mesh, M)

TF   = fe.TrialFunction(FS)
(w, q, y) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
(u, p, nuhat) = fe.split(TFsol)

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

# Model equations
d       = fe.Constant(10)
X       = nuhat/mu
fv1     = (X**3)/((X**3)+(cv1**3))
nutr    = fv1*nuhat
fv2     = 1-(X/(1+(X*fv1)))
RotTnsr = 0.5*(fe.grad(u)-fe.grad(u).T)
Omega   = fe.sqrt(2*fe.inner(RotTnsr,RotTnsr))
Shat    = Omega+fv2*nuhat/((k**2)*(d**2))
Shat    = fe.conditional(fe.le(Shat,0.3*Omega), 0.3*Omega, Shat)
r       = nuhat/(Shat*(k**2)*(d**2))
#r       = fe.conditional(fe.le(Shat,0),10,nuhat/(Shat*(k**2)*(d**2)))
ft2     = ct3*fe.exp(-ct4*(X**2))
g       = r+cw2*((r**6)-r)
fw      = g*((1+(cw3**6))/((g**6)+(cw3**6)))**(1.0/6.0)

S = 0.5*(fe.grad(u)+fe.grad(u).T)      # Strain rate tensor
I = fe.Identity(Mesh.topology().dim()) # Identity tensor
Cst = -(p)*I + 2*(mu)*S                # Cauchy stress tensor

dx = fe.dx(metadata={"quadrature_degree":4})

RANS = fe.inner(w,fe.grad(rho*u)*u)*dx    \
     - fe.dot(w,(rho*b))*dx               \
     + fe.inner(fe.grad(w),Cst)*dx        \
     + 2*fe.inner(fe.grad(w),(nutr*S))*dx

Conti = q*fe.div(rho*u)*dx

MVisc = y*fe.dot(u,fe.grad(nuhat))*dx                                 \
      - y*(cb1*(1-ft2)*Shat*nuhat)*dx                                 \
      + y*(((cw1*fw)-(cb1/(k**2)))*(nuhat/d)**2)*dx                   \
      + (1/sigma)*fe.dot(fe.grad(y),(nuhat+mu/rho)*fe.grad(nuhat))*dx \
      - y*(cb2/sigma)*fe.dot(fe.grad(nuhat),fe.grad(nuhat))*dx

#WeakForm = RANS + Conti + MVisc
WeakForm = RANS + Conti + MVisc


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

BCs = [EntryTopBC, NozzleBC, EntryBottomBC, ExitBC, ExitBCnuhat, NozzleBCnuhat]

InitialVel   = fe.interpolate(fe.Constant((fe.Constant(0.01), fe.Constant(0.01))), FS.sub(0).collapse()) 
InitialPres  = fe.interpolate(fe.Constant(0.01), FS.sub(1).collapse())
Initialnuhat = fe.interpolate(fe.Constant(0.01), FS.sub(2).collapse())

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