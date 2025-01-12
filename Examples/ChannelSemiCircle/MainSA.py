# Julian Castrillon
# Turbulent Spalart-Allmaras

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

Mesh = fe.Mesh('Mesh/2D Mesh/2DChannelSemiCircle.xml')

# Fluid Properties
mu  = fe.Constant(0.1)      # Dynamic viscosity  [kg/ms]
rho = fe.Constant(10)       # Density            [kg/m3]
b   = fe.Constant((0.0, 0.0)) # Body accelerations [m/s2]

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
#d = fe.Constant(0.05)

dspace = fe.FunctionSpace(Mesh, 'Lagrange', 1)
d = fe.Function(dspace)
fe.File('Results/Distance.xml') >> d

X       = nuhat/mu
fv1     = (X**3)/((X**3)+(cv1**3))
nutr    = fv1*nuhat
fv2     = 1-(X/(1+(X*fv1)))
RotTnsr = 0.5*(fe.grad(u)-fe.grad(u).T)
Omega   = fe.sqrt(2*fe.inner(RotTnsr,RotTnsr))
#Shat    = Omega+fv2*nuhat/((k**2)*(d**2))
Shat    = fe.conditional(fe.le(Omega+fv2*nuhat/((k**2)*(d**2)),0.3*Omega), 0.3*Omega, Omega+fv2*nuhat/((k**2)*(d**2)))
#r       = nuhat/(Shat*(k**2)*(d**2))
r       = fe.conditional(fe.le(nuhat/(Shat*(k**2)*(d**2)),10),10,nuhat/(Shat*(k**2)*(d**2)))
ft2     = ct3*fe.exp(-ct4*(X**2))
g       = r+cw2*((r**6)-r)
fw      = g*((1+(cw3**6))/((g**6)+(cw3**6)))**(1.0/6.0)


S = 0.5*(fe.grad(u)+fe.grad(u).T)      # Strain rate tensor
I = fe.Identity(Mesh.topology().dim()) # Identity tensor
Cst = -(p)*I + 2*(mu)*S                # Cauchy stress tensor

dx = fe.dx(metadata={"quadrature_degree":4})

RMnt = fe.inner(w,fe.grad(rho*u)*u)*dx    \
     - fe.dot(w,(rho*b))*dx               \
     + fe.inner(fe.grad(w),Cst)*dx        \
     + 2*fe.inner(fe.grad(w),(nutr*S))*dx

Cnt = q*fe.div(rho*u)*dx

MVisc = y*fe.dot(u,fe.grad(nuhat))*dx                                 \
      - y*(cb1*(1-ft2)*Shat*nuhat)*dx                                 \
      + y*(((cw1*fw)-(cb1/(k**2)))*(nuhat/d)**2)*dx                   \
      + (1/sigma)*fe.dot(fe.grad(y),(nuhat+mu/rho)*fe.grad(nuhat))*dx \
      - y*(cb2/sigma)*fe.dot(fe.grad(nuhat),fe.grad(nuhat))*dx

WeakForm = RMnt + Cnt + MVisc

J = fe.derivative(WeakForm, TFsol, TF)  

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DChannelSemiCircle_facet_region.xml')

Entry      = 1
BottomWall = 2
Exit       = 3
TopWall    = 4

NoSlip    = fe.Constant((0, 0))
POut      = fe.Constant(0)
#InletFlow = fe.Expression(['2*((4*x[1]*(0.1-x[1]))/(0.1*0.1))','0'], degree=2)
InletFlow = fe.Constant((20, 0))
nuhatOut  = fe.Constant(3*mu)

EntryBC           = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
BottomWallBC      = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
ExitBC            = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
TopWallBC         = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)

TopWallBCnuhat    = fe.DirichletBC(FS.sub(2), POut,      DomainBoundaries, TopWall)
BottomWallBCnuhat = fe.DirichletBC(FS.sub(2), POut,      DomainBoundaries, BottomWall)

EntryBCnuhat      = fe.DirichletBC(FS.sub(2), nuhatOut,  DomainBoundaries, Entry)
ExitBCnuhat       = fe.DirichletBC(FS.sub(2), nuhatOut,  DomainBoundaries, Exit)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, TopWallBCnuhat, BottomWallBCnuhat, EntryBCnuhat, ExitBCnuhat]

InitialVel   = fe.interpolate(fe.Constant((fe.Constant(0.01), fe.Constant(0.01))), FS.sub(0).collapse()) 
InitialPres  = fe.interpolate(fe.Constant(0.01), FS.sub(1).collapse())
Initialnuhat = fe.interpolate(fe.Constant(0.01), FS.sub(2).collapse())

fe.assign(TFsol.sub(0), InitialVel)
fe.assign(TFsol.sub(1), InitialPres)
fe.assign(TFsol.sub(2), Initialnuhat)

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver'][ "linear_solver"] = "mumps"
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['maximum_iterations']   = 1
Parameters['newton_solver']['relaxation_parameter'] = 0.5
Parameters['newton_solver']['error_on_nonconvergence'] = False

URFu = 0.9
URFp = 0.9
URFnuhat = 0.9

TFsol_prev = fe.Function(FS)
(u_prev, p_prev, nuhat_prev) = TFsol_prev.split()
TFsol_relaxed = fe.Function(FS)
(u_relaxed, p_relaxed, nuhat_relaxed) = TFsol_relaxed.split()

Total_iterations = 8
for i in range(Total_iterations):
      print(f"\nIteration {i+1}: \n")
      Solver.solve()
      
      (u_star, p_star, nuhat_star) = TFsol.split()
      
      u_relaxed.vector()[:] = URFu*u_star.vector()[:] + (1-URFu)*u_prev.vector()[:]
      p_relaxed.vector()[:] = URFp*p_star.vector()[:] + (1-URFp)*p_prev.vector()[:]
      nuhat_relaxed.vector()[:] = URFnuhat*nuhat_star.vector()[:] + (1-URFp)*nuhat_prev.vector()[:]

      fe.assign(TFsol.sub(0), u_relaxed)
      fe.assign(TFsol.sub(1), p_relaxed)
      fe.assign(TFsol.sub(1), nuhat_relaxed)

      u_prev = u_relaxed
      p_prev = p_relaxed
      nuhat_prev = nuhat_relaxed

      print(g)

      if i == Total_iterations - 1:
            print("Reached maximum limit of iterations")
            break

FS2 = fe.FunctionSpace(Mesh, 'Lagrange', 1)

L        = fe.Constant(0.1)
uMag     = fe.sqrt(fe.dot(u,u))
Re       = (rho*uMag*L)/mu
Reynolds = fe.project(Re, FS2)

ShearTnsr    = (2*mu*S) 
ShearStressV = fe.sqrt(fe.inner(ShearTnsr,ShearTnsr))
ShearStress  = fe.project(ShearStressV, FS2)

h      = fe.CellDiameter(Mesh)
Pe     = uMag*h/(2*mu)
Peclet = fe.project(Pe, FS2)

(Velocity, Pressure, ViscosityHat) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')
ViscosityHat.rename('KineticEnergy','KineticEnergy')
Reynolds.rename('Reynolds','Reynolds')
ShearStress.rename('ShearStress','ShearStress')
Peclet.rename('Peclet','Peclet')
DomainBoundaries.rename('DomainBoundaries','DomainBoundaries')

fe.File('Results Re 10t6/Peclet.pvd')       << Peclet
fe.File('Results Re 10t6/Velocity.xml')     << Velocity
fe.File('Results Re 10t6/Pressure.xml')     << Pressure
fe.File('Results Re 10t6/Velocity.pvd')     << Velocity
fe.File('Results Re 10t6/Pressure.pvd')     << Pressure
fe.File('Results Re 10t6/ViscosityHat.pvd') << ViscosityHat
fe.File('Results Re 10t6/Reynolds.xml')     << Reynolds
fe.File('Results Re 10t6/ShearStress.xml')  << ShearStress
fe.File('Results Re 10t6/Velocity.pvd')     << Velocity
fe.File('Results Re 10t6/Pressure.pvd')     << Pressure
fe.File('Results Re 10t6/Reynolds.pvd')     << Reynolds
fe.File('Results Re 10t6/ShearStress.pvd')  << ShearStress
fe.File('Results Re 10t6/Boundaries.pvd')   << DomainBoundaries






