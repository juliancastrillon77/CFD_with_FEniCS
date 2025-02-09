# Julian Castrillon
# Channel heating

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

Mesh = fe.Mesh('Mesh/2D Mesh/2DCircuit.xml')

URFu = 1 # Under relaxation factor for velocity
URFp = 1 # Under relaxation factor for pressure
URFT = 1 # Under relaxation factor for temperature

## Air Properties at 25ºC 1atm
#mu  = fe.Constant(0.00001825) # Dynamic viscosity               [kg/ms]
mu  = fe.Constant(0.001825)    # Dynamic viscosity               [kg/ms]
rho = fe.Constant(1.225)       # Density                         [kg/m3]
k   = fe.Constant(0.02551)     # Thermal conductivity            [W/mK]
cp  = fe.Constant(20.7)        # Heat capacity constant pressure [J/kgK]
cv  = fe.Constant(718)         # Heat capacity constant volume   [J/kgK]
RR  = fe.Constant(287)         # Gas constant                    [J/kgK]
b   = fe.Constant((0, 0))      # Body accelerations              [m/s2]

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
Temp = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres, Temp])
FS   = fe.FunctionSpace(Mesh, M)

TF     = fe.TrialFunction(FS)
(w, q, s) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
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
    
WeakForm = Mnt + Cnt + Egy

J = fe.derivative(WeakForm, TFsol, TF)  

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DCircuit_facet_region.xml')

Entry      = 1
BottomWall = 2
Exit       = 3
TopWall    = 4
Circle1T = 5
Circle1B = 6
Circle2T = 7
Circle2B = 8
Circle3T = 9
Circle3B = 10
Circle4T = 11
Circle4B = 12
Circle5T = 13
Circle5B = 14
Circle6T = 15
Circle6B = 16
Circle7T = 17
Circle7B = 18
Circle8T = 19
Circle8B = 20

NoSlip    = fe.Constant((0, 0))
POut      = fe.Constant(0)
InletFlow = fe.Expression(['2*((4*x[1]*(0.1-x[1]))/(0.1*0.1))','0'], degree=2)
Temp1     = fe.Constant(400)
Temp2     = fe.Constant(273)

EntryBC      = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
BottomWallBC = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
ExitBC       = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
TopWallBC    = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)
Circle1TBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle1T)
Circle1BBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle1B)
Circle2TBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle2T)
Circle2BBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle2B)
Circle3TBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle3T)
Circle3BBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle3B)
Circle4TBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle4T)
Circle4BBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle4B)
Circle5TBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle5T)
Circle5BBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle5B)
Circle6TBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle6T)
Circle6BBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle6B)
Circle7TBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle7T)
Circle7BBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle7B)
Circle8TBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle8T)
Circle8BBC = fe.DirichletBC(FS.sub(0), NoSlip, DomainBoundaries, Circle8B)

EntryBCT      = fe.DirichletBC(FS.sub(2), Temp2, DomainBoundaries, Entry)
BottomWallBCT = fe.DirichletBC(FS.sub(2), Temp2, DomainBoundaries, BottomWall)
TopWallBCT    = fe.DirichletBC(FS.sub(2), Temp2, DomainBoundaries, TopWall)
Circle1TBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle1T)
Circle1BBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle1B)
Circle2TBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle2T)
Circle2BBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle2B)
Circle3TBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle3T)
Circle3BBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle3B)
Circle4TBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle4T)
Circle4BBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle4B)
Circle5TBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle5T)
Circle5BBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle5B)
Circle6TBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle6T)
Circle6BBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle6B)
Circle7TBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle7T)
Circle7BBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle7B)
Circle8TBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle8T)
Circle8BBCT   = fe.DirichletBC(FS.sub(2), Temp1, DomainBoundaries, Circle8B)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, Circle1TBC, Circle1BBC, Circle2TBC, \
       Circle2BBC, Circle3TBC, Circle3BBC, Circle4TBC, Circle4BBC, Circle5TBC, Circle5BBC, \
       Circle6TBC, Circle6BBC, Circle7TBC, Circle7BBC, Circle8TBC, Circle8BBC, \
       EntryBCT, Circle1TBCT, Circle1BBCT, Circle2TBCT, \
       Circle2BBCT, Circle3TBCT, Circle3BBCT, Circle4TBCT, Circle4BBCT, Circle5TBCT, Circle5BBCT, \
       Circle6TBCT, Circle6BBCT, Circle7TBCT, Circle7BBCT, Circle8TBCT, Circle8BBCT]

InitialVel  = fe.interpolate(fe.Constant((fe.Constant(1), fe.Constant(0.01))), FS.sub(0).collapse()) 
InitialPres = fe.interpolate(fe.Constant(10), FS.sub(1).collapse())
InitialTemp = fe.interpolate(fe.Constant(273),  FS.sub(1).collapse())

fe.assign(TFsol.sub(0), InitialVel)
fe.assign(TFsol.sub(1), InitialPres)
fe.assign(TFsol.sub(2), InitialTemp)

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['maximum_iterations']   = 1
Parameters['newton_solver']['error_on_nonconvergence'] = False

TFsol_prev = fe.Function(FS)
(u_prev, p_prev, T_prev) = TFsol_prev.split()
TFsol_relaxed = fe.Function(FS)
(u_relaxed, p_relaxed, T_relaxed) = TFsol_relaxed.split()

Total_iterations = 7
for i in range(Total_iterations):
      print(f"\nIteration {i+1}: \n")
      Solver.solve()
      
      (u_star, p_star, T_star) = TFsol.split()
      
      u_relaxed.vector()[:] = URFu*u_star.vector()[:] + (1-URFu)*u_prev.vector()[:]
      p_relaxed.vector()[:] = URFp*p_star.vector()[:] + (1-URFp)*p_prev.vector()[:]
      T_relaxed.vector()[:] = URFT*T_star.vector()[:] + (1-URFT)*T_prev.vector()[:]

      fe.assign(TFsol.sub(0), u_relaxed)
      fe.assign(TFsol.sub(1), p_relaxed)
      fe.assign(TFsol.sub(2), T_relaxed)

      u_prev = u_relaxed
      p_prev = p_relaxed
      T_prev = T_relaxed

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

(Velocity, Pressure, Temperature) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')
Temperature.rename('Temperature','Temperature')
Reynolds.rename('Reynolds','Reynolds')
ShearStress.rename('ShearStress','ShearStress')
Peclet.rename('Peclet','Peclet')
DomainBoundaries.rename('DomainBoundaries','DomainBoundaries')

fe.File('ResultsCircuitHT/Peclet.pvd')      << Peclet
fe.File('ResultsCircuitHT/Velocity.xml')    << Velocity
fe.File('ResultsCircuitHT/Pressure.xml')    << Pressure
fe.File('ResultsCircuitHT/Reynolds.xml')    << Reynolds
fe.File('ResultsCircuitHT/ShearStress.xml') << ShearStress
fe.File('ResultsCircuitHT/Velocity.pvd')    << Velocity
fe.File('ResultsCircuitHT/Pressure.pvd')    << Pressure
fe.File('ResultsCircuitHT/Temperature.pvd') << Temperature
fe.File('ResultsCircuitHT/Reynolds.pvd')    << Reynolds
fe.File('ResultsCircuitHT/ShearStress.pvd') << ShearStress
fe.File('ResultsCircuitHT/Boundaries.pvd')  << DomainBoundaries

