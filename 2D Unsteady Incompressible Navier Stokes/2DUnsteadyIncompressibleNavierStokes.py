# Julian Castrillon
# 2D Unsteady Incompressible Navier Stokes Equations

import os
import numpy as np
import fenics as fe

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

fe.set_log_level(fe.LogLevel.INFO)
fe.parameters['form_compiler']['representation']     = 'uflacs'
fe.parameters['form_compiler']['optimize']           = True
fe.parameters['form_compiler']['cpp_optimize']       = True
fe.parameters["form_compiler"]["cpp_optimize_flags"] = '-O2 -funroll-loops'

Mesh = fe.Mesh('Mesh/2D Mesh/2DMesh.xml')

dt    = 0.01
t_end = 11

# Fluid Properties
mu = fe.Constant(1/500)  # Dynamic viscosity  [kg/ms]
rho = fe.Constant(1)     # Density            [kg/m3]
b  = fe.Constant((0, 0)) # Body accelerations [m/s2]

theta = fe.Constant(0.5)      
n     = fe.FacetNormal(Mesh)   
h     = fe.CellDiameter(Mesh) 

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres])
FS   = fe.FunctionSpace(Mesh, M)

TF     = fe.TrialFunction(FS)
(w, q) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
(u, p) = fe.split(TFsol)

WeakFormB =   rho*fe.inner(w,fe.grad(u)*u)*fe.dx       \
            - fe.div(w)*p*fe.dx                        \
            - rho*fe.dot(w,b)*fe.dx                    \
            + mu*fe.inner(fe.grad(w),fe.grad(u))*fe.dx \
            - q*fe.div(u)*fe.dx                        \

TFsol0   = fe.Function(FS)
(u0, p0) = fe.split(TFsol0)

WeakFormF =   rho*fe.inner(w,fe.grad(u0)*u0)*fe.dx      \
            - fe.div(w)*p*fe.dx                         \
            - rho*fe.dot(w,b)*fe.dx                     \
            + mu*fe.inner(fe.grad(w),fe.grad(u0))*fe.dx \
            - q*fe.div(u0)*fe.dx                        \
                           
WeakForm = fe.inner((u-u0)/dt,w)*fe.dx + (theta)*WeakFormB + (1-theta)*WeakFormF

### Streamwise Upwind Petrov Galerkin (SUPG) stabilization for convection #######################
#vnorm = fe.sqrt(fe.dot(u0,u0))
#tau   = ((1/(theta*dt))**2+(2*vnorm/h)**2+9*(4*mu/h**2)**2)**(-0.5)
#R =  idt*(u-u0)+theta*(fe.grad(u)*u+fe.grad(p)-mu*fe.div(fe.grad(u))-b) \
#    +(1-theta)*(fe.grad(u0)*u0+fe.grad(p)-mu*fe.div(fe.grad(u0))-b)
#SUPG = tau*fe.inner(fe.grad(w)*u0,R)*fe.dx(metadata={'quadrature_degree':4})
#WeakForm += SUPG

### Petrov Galerkin Pressure Stabilzation (PSPG) stabilization for pressure field // The Ladyzhenskaya-Babuska-Brezzi condition not met
#PSPG = -tau*fe.inner(fe.grad(q),R)*fe.dx(metadata={'quadrature_degree':4})
#WeakForm += PSPG

J = fe.derivative(WeakForm, TFsol, TF)

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DMesh_facet_region.xml')

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

EntryBC         = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
BottomWallBC    = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
ExitBC          = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
TopWallBC       = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)
CircleBC        = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, Circle)
TriangleLeftBC  = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleLeft)
TriangleRightBC = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TriangleRight)

BCs = [EntryBC, BottomWallBC, ExitBC, TopWallBC, CircleBC, TriangleLeftBC, TriangleRightBC]

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['maximum_iterations']   = 4
Parameters['newton_solver']['relaxation_parameter'] = 1.0

FileVel  = fe.File('Results/Velocity.pvd')
FilePres = fe.File('Results/Pressure.pvd')

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