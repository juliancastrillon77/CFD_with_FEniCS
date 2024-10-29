# Julian Castrillon
# 3D Steady Incompressible

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

Mesh = fe.Mesh('Mesh/3D Mesh/3DMesh.xml')

vi  = 0.2                        # Inlet velocity [m/s]
mu  = fe.Constant(1)             # Dynamic viscosity  [kg/ms]
rho = fe.Constant(10)            # Density            [kg/m3]
b   = fe.Constant((0, -9.81, 0)) # Body accelerations [m/s2]

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres])
FS   = fe.FunctionSpace(Mesh, M)

TF     = fe.TrialFunction(FS)
(w, q) = fe.TestFunctions(FS)

TFsol  = fe.Function(FS)
(u, p) = fe.split(TFsol)

WeakForm =   rho*fe.inner(w,fe.grad(u)*u)*fe.dx       \
           - fe.div(w)*p*fe.dx                        \
           - rho*fe.dot(w,b)*fe.dx                    \
           + mu*fe.inner(fe.grad(w),fe.grad(u))*fe.dx \
           + q*fe.div(u)*fe.dx

J = fe.derivative(WeakForm, TFsol, TF)  

DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/3D Mesh/3DMesh_facet_region.xml')

LeftWall            = 46
RightWall           = 47
Entry               = 48
Exit                = 49
BottomWall          = 50
TopWall             = 51
SphereForwardRight  = 52
SphereForwardLeft   = 53
SphereBackwardLeft  = 54
SphereBackwardRight = 55
ConeForwardRight    = 56
ConeForwardLeft     = 57
ConeBackwardLeft    = 58
ConeBackwardRight   = 59

NoSlip    = fe.Constant((0, 0, 0))
POut      = fe.Constant(0)
InletFlow = fe.Constant((vi, 0, 0))

LeftWallBC            = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, LeftWall)
RightWallBC           = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, RightWall)
EntryBC               = fe.DirichletBC(FS.sub(0), InletFlow, DomainBoundaries, Entry)
ExitBC                = fe.DirichletBC(FS.sub(1), POut,      DomainBoundaries, Exit)
BottomWallBC          = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, BottomWall)
TopWallBC             = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, TopWall)
SphereForwardRightBC  = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, SphereForwardRight)
SphereForwardLeftBC   = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, SphereForwardLeft)
SphereBackwardLeftBC  = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, SphereBackwardLeft)
SphereBackwardRightBC = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, SphereBackwardRight)
ConeForwardRightBC    = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, ConeForwardRight)
ConeForwardLeftBC     = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, ConeForwardLeft)
ConeBackwardLeftBC    = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, ConeBackwardLeft)
ConeBackwardRightBC   = fe.DirichletBC(FS.sub(0), NoSlip,    DomainBoundaries, ConeBackwardRight)

BCs = [LeftWallBC, RightWallBC, EntryBC, ExitBC, BottomWallBC, TopWallBC, SphereForwardRightBC, \
       SphereForwardLeftBC, SphereBackwardLeftBC, SphereBackwardRightBC, ConeForwardRightBC,    \
       ConeForwardLeftBC, ConeBackwardLeftBC, ConeBackwardRightBC]

Problem = fe.NonlinearVariationalProblem(WeakForm, TFsol, BCs, J)
Solver  = fe.NonlinearVariationalSolver(Problem)

Parameters = Solver.parameters
Parameters['newton_solver']['absolute_tolerance']   = 1e-8
Parameters['newton_solver']['relative_tolerance']   = 1e-7
Parameters['newton_solver']['maximum_iterations']   = 7
Parameters['newton_solver']['relaxation_parameter'] = 1.0

Solver.solve()

(Velocity, Pressure) = TFsol.split(deepcopy = True)
Velocity.rename('Velocity','Velocity')
Pressure.rename('Pressure','Pressure')

FileVel  = fe.File('Results/Velocity.pvd')
FilePres = fe.File('Results/Pressure.pvd')
FileVel  << Velocity
FilePres << Pressure