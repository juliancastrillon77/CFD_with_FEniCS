"""A demo FEniCS script for solvin an unsteady version of Reynolds Averaged
Navier Stokes equations using Q2P1 elements and a one equation Spalart Allmaras
model for turbulence closure. There are no stabilization terms implemented.
Additionally, for a more complete version of the model, a Boussinesq form of
the equations, including energy equation has to be included in the full system.

This particular demo pertains to jet flow in a 2D channel geometry where an
inlet jet is specified in a region at the center of the inlet.

Last Revised:
-------------
Spring 2022

Disclaimer:
-----------
Developed for computational fluid dynamics class taught at the Paul M Rady
Department of Mechanical Engineering at the University of Colorado Boulder
by Chayut Teeraratkul and Prof. Debanjan Mukherjee.

All inquiries addressed to Prof. Mukherjee at debanjan@Colorado.Edu
"""

import fenics as fe
import sys

#----------------------------------------------
# Define compiler specific optimization options
#----------------------------------------------
fe.set_log_level(fe.LogLevel.INFO)
fe.parameters['form_compiler']['representation']    = 'uflacs'
fe.parameters['form_compiler']['optimize']          = True
fe.parameters['form_compiler']['cpp_optimize']      = True
fe.parameters["form_compiler"]["cpp_optimize_flags"]= '-O2 -funroll-loops'

#------------------------------------------------------------
# Define a user defined expression to implement a jet inlet
# boundary condition that can be assigned to a region on the
# boundary between a y=a and a y=b
#------------------------------------------------------------
class INLET_JET(fe.UserExpression):
    '''
    This has value = (U,0) when a < x[1] < b and 0 everywhere else
    '''
    def __init__(self, U, a, b, **kwargs):
        self.U = U
        self.a = a
        self.b = b
        super().__init__(**kwargs)
    def eval(self, value, x):
        if self.a < x[1] and x[1] < self.b:
            value[0] = self.U
        return value
    def value_shape(self):
        return (2, )

#-------------------------------------------------------------------
# Define a user defined expression to implement a boundary condition
# that specifies a non-zero eddy viscosity at the jet inlet
# segment of the boundary
#-------------------------------------------------------------------
class EDV_INLET(fe.UserExpression):
    def __init__(self, nu0, a, b, **kwargs):
        self.a = a
        self.b = b
        self.nu0 = nu0
        super().__init__(**kwargs)
    def eval(self, value, x):
        if self.a < x[1] and x[1] < self.b:
            value[0] = self.nu0
        return value
    def value_shape(self):
        return ()

#--------------------------------------------------
# Physical parameters and inputs for the simulation
#--------------------------------------------------
nu      = 1.0
rho0    = 1.0
dt      = 1.0e-3
jet_Y0  = 0.8
jet_Y1  = 1.2
jet_U   = 1000.00
theta_value = 0.5

c_b1_value  = 0.1355
sig_value   = 2.0/3.0
c_b2_value  = 0.622
kt_value    = 0.41
c_w2_value  = 0.30
c_w3_value  = 2.0
c_v1_value  = 7.1
c_t3_value  = 1.2
c_t4_value  = 0.5

ev_0    = 3.0*nu
isNonZeroInletEddyVisc = True

isWriteVelComponentWise = True

outfile_u   = "output-2/SA_u.pvd"
outfile_v   = "output-2/SA_v.pvd"
outfile_p   = "output-2/SA_p.pvd"
outfile_nh  = "output-2/SA_nh.pvd"
outfile_nt  = "output-2/SA_nt.pvd"

zero    = fe.Constant(0.0)
zero_v  = fe.Constant((0.0,0.0))

#-------------------------------------------------------------------
# Define a rectangular mesh using the FEniCS internal meshing method
# Mesh dimensions and spacing are specified internally.
# This avoids using external Gmsh meshes, but can be used only for
# simple domain geometries.
# DX and DY are the mesh sizing intervals in each direction.
# NEX and NEY are the number of elements in each direction.
#-------------------------------------------------------------------
WIDTH   = 10
HEIGHT  = 2
DX      = 0.075
DY      = DX
NEX     = int(WIDTH/DX)
NEY     = int(HEIGHT/DY)
mesh    = fe.RectangleMesh(fe.Point(0,0), fe.Point(WIDTH, HEIGHT), NEX, NEY)

#------------------------------------------
# Define the domain boundaries for the mesh
#------------------------------------------
def left(x, on_boundary):
    return on_boundary and x[0] < fe.DOLFIN_EPS
def right(x, on_boundary):
    return on_boundary and WIDTH - x[0] < fe.DOLFIN_EPS
def bottom(x, on_boundary):
    return on_boundary and x[1] < fe.DOLFIN_EPS
def top(x, on_boundary):
    return on_boundary and HEIGHT - x[1] < fe.DOLFIN_EPS

#-------------------------------------------------------------------
# Define function spaces
# For turbulence model, we will need to define spaces for velocity
# and pressure as well as the eddy viscosity (for the one equation
# Spallart Almaras model version)
#-------------------------------------------------------------------
V   = fe.VectorElement("Lagrange", mesh.ufl_cell(), 2)
P   = fe.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
N   = fe.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
M   = fe.MixedElement([V, P, N])
W   = fe.FunctionSpace(mesh, M)

W_u     = W.sub(0)
W_p     = W.sub(1)
W_nh    = W.sub(2)
(w, q, y)  = fe.TestFunctions(W)

W0 = fe.Function(W)
Wn = fe.Function(W)

#------------------------------------------------------------------------
# Defining small numerical entries to be used for tolerance and precision
# and avoid divide by zero error due to the eddy viscosity model
#------------------------------------------------------------------------
eps     = fe.interpolate(fe.Constant(fe.DOLFIN_EPS), W.sub(1).collapse())
eps_v   = fe.interpolate(fe.Constant( (fe.DOLFIN_EPS, fe.DOLFIN_EPS) ), W.sub(0).collapse())

#-------------------------------------------------------------------
# Assigning initial conditions.
# For this implementation, setting non-zero initial conditions helps
# with avoiding divide by zero error in the turbulence modeling
#-------------------------------------------------------------------
fe.assign(W0.sub(0), eps_v)
fe.assign(W0.sub(1), eps)
fe.assign(W0.sub(2), eps)
fe.assign(Wn.sub(0), eps_v)
fe.assign(Wn.sub(1), eps)
fe.assign(Wn.sub(2), eps)

#---------------------------
# Define boundary conditions
#---------------------------

#------------------------------------------------
# Specify velocity boundary conditions first
#------------------------------------------------
inlet_jet       = INLET_JET(jet_U, jet_Y0, jet_Y1)
bc_vel_outlet   = fe.DirichletBC(W_p, zero, right)
bc_vel_inlet    = fe.DirichletBC(W_u, inlet_jet, left)
bc_vel_wall_t   = fe.DirichletBC(W_u, zero_v, top)
bc_vel_wall_b   = fe.DirichletBC(W_u, zero_v, bottom)

#------------------------------------------------
# Specify eddy viscosity boundary conditions next
#------------------------------------------------
bc_ev_wall_t = fe.DirichletBC( W_nh, zero, top )
bc_ev_wall_b = fe.DirichletBC( W_nh, zero, bottom )
bc_ev_outlet = fe.DirichletBC( W_nh, fe.Constant(3.0*nu), right )

if isNonZeroInletEddyVisc == True:
    evInlet     = EDV_INLET(ev_0,jet_Y0,jet_Y1)
    bc_ev_inlet = fe.DirichletBC( W_nh, evInlet, left )

#---------------------------------------------------------------------
# Collect all boundary conditions together for integrating into solver
#---------------------------------------------------------------------
bc_set  = []
bc_set.append( bc_vel_outlet )
bc_set.append( bc_vel_inlet )
bc_set.append( bc_vel_wall_t )
bc_set.append( bc_vel_wall_b )
bc_set.append( bc_ev_wall_t )
bc_set.append( bc_ev_wall_b )
bc_set.append( bc_ev_outlet )

if isNonZeroInletEddyVisc == True:
    bc_set.append( bc_ev_inlet )

#----------------------------------------------------------
# Turbulence Modeling: Resolution of wall distance for the
# turbulence model. Here it is a simple check of heightwise
# distance between top and bottom wall y and 1-y etc.
#----------------------------------------------------------
class dist_closest_wall(fe.UserExpression):
    '''
    This is just distance from top or bottom wall
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    def eval(self, value, x):
        if x[1] > HEIGHT/2:
            value[0] = HEIGHT - x[1]
        else:
            value[0] = x[1]
        return value
    def value_shape(self):
        return ()

d = dist_closest_wall()
d = fe.interpolate(d, W_nh.collapse())

#---------------------------------------------------------
# Turbulence Modeling: Definition of all the SA parameters
# that are not dependent on the solution variables
#---------------------------------------------------------
c_b1 = fe.Constant(c_b1_value)
sig  = fe.Constant(sig_value)
c_b2 = fe.Constant(c_b2_value)
kt   = fe.Constant(kt_value)
c_w2 = fe.Constant(c_w2_value)
c_w3 = fe.Constant(c_w3_value)
c_v1 = fe.Constant(c_v1_value)
c_t3 = fe.Constant(c_t3_value)
c_t4 = fe.Constant(c_t4_value)
c_w1 = c_b1/(kt**2) + (1.0 + c_b2)/sig

#----------------------------------------------------------
# Turbulence Modeling: Definition of all the SA parameters
# that are based on the solution variables
#----------------------------------------------------------
MIN = lambda a,b: (a+b - abs(a-b)) / fe.Constant(2)

def SA_params(sol):
    '''
    Return a dictionary of SA parameters dependent on the solution variables
    '''
    (u,p,nh) = fe.split(sol)
    Sij     = 0.5*( fe.grad(u) - fe.grad(u).T )
    Omega   = fe.sqrt( 2 * fe.inner(Sij, Sij) )
    chi     = nh/nu
    f_v1    = chi**3 / ( chi**3 + c_v1**3 )
    f_v2    = 1 - (chi / (1+chi*f_v1))
    f_t2    = c_t3 * fe.exp( -c_t4*chi**2 )
    Sh_v    = Omega + nh/(kt**2 + d**2) * f_v2
    Sh      = fe.conditional( fe.le(Sh_v, 0), 0.3*Omega, Sh_v)
    r       = fe.conditional( fe.le(Sh,0), 10, MIN(10, nh/(Sh*kt**2*d**2) ) )
    g       = r + c_w2*( r**6 - r )
    f_w     = g* ( (1 + c_w3**6)/(g**6 + c_w3**6)  ) ** (1./6.)

    param = {}
    param['f_v1'] = f_v1
    param['f_t2'] = f_t2
    param['S_hat'] = Sh
    param['f_w'] = f_w
    return param

#--------------------------------------------------------
# Define the combined weak form for a monolithic solution
#--------------------------------------------------------

#-----------------------------------
# Defining quadrature specifications
#-----------------------------------
quad_deg    = 4
dx          = fe.dx(metadata={"quadrature_degree": quad_deg})

#------------------------------------------------------------
# Theta-Galerkin: Parameters at previous and current timestep
#------------------------------------------------------------
param_0 = SA_params(W0)
param_n = SA_params(Wn)

(u0,p0,nh0) = fe.split(W0)
(un,pn,nhn) = fe.split(Wn)

nu_t0  = param_0['f_v1']*nh0
nu_tn  = param_n['f_v1']*nhn

def symmetric_strain_rate(u):
    return 0.5*(fe.grad(u) + fe.grad(u).T)

def stress(u,p,nu_t):
    strain = symmetric_strain_rate(u)
    I = fe.Identity(mesh.topology().dim())
    return -(p/rho0)*I + 2*(nu+nu_t)*strain

#-----------------------------------------------------------
# Theta-Galerkin: Contribution from continuity/mass balance
#-----------------------------------------------------------
MASS_eq = ( q*fe.div(un) )*dx

#--------------------------------------------------------
# Theta-Galerkin: Contribution from RANS momentum balance
#--------------------------------------------------------
theta = fe.Constant(theta_value)
T0 = stress(u0,pn,nu_t0)
Tn = stress(un,pn,nu_tn)


F0 = fe.inner( T0, fe.grad(w) ) + fe.inner( fe.grad(u0)*u0, w )



Fn = fe.inner( Tn, fe.grad(w) ) + fe.inner( fe.grad(un)*un, w )
MMT_eq = ( fe.inner( w, (un-u0)/dt ) + (1-theta)*F0 + theta*Fn )*dx

#-----------------------------------------------------------
# Theta-Galerkin: Contributions from Eddy viscosity equation
#-----------------------------------------------------------
def ev_rhs(u,nh,param):
    rhs  = y*fe.dot(u, fe.grad(nh))
    rhs -= y*c_b1*( 1 - param['f_t2'] )*param['S_hat']*nh
    rhs += y*(c_w1*param['f_w'] - c_b1/kt**2)*(nh/d)**2
    rhs += 1/sig * fe.dot( fe.grad(y), (nu+nh)*fe.grad(nh) )
    rhs -= c_b2/sig * y * fe.dot(fe.grad(nh), fe.grad(nh))
    return rhs

Y0 = ev_rhs(u0,nh0,param_0)
Yn = ev_rhs(un,nhn,param_n)

EVS_eq  = ( y*(nhn-nh0)/dt + (1.0 - theta)*Y0 + theta*Yn )*dx

#--------------------------------------------------------------------------
# Setting up the complete monolithic non-linear solver using Newton-Raphson
#--------------------------------------------------------------------------
weakform    = MASS_eq + MMT_eq + EVS_eq
J           = fe.derivative(weakform, Wn, fe.TrialFunction(W))
problem     = fe.NonlinearVariationalProblem(weakform, Wn, bc_set, J=J)
solver      = fe.NonlinearVariationalSolver(problem)

prm = solver.parameters
prm['newton_solver']['absolute_tolerance'] = 1E-8
prm['newton_solver']['relative_tolerance'] = 1E-7
prm['newton_solver']['maximum_iterations'] = 25
prm['newton_solver']['relaxation_parameter'] = 1.0

#--------------------------------
# The time-stepping solution loop
#--------------------------------

t = 0
ufid = fe.File(outfile_u)
if isWriteVelComponentWise == True:
    vfid = fe.File(outfile_v)
pfid = fe.File(outfile_p)
nhfid = fe.File(outfile_nh)
ntfid = fe.File(outfile_nt)

while t < 0.1:
    t+=dt
    solver.solve()
    (u_sol,p_sol,nh_sol) = Wn.split(deepcopy=True)

    fv1_sol = ((nh_sol/nu)**3)/(c_v1**3 + (nh_sol/nu)**3)
    nu_t    = nh_sol*fv1_sol
    nu_t    = fe.project(nu_t, fe.FunctionSpace(mesh,'Lagrange',1))

    if isWriteVelComponentWise == True:

        u1_sol = fe.project(fe.dot(u_sol,fe.Constant((1.0,0.0))), fe.FunctionSpace(mesh,'Lagrange',2))
        u2_sol = fe.project(fe.dot(u_sol,fe.Constant((0.0,1.0))), fe.FunctionSpace(mesh,'Lagrange',2))

        u1_sol.rename('uvel','velocity')
        u2_sol.rename('vvel','velocity')
        p_sol.rename('pres','pressure')
        nh_sol.rename('nuh','eddyvisc')
        nu_t.rename('nut','eddyvisc')

        ufid << u1_sol
        vfid << u2_sol
        pfid << p_sol
        nhfid << nh_sol
        ntfid << nu_t

    else:

        u_sol.rename('vel','velocity')
        p_sol.rename('pres','pressure')
        nh_sol.rename('nuh','eddyvisc')
        nu_t.rename('nut','eddyvisc')

        ufid << u_sol
        pfid << p_sol
        nhfid << nh_sol
        ntfid << nu_t

    W0.assign(Wn)
