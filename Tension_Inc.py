# FEA LAGRANGE M CIRCLE
##############################
##### CIRCLE #########
##############################
##############################


import numpy as np

from petsc4py import PETSc
from dolfinx import default_scalar_type, log
import dolfinx
from dolfinx.fem import (dirichletbc, Expression, Constant, Function, 
                         functionspace, locate_dofs_topological,
                         petsc, form, create_matrix, create_vector)
from dolfinx.mesh import locate_entities_boundary
from ufl import (TestFunction, TrialFunction, Measure, Identity, variable, dx, 
                 grad, inner, tr, det, diff, derivative, SpatialCoordinate, cross, atan2)
#from dolfinx.fem.petsc import NonlinearProblem
#from dolfinx.nls.petsc import NewtonSolver
from mpi4py import MPI
import matplotlib.pyplot as plt
from basix.ufl import element, mixed_element 

from time import *

from mesh_gen import generate_mesh_dipole

        
        
def point_eval(xpts, ypts, ufunc):
    from dolfinx import geometry
    bb_tree = geometry.bb_tree(domain, domain.topology.dim)
    cells = []
    points_on_proc = []
    # Find cells whose bounding-box collide with the the points
    points = np.zeros((3, len(xpts)))
    points[0] = xpts
    points[1] = ypts
    cell_candidates = geometry.compute_collisions_points(bb_tree, points.T)
    # Choose one of the cells that contains the point
    colliding_cells = geometry.compute_colliding_cells(domain, cell_candidates, points.T)
    for i, point in enumerate(points.T):
        if len(colliding_cells.links(i)) > 0:
            points_on_proc.append(point)
            cells.append(colliding_cells.links(i)[0])

    points_on_proc = np.array(points_on_proc, dtype=np.float64)
    u_values = ufunc.eval(points_on_proc, cells)  
    
    return u_values

def radius_ellipse(a, b, α, γ):
    return 1.0/ np.sqrt((np.cos(α)**2 *np.cos(γ)**2)/a**2 + (np.cos(γ)**2 *np.sin(α)**2)/b**2 \
        + (2* np.cos(α)* np.cos(γ)* np.sin(α)* np.sin(γ))/a**2 - (2 *np.cos(α) *np.cos(γ)* \
            np.sin(α)* np.sin(γ))/b**2 + (np.cos(α)**2 *np.sin(γ)**2)/b**2 + (np.sin(α)**2 *\
                np.sin(γ)**2)/a**2)
    
def x_ellipse(a, b, α, γ):
    return radius_ellipse(a, b, α, γ)*np.cos(γ)

def y_ellipse(a, b, α, γ):
    return radius_ellipse(a, b, α, γ)*np.sin(γ)

def θapproxtension(a, b, λv, α):
    return np.arctan2((a**4 - b**4)*(-(1/λv) + λv)* np.sin(2*α),\
        (a**4 + 6* a**2 *b**2 + b**4) *(1/λv + λv) + (a**4 -\
            b**4)* (-(1/λv) + λv)* np.cos(2*α))

def θapproxShear(a, b, λv, α):
    return np.arctan2(λv* (a**4 + 6* a**2 *b**2 + b**4 + (-a**4 + b**4)* np.cos(2*α)) ,
                      2* (a**4 + 6* a**2 *b**2 + b**4) + (a**4 - b**4)* λv* np.sin(2*α)) 


# Start time
st_time = time()

a = 0.2
b = 0.1
c = np.sqrt(a**2 - b**2) # Coordinate of focal point
d = 0.015 # Size of mesh
R = 3          # Radius of dom
N = 1          # Number of dipole
xL = 3.    # Length of elastomer in x-dir
yL = 3.    # Length of elastomer in y-dir
xC1 = 0.4      # Distance between charge magnets in x-dir
yC1 = 0.4      # Distance between elastomers in y-dir
angl = np.pi/4

M1 = 1.0
M2 = 0.


MaxStr = 2.
incB = 0.001


domain_size = [3.]
xl=3.
domain_size = [3.]
element_size = [0.01, 0.015]


for x_l in domain_size:
    z = 1.
    iteration = 0
    frame = 0.0
    G = 1.0

    mesh_id, vol = generate_mesh_dipole(xLgth=x_l, yLgth=x_l, xDip=xC1, 
                                yDip=yC1, nDip=N, d=d, a=a, b=b, angle=angl, plot=False)

    print("Volume fraction of inclusion:", vol)

    with dolfinx.io.XDMFFile(MPI.COMM_WORLD, mesh_id, "r") as xdmf:
        domain = xdmf.read_mesh(name="mesh")
        ct = xdmf.read_meshtags(domain, name="cell_tags")
        
    DG = dolfinx.fem.functionspace(domain, ("DG", 0))
    material_tags = np.unique(ct.values)
    mu = Function(DG) 
    mu.x.array[ct.find(0)] = 1.0    # Matrix
    mu.x.array[ct.find(1)] = 1.e5    # Inclusion (high stiffness)


    tdim = domain.topology.dim
    num_cells = domain.topology.index_map(tdim).size_global
    cells = np.arange(num_cells, dtype=np.int32)

    # Cell size for each local cell
    h = domain.h(tdim, cells)

    # Average element size on this MPI rank
    h_avg_local = np.mean(h)
    print("Average local h =", h_avg_local)

    ## Boundary condition
    # Tension
    class MacroF():
        def __init__(self, t):
            self.t = t

        def __call__(self, x):
            vals = np.zeros((2, x.shape[1]), dtype=PETSc.ScalarType)
            vals[0] = (self.t - 1.0)*x[0]
            vals[1] = ((1.0/self.t) - 1.0)*x[1]
            return vals

 

    V = functionspace(domain, ("Lagrange", 2, (domain.geometry.dim, )))
    Q = functionspace(domain, ("Lagrange", 1))

    tdim = domain.topology.dim
    facets = locate_entities_boundary(domain, tdim - 1, lambda x: np.full(x.shape[1], True))
    dofs = locate_dofs_topological(V, tdim - 1, facets)

    # Boundary Condition
    Disp_macro = MacroF(z)

    u_D = Function(V)
    u_D.interpolate(Disp_macro)
    bc = dirichletbc(u_D, dofs)

    uh, δu, du = Function(V), TestFunction(V), TrialFunction(V)
    ph, δp, dp = Function(Q), TestFunction(Q), TrialFunction(Q)

    # Spatial dimension
    d = len(uh)

    # Identity tensor
    I = variable(Identity(d))

    # Deformation gradient
    F = variable(I + grad(uh))

    # Right Cauchy-Green tensor
    C = variable(F.T * F)

    # Invariants of deformation tensors
    Ic = variable(tr(C))
    Jac = variable(det(F))

    Constraint = variable((1/2)*(C - I))

    # Elasticity parameters
    # Stored strain energy density (incompressible neo-Hookean model)
    metadata = {"quadrature_degree": 16}
    ds = Measure('ds', domain=domain, metadata=metadata)
    dx = Measure("dx", domain=domain, subdomain_data=ct, metadata=metadata)


    K = Constant(domain, 1.e5)
    γ = Constant(domain, 1.e6)
    mu2 = Constant(domain, 1.e0)
    psi_m = (mu2 / 2) * (Ic - 2)*dx + (K/2)*(Jac - 1)**2*dx
    psi_i = γ*inner(Constraint, Constraint)*dx(1)
    psi = psi_m + psi_i
    # Define form F (we want to find u such that F(u) = 0)
    Ψ = derivative(psi, uh, δu)

    δΨ = derivative(Ψ, uh, du)


    ##############  NEWTON SOLVER  #####################
    # Create PETSc forms from FEniCS forms
    # Create nonlinear problem

    petsc_options = {
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
            "ksp_error_if_not_converged": True,
            "snes_monitor": None,
            "snes_monitor": None,
            "snes_error_if_not_converged": True
        }

    problem = dolfinx.fem.petsc.NonlinearProblem(
        Ψ,                       # residual forms: [derivative(psi, uh, δu), derivative(psi, ph, δp)]
        [uh],                # solution functions
        J=δΨ,                    # Jacobian forms: nested 2x2 list
        bcs=[bc],
        petsc_options_prefix="hyperelasticity_lm_",
        petsc_options={
        "snes_type": "newtonls",
        "snes_linesearch_type": "basic",
        "snes_monitor": None,
        "snes_max_it": 100,
        "snes_atol": 1e-8,
        "snes_rtol": 1e-9,
        "snes_stol": 0.0,
        "snes_converged_reason": None,
        
        "ksp_type": "preonly",
        "ksp_converged_reason": None,
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",

        # MUMPS tuning (often helps for large-ish problems)
        "mat_mumps_icntl_14": 100,   # workspace
        },
    )

    snes = problem.solver


    #log.set_log_level(log.LogLevel.INFO)
    tval0 = -1.5

    from math import degrees
    lbl = str(round(degrees(angl)))

    sol_id = "Sol/Affine_micromagnet_sol_ellipse_tensionLM"+ lbl + ".bp"
    from dolfinx.io import VTXWriter
    writer = VTXWriter(domain.comm, sol_id, [uh], "BP4")

        
    rot = []
    rot2 = []
    Mrot = []
    stretch = []


    uhc = uh.x.petsc_vec.copy()
    uhc.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

    while z <= MaxStr:
        # Update Diriclet boundary condition
        Disp_macro.t = z
        u_D.interpolate(Disp_macro)
        print("================  lambda = ", z, "====================")
        try:
            # Solve the first iterations inaccurately
            dolfinx.fem.petsc.assign([uh], problem.x)
            problem.solve()
            assert snes.getConvergedReason() > 0
            # Pull the SNES solution vector back into the Function objects
            dolfinx.fem.petsc.assign(problem.x, [uh])
        except:
            print("Did not converge")
            #uh.x.destroy()
            snes.destroy()
            break
        
        if iteration % 10. == 0.0:
            #beta = np.linspace(0., 2*np.pi, 150)  , angl + np.pi
            beta = np.array([angl, np.pi + angl])
            xpts = x_ellipse(0.5*a, 0.5*b, angl, beta)
            ypts = y_ellipse(0.5*a, 0.5*b, angl, beta)
            P_org = np.array([xpts[0] - xpts[1], ypts[0] - ypts[1]])
            #P_org = np.array([xpts[0], ypts[0]])
            

            uh.x.scatter_forward()
            Disp = point_eval(xpts, ypts, uh)
            Disp2 = point_eval(np.array([0.]), np.array([0.]), uh)
            new_D = np.array_split(Disp, 2)
            
            Disp3 = (Disp - Disp2)[0]
            
            P_def = np.array([P_org[0] + Disp3[0], P_org[1] + Disp3[1]])
            
            cost = np.dot(P_org, P_def) /(np.linalg.norm(P_org)*np.linalg.norm(P_def))
            x = SpatialCoordinate(domain)
            
            cross_pr = (uh[1] + x[1])*x[0] - (uh[0] + x[0])*x[1]
            dot_pr = (uh[0] + x[0])*x[0] + (uh[1] + x[1])*x[1]
            angle_num = atan2(cross_pr, dot_pr)
            
            #xdmf_s.write_function(uh, z)
            
            rot.append(np.arccos(cost))
            
            Mrot.append(np.arctan2((M2 *np.cos(angl - np.arccos(cost)) + M1 *np.sin(angl - np.arccos(cost))) , 
                                    (M1* np.cos(angl - np.arccos(cost)) - M2 *np.sin(angl - np.arccos(cost))))-angl)
            #print("Changeinlength",np.linalg.norm(P_org)-np.linalg.norm(P_def))
            #print("rot= ", np.arccos(cost), "analytical= ", θapproxShear(a, b, z, angl) )
            stretch.append(z)
            
            I = Identity(len(uh))
            # Deformation gradient
            Fp = (I + grad(uh))
            
            Cp = (Fp.T * Fp)
            
            
            WEng = functionspace(domain, ("CG", 1))
            Eg = Function(WEng)
            en = (mu / 2) * (tr(Cp) - 2) 
            Eg_expr = Expression(en, WEng.element.interpolation_points)
            Eg.interpolate(Eg_expr)
            Eg.name = "Energy"
            #xdmf_s.write_function(Eg, z)
            
            WCp = functionspace(domain, ("CG", 1))
            FCp = Function(WCp)
            Cp_expr = Expression(angle_num, WCp.element.interpolation_points)
            FCp.interpolate(Cp_expr)
            FCp.name = "Angle_numerical"
            #xdmf_s.write_function(FCp, z)
            
        iteration += 1
        z += incB
        
    dataN = {
                "Tensile load": stretch,
                "Numerical": Mrot
            }

    import pandas as pd
    vol_f = round(vol*100, 1)
    df = pd.DataFrame(dataN)
    df.to_excel("Tension_Pen_vol_f_" + str(vol) +"_45_another.xlsx", index=False)  # Uses openpyxl automatically


# Create DataFrame

    
def θpolardecomp(a, b, λv):
    return np.arctan2(-(np.sqrt(2 + λv**2 - λv* np.sqrt(4 + λv**2)) - np.sqrt(2 \
        + λv* (λv + np.sqrt(4 + λv**2))))/np.sqrt(2), ((λv + np.sqrt(4 + λv**2))*\
            np.sqrt(2 + λv**2 - λv* np.sqrt(4 + λv**2)) + (-λv + np.sqrt(4 + λv**2))*\
                np.sqrt(2 + λv* (λv + np.sqrt(4 + λv**2))))/(2* np.sqrt(2) *np.sqrt(4 + λv**2)))

λ = np.linspace(1., z, 100)

θ = []
θapp = []
θPD = []
Diffθapp = []
for τ in λ:
    θapp.append(np.arctan2(M2 *np.cos(angl - θapproxtension(a, b, τ, angl)) + M1 *np.sin(angl - θapproxtension(a, b, τ, angl)) , 
                           M1* np.cos(angl - θapproxtension(a, b, τ, angl)) - M2 *np.sin(angl - θapproxtension(a, b, τ, angl))) - angl)
    
