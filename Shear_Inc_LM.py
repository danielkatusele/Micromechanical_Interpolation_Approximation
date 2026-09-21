# FEA LAGRANGE M CIRCLE
##############################
##### CIRCLE #########
##############################
##############################


import numpy as np

from petsc4py import PETSc
from dolfinx import default_scalar_type, log
import dolfinx
import ufl
from dolfinx.fem import (dirichletbc, Expression, Constant, Function, 
                         functionspace, locate_dofs_topological,
                         petsc, form, assemble_scalar)
from dolfinx.mesh import locate_entities_boundary
from ufl import (TestFunction, TrialFunction, Measure, Identity, variable, dx, 
                 grad, inner, tr, det, diff, derivative, MixedFunctionSpace)
#from dolfinx.fem.petsc import NonlinearProblem
#from dolfinx.nls.petsc import NewtonSolver
from mpi4py import MPI
import matplotlib.pyplot as plt
from basix.ufl import element, mixed_element 

from time import *

from mesh_gen import generate_mesh_dipole

# Class for interfacing with the SNES
import ufl
import typing
import dolfinx
import multiphenicsx.fem
import multiphenicsx.fem.petsc


#plt.rcParams['text.usetex'] = False
#plt.rcParams['font.family'] = 'Times'
plt.rcParams["font.size"] = 12.0
plt.rcParams["xtick.top"] = True 
plt.rcParams["ytick.right"] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams["legend.fancybox"] = False
plt.rcParams["legend.framealpha"] = 1
plt.margins(0.0)

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
d = 0.05 # Size of mesh
R = 3          # Radius of dom
N = 1          # Number of dipole
xL = 5.    # Length of elastomer in x-dir
yL = 5.    # Length of elastomer in y-dir
xC1 = 0.4      # Distance between charge magnets in x-dir
yC1 = 0.4      # Distance between elastomers in y-dir
#angl = [0., np.pi/6, np.pi/3., np.pi/4, np.pi/6, np.pi]
#colors = ["purple", "green", "red", "blue", "orange", "brown"]

M1 = 1.0
M2 = 0.

#angl = [np.pi/6]
angl = np.pi/3
colors = ["blue"]
label = [ "$\pi / 4$", "$\pi / 3$", "$\pi / 2$"]

MaxStr = 0.02
incB = 0.01

 
domain_size = [3.]
xl=3.
domain_size = [3.]
element_size = [0.015, 0.04, 0.1]

for el_s in element_size:
    z = 0.
    iteration = 0
    frame = 0.0
    G = 1.0
    
    mesh_id, vol = generate_mesh_dipole(xLgth=xL, yLgth=xL, xDip=xC1, 
                                yDip=yC1, nDip=N, d=el_s, a=a, b=b, angle=angl, plot=False)

    print("Volume fraction of inclusion:", vol)

    with dolfinx.io.XDMFFile(MPI.COMM_WORLD, mesh_id, "r") as xdmf:
        domain = xdmf.read_mesh(name="mesh")
        ct = xdmf.read_meshtags(domain, name="cell_tags")

    tdim = domain.topology.dim

    num_cells = domain.topology.index_map(tdim).size_local
    cells = np.arange(num_cells, dtype=np.int32)

    # Cell size for each local cell
    h = domain.h(tdim, cells)

    # Average element size on this MPI rank
    h_avg_local = np.mean(h)
    print("Average local h =", h_avg_local)

    num_cells_local = domain.topology.index_map(tdim).size_local
    num_cells_global = domain.topology.index_map(tdim).size_global

    if domain.comm.rank == 0:
        print(f"Number of elements: {num_cells_global}")

    
        
    DG = dolfinx.fem.functionspace(domain, ("DG", 0))
    mu = Function(DG)
    mu.x.array[:] = 1.0          # matrix
    mu.x.array[ct.find(1)] = 1e6 # inclusion
    mu.x.scatter_forward()

    # Boundary condition
    #class MacroF():
    #    def __init__(self, t):
    #        self.t = t
    #
    #    def __call__(self, x):
    #        vals = np.zeros((2, x.shape[1]), dtype=PETSc.ScalarType)
    #        vals[0] = (self.t - 1.0)*x[0]
    #        vals[1] = ((1.0/self.t) - 1.0)*x[1]
    #        return vals

    class MacroF():
        def __init__(self, t):
            self.t = t
#
        def __call__(self, x):
            vals = np.zeros((2, x.shape[1]), dtype=PETSc.ScalarType)
            vals[0] = self.t*x[1]
            vals[1] = 0.
            return vals


    V = functionspace(domain, ("Lagrange", 2, (domain.geometry.dim, )))
    Q = functionspace(domain, ("Lagrange", 1))

    tdim = domain.topology.dim
    facets = locate_entities_boundary(domain, tdim - 1, lambda x: np.full(x.shape[1], True))
    dofs = locate_dofs_topological(V, tdim - 1, facets)

    print("Number of dofs on boundary:", len(dofs))


    Disp_macro = MacroF(z)
    u_D = Function(V)
    u_D.interpolate(Disp_macro)
    bc = dirichletbc(u_D, dofs)

    uh, δu, du = Function(V), TestFunction(V), TrialFunction(V)
    ph, δp, dp = Function(Q), TestFunction(Q), TrialFunction(Q)


    # Identity tensor
    I = variable(Identity(len(uh)))

    # Deformation gradient
    F = variable(I + grad(uh))

    # Right Cauchy-Green tensor
    C = variable(F.T * F)

    # Invariants of deformation tensors
    Ic = variable(tr(C))
    Jac = variable(det(F))

    # Elasticity parameters
    # Stored strain energy density (incompressible neo-Hookean model)
    metadata = {"quadrature_degree": 8}
    ds = Measure('ds', domain=domain, metadata=metadata)
    dx = Measure("dx", domain=domain, subdomain_data=ct, metadata=metadata)


    K = Constant(domain, 1.e5)
    psi = (mu / 2) * (Ic - 2)*dx + ph*(Jac - 1)*dx

    # Define form F (we want to find u such that F(u) = 0)
    Ψ = [derivative(psi, uh, δu), derivative(psi, ph, δp)]

    δΨ = [[derivative(Ψ[0], uh, du), derivative(Ψ[0], ph, dp)],
        [derivative(Ψ[1], uh, du), derivative(Ψ[1], ph, dp)]]
    
    import traceback
            
    ##############  NEWTON SOLVER  #####################
    # Create PETSc forms from FEniCS forms
    # Create nonlinear problem

    # Create problem
    # Solve
    petsc_options = {
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "ksp_error_if_not_converged": True,
        "snes_monitor": None,
        "snes_error_if_not_converged": True
    }

    problem = dolfinx.fem.petsc.NonlinearProblem(
        Ψ,                       # residual forms: [derivative(psi, uh, δu), derivative(psi, ph, δp)]
        [uh, ph],                # solution functions
        J=δΨ,                    # Jacobian forms: nested 2x2 list
        bcs=[bc],
        petsc_options_prefix="hyperelasticity_lm_",
        petsc_options={
            "snes_rtol": 1.0e-7,
            "snes_atol": 1.0e-8,
            "snes_max_it": 250,
            "snes_linesearch_type": "bt",
            "snes_monitor": None,
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
        },
    )

    snes = problem.solver


    max_iterations = 25
    normed_diff = 0
    tol = 1e-5
    m = 15
    u_prev = dolfinx.fem.Function(V)
    diff = dolfinx.fem.Function(V)



    #log.set_log_level(log.LogLevel.INFO)
    tval0 = -1.5

    from math import degrees
    
    orient = str(round(degrees(angl)))
        
    rot = []
    Mrot = []
    rot2 = []
    stretch = []
    strore_energy = []

    while z <= MaxStr:
        # Update Diriclet boundary condition
        Disp_macro.t = z
        u_D.interpolate(Disp_macro)
        print("================  lambda = ", z, "====================")
        try:

            dolfinx.fem.petsc.assign([uh, ph], problem.x)
            problem.solve()
            assert snes.getConvergedReason() > 0
            # Pull the SNES solution vector back into the Function objects
            dolfinx.fem.petsc.assign(problem.x, [uh, ph])
        
            #writer.write(z)
        except Exception as e:
            print("Exception during solve:", repr(e))
            print("Did not converge")
            #solution.destroy()
            #snes.destroy()
            break
        
        if iteration % 1 == 0.0:
            #beta = np.linspace(0., 2*np.pi, 150)
            beta = np.array([angl, angl + np.pi])
            xpts = x_ellipse(0.5*a, 0.5*b, angl, beta)
            ypts = y_ellipse(0.5*a, 0.5*b, angl, beta)
            P_org = np.array([xpts[0] - xpts[1], ypts[0] - ypts[1]])
            Disp = point_eval(xpts, ypts, uh)
            new_D = np.array_split(Disp, 2)
            
            Disp2 = (new_D[0] - new_D[1])[0]
            P_def = np.array([(xpts[0] - xpts[1]) + Disp2[0], (ypts[0] - ypts[1]) + Disp2[1]])
            print("Changeinlength",np.linalg.norm(P_org)-np.linalg.norm(P_def))
            cost = np.dot(P_org, P_def) /(np.linalg.norm(P_org)*np.linalg.norm(P_def))
            rot.append(np.arccos(cost))
            stretch.append(z)
            
            Mrot.append(np.arctan2((M2 *np.cos(angl - np.arccos(cost)) + M1 *np.sin(angl - np.arccos(cost))) , 
                                   (M1* np.cos(angl - np.arccos(cost)) - M2 *np.sin(angl - np.arccos(cost))))-angl)
            
            
            # Example data
            
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
            
            # Compute the total energy
            L2_error = form(Eg * ufl.dx)
            error_local = assemble_scalar(L2_error)
            
            WCp = functionspace(domain, ("CG", 1, (2,2)))
            FCp = Function(WCp)
            Cp_expr = Expression(Cp.T * Cp, WCp.element.interpolation_points)
            FCp.interpolate(Cp_expr)
            FCp.name = "Rotation"
            #xdmf_s.write_function(FCp, z)
            
        iteration += 1
        z += incB

    plt.plot(stretch, rot,"--", label="Numerical, elsize= {vol_f:.1f}%".format(vol_f=vol), linewidth=3., zorder=3)
    dataN = {
            "Shear load": stretch,
            "Numerical": rot}

        # Create DataFrame
    import pandas as pd
    df = pd.DataFrame(dataN)
    
    vol_f = round(vol*100, 1)
    df.to_excel("Shear_LM_n_el" + str(num_cells_global) + ".xlsx", index=False)  # Uses openpyxl automatically
    

    #xdmf_s.close()
       
       
    
    
def θapproxShear(a, b, λv, α):
    return np.arctan2(λv* (a**4 + 6* a**2 *b**2 + b**4 + (-a**4 + b**4)* np.cos(2*α)) ,
                    2* (a**4 + 6* a**2 *b**2 + b**4) + (a**4 - b**4)* λv* np.sin(2*α)) 


    
        
def θpolardecomp(a, b, λv):
    return np.arctan2(-(np.sqrt(2 + λv**2 - λv* np.sqrt(4 + λv**2)) - np.sqrt(2 \
        + λv* (λv + np.sqrt(4 + λv**2))))/np.sqrt(2), ((λv + np.sqrt(4 + λv**2))*\
            np.sqrt(2 + λv**2 - λv* np.sqrt(4 + λv**2)) + (-λv + np.sqrt(4 + λv**2))*\
                np.sqrt(2 + λv* (λv + np.sqrt(4 + λv**2))))/(2* np.sqrt(2) *np.sqrt(4 + λv**2)))

λ = np.linspace(0., z, 100)


θ = []
θapp = []
θPD = []
Diffθapp = []
for τ in λ:
    orient = θapproxShear(a, b, τ, angl)*180/np.pi
    θapp.append(orient)
    θPD.append(θpolardecomp(a, b, τ))
    

    
    