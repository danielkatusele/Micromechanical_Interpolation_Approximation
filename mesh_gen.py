import gmsh
import numpy as np

from math import sqrt


###### Plotting mesh ###################
def mesh_plotter(gmsh):
    import sys
    print(sys.argv)
    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()


def generate_mesh_dipole(xLgth, yLgth, xDip, yDip, nDip, d, a, b, angle, plot=False):

    gmsh.initialize()
    r_dip = (1/2)*sqrt(xDip**2 + yDip**2) # location of center of dipole
    
    gdim = 2  # Geometric dimension of the mesh


    # Define geometry for iron cylinder
    p1 = gmsh.model.occ.addPoint(xLgth/2, yLgth/2, 0)
    p2 = gmsh.model.occ.addPoint(-xLgth/2, yLgth/2, 0)
    p3 = gmsh.model.occ.addPoint(-xLgth/2, -yLgth/2, 0)
    p4 = gmsh.model.occ.addPoint(xLgth/2, -yLgth/2, 0)
    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(p2, p3)
    l3 = gmsh.model.occ.addLine(p3, p4)
    l4 = gmsh.model.occ.addLine(p4, p1)
    rec = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
    material = gmsh.model.occ.addPlaneSurface([rec])
    gmsh.model.occ.synchronize()
    

    # Define geometry for background

    #bd = gmsh.model.occ.addDisk(0, 0, 0, R, R)
    #background = gmsh.model.occ.fragment([(2, bd)],[(2, material)], removeTool=False)
    #print("background= ", background )
    #gmsh.model.occ.synchronize()


    # Define the 4 dipoles starting with the 1st quadrant
    a_o = [(np.pi/4 + i*2*np.pi/nDip) for i in range(nDip)]
    #a_o = [(i*2*np.pi/N) for i in range(N)]
    #chg_o = [(2, gmsh.model.occ.addDisk(r_chargeo*np.cos(v), r_chargeo*np.sin(v), 0, r, r)) for v in a_o]
    if nDip == 1:
        dipole = [(2, gmsh.model.occ.addDisk(0, 0, 0, a, b)) for v in a_o]
    else:
        dipole = [(2, gmsh.model.occ.addDisk(r_dip*np.cos(v), r_dip*np.sin(v), 0, a, b)) for v in a_o]

    for t in range(len(a_o)):
        if nDip == 1:
            gmsh.model.occ.rotate([(dipole[t][0], dipole[t][1])], 0, 0, 0, 0, 0, 1., np.tan(a_o[t])*angle)
            pass
        else:
            gmsh.model.occ.rotate([(dipole[t][0], dipole[t][1])], r_dip*np.cos(a_o[t]), r_dip*np.sin(a_o[t]), 0, 0, 0, 1., np.tan(a_o[t])*np.pi/4)
    


    # Define the inner charges
    #a_i = [(np.pi/4 + i*2*np.pi/N) for i in range(N)]
    #chg_i = [(2, gmsh.model.occ.addDisk(r_chargei*np.cos(v), r_chargei*np.sin(v), 0, r, r)) for v in a_i]
    gmsh.model.occ.synchronize()
    # Resolve all boundaries of the different wires in the background domain
    all_surfaces = []
    all_surfaces.extend(dipole)
    #all_surfaces.extend(chg_i)
    mat_domain = gmsh.model.occ.fragment([(2, material)], all_surfaces)
    #whole_domain = (background[0] + mat_domain[0], background[1][0] + mat_domain[1][0])
    whole_domain = mat_domain
    gmsh.model.occ.synchronize()

    # For full-quad mesh
    #for c in gmsh.model.getEntities(1):
    #    gmsh.model.mesh.setTransfiniteCurve(c[1], 33)
    #for s in gmsh.model.getEntities(2):
    #    gmsh.model.mesh.setRecombine(s[0], s[1])
        #gmsh.model.mesh.setSmoothing(s[0], s[1], 100)


    # Create physical markers for the different wires.
    # We use the following markers:
    # - Vacuum: 0
    # - Material: 1
    # - Outer charges: $[2,3,\dots,N+1]$
    # - Inner charges: $[N+2,\dots, 2\cdot N+1]
    dipole_tag = 1
    background_surfaces = []
    other_surfaces = []
    material_surfaces = []

    volume_f = 0
    area_d = np.pi*a*b

    for dom in whole_domain[0]:
        com = gmsh.model.occ.getCenterOfMass(dom[0], dom[1])
        mass = gmsh.model.occ.getMass(dom[0], dom[1])
        print("dom= ", dom, "mass= ", mass, "area_d= ", area_d, "com= ", com)
        
        
        
        # Identify material by its mass
        #if np.isclose(mass, (xLgth*yLgth - 2*N*np.pi*r**2)):
        #    gmsh.model.addPhysicalGroup(dom[0], [dom[1]], tag=0)
        #    material_surfaces.append(dom)  
        
         
            # Identify the outer charges by their center of mass
        if np.less_equal(np.abs(mass- area_d), 1e-6):
            gmsh.model.addPhysicalGroup(dom[0], [dom[1]], dipole_tag)
            dipole_tag +=1
            other_surfaces.append(dom)
        # Identify the inner charges by their center of mass
            
        # Identify the background material by its center of mass
        elif np.greater(mass, area_d):
        #    background_surfaces.append(dom[1]) 
            gmsh.model.addPhysicalGroup(dom[0], [dom[1]], tag=0)
            material_surfaces.append(dom)
            volume_f += mass
            
    print("area_d= ", area_d, "volume_f= ", volume_f, "total= ", area_d + volume_f)
    volume_fraction = area_d/(area_d  + volume_f)
    


            
    # Add marker for the vacuum
    #gmsh.model.addPhysicalGroup(2, background_surfaces, tag=0)



    # Create mesh resolution that is fine around the wires and
    # iron cylinder, coarser the further away you get
    
    gmsh.model.mesh.field.add("Distance", 1)
    edges = gmsh.model.getBoundary(other_surfaces, oriented=False)
    #edge2 = gmsh.model.getBoundary(material_surfaces)
    #edge3 = edges + edge2
    gmsh.model.mesh.field.setNumbers(1, "EdgesList", [e[1] for e in edges])
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "IField", 1)
    gmsh.model.mesh.field.setNumber(2, "LcMin", d/4)
    gmsh.model.mesh.field.setNumber(2, "LcMax", 1.5*d)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.9*a)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 3.5*a)
    gmsh.model.mesh.field.add("Min", 5)
    gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(5)
    # Generate mesh
    # For full-triangular mesh 
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    
    # For full-quad mesh
    #gmsh.option.setNumber('Mesh.RecombineAll', 1)
    #gmsh.option.setNumber('Mesh.RecombinationAlgorithm', 3)
    #gmsh.option.setNumber('Mesh.ElementOrder', 1)

    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.mesh.generate(gdim)   
    gmsh.model.mesh.optimize("Netgen")
    gmsh.model.mesh.removeDuplicateNodes() 
    gmsh.model.mesh.removeDuplicateElements()

    if plot:
        mesh_plotter(gmsh)

    from dolfinx.io import gmsh as gmshio
    from mpi4py import MPI
    model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    mesh_data = gmshio.model_to_mesh(gmsh.model,mesh_comm,  model_rank, gdim=2)
    domn = mesh_data.mesh
    ct = mesh_data.cell_tags
    ft = mesh_data.facet_tags
    ph_d = mesh_data.physical_groups
    gmsh.finalize()

    from math import degrees
    lbl = str(round(degrees(angle)))

    if a !=  b:
        mesh_id = "mesh/mesh_ellipse" + str(nDip) + "orient" + lbl + "size_" + str(xLgth) + ".xdmf"
    else:
        mesh_id = "mesh/mesh_circle" + str(nDip) + "size_" + str(xLgth) + ".xdmf"

    import dolfinx
    points = domn.geometry
    with dolfinx.io.XDMFFile(MPI.COMM_WORLD, mesh_id, "w") as xdmf:
        xdmf.write_mesh(domn)
        xdmf.write_meshtags(ct, domn.geometry)
    
    print("Mesh generated with volume fraction: ", volume_fraction*100, "%")
    return mesh_id, volume_fraction*100

def generate_mesh_rigidInc(xLgth, yLgth, xDip, yDip, nDip, d, a, b, angle=np.pi/4, plot=False):

    gmsh.initialize()
    r_dip = (1/2)*sqrt(xDip**2 + yDip**2) # location of center of dipole
    
    gdim = 2  # Geometric dimension of the mesh


    # Define geometry for iron cylinder
    p1 = gmsh.model.occ.addPoint(xLgth/2, yLgth/2, 0)
    p2 = gmsh.model.occ.addPoint(-xLgth/2, yLgth/2, 0)
    p3 = gmsh.model.occ.addPoint(-xLgth/2, -yLgth/2, 0)
    p4 = gmsh.model.occ.addPoint(xLgth/2, -yLgth/2, 0)
    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(p2, p3)
    l3 = gmsh.model.occ.addLine(p3, p4)
    l4 = gmsh.model.occ.addLine(p4, p1)
    rec = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
    material = gmsh.model.occ.addPlaneSurface([rec])
    gmsh.model.occ.synchronize()

    # Define geometry for background

    #bd = gmsh.model.occ.addDisk(0, 0, 0, R, R)
    #background = gmsh.model.occ.fragment([(2, bd)],[(2, material)], removeTool=False)
    #print("background= ", background )
    #gmsh.model.occ.synchronize()


    # Define the 4 dipoles starting with the 1st quadrant
    a_o = [(np.pi/4 + i*2*np.pi/nDip) for i in range(nDip)]
    #a_o = [(i*2*np.pi/N) for i in range(N)]
    #chg_o = [(2, gmsh.model.occ.addDisk(r_chargeo*np.cos(v), r_chargeo*np.sin(v), 0, r, r)) for v in a_o]
    if nDip == 1:
        dipole = [(2, gmsh.model.occ.addDisk(0, 0, 0, a, b)) for v in a_o]
    else:
        dipole = [(2, gmsh.model.occ.addDisk(r_dip*np.cos(v), r_dip*np.sin(v), 0, a, b)) for v in a_o]

    for t in range(len(a_o)):
        if nDip == 1:
            gmsh.model.occ.rotate([(dipole[t][0], dipole[t][1])], 0, 0, 0, 0, 0, 1., np.tan(a_o[t])*angle)
            pass
        else:
            gmsh.model.occ.rotate([(dipole[t][0], dipole[t][1])], r_dip*np.cos(a_o[t]), r_dip*np.sin(a_o[t]), 0, 0, 0, 1., np.tan(a_o[t])*np.pi/4)
    


    # Define the inner charges
    #a_i = [(np.pi/4 + i*2*np.pi/N) for i in range(N)]
    #chg_i = [(2, gmsh.model.occ.addDisk(r_chargei*np.cos(v), r_chargei*np.sin(v), 0, r, r)) for v in a_i]
    gmsh.model.occ.synchronize()
    # Resolve all boundaries of the different wires in the background domain
    all_surfaces = []
    all_surfaces.extend(dipole)
    #all_surfaces.extend(chg_i)
    #mat_domain = gmsh.model.occ.cut([(2, material)], all_surfaces)
    mat_domain = gmsh.model.occ.fragment([(2, material)], all_surfaces)
    whole_domain = mat_domain
    gmsh.model.occ.synchronize()

    # For full-quad mesh
    #for c in gmsh.model.getEntities(1):
    #    gmsh.model.mesh.setTransfiniteCurve(c[1], 33)
    #for s in gmsh.model.getEntities(2):
    #    gmsh.model.mesh.setRecombine(s[0], s[1])
        #gmsh.model.mesh.setSmoothing(s[0], s[1], 100)


    # Create physical markers for the different wires.

    material_surfaces = []; inclusion_surfaces = []
    inclusion_facet = []; domain_bdry =[]
    incl_tag = 1
    surf = gmsh.model.getEntities(dim=gdim)
    for dom in surf:
        com = gmsh.model.occ.getCenterOfMass(dom[0], dom[1])
        mass = gmsh.model.occ.getMass(dom[0], dom[1])
        area_d = 2*np.pi*a*b
        
        
        # Identify material by its mass
        #if np.isclose(mass, (xLgth*yLgth - 2*N*np.pi*r**2)):
        #    gmsh.model.addPhysicalGroup(dom[0], [dom[1]], tag=0)
        #    material_surfaces.append(dom)  
         
            # Identify the outer charges by their center of mass
        if np.less(mass, area_d):
            gmsh.model.addPhysicalGroup(dom[0], [dom[1]], incl_tag)
            incl_tag +=1
            inclusion_surfaces.append(dom)
        # Identify the inner charges by their center of mass
            
        # Identify the background material by its center of mass
        elif np.greater(mass, area_d):
        #    background_surfaces.append(dom[1])
            gmsh.model.addPhysicalGroup(dom[0], [dom[1]], tag=0)

    gmsh.model.occ.synchronize()
        
    bdry = gmsh.model.getEntities(dim=gdim - 1)
    for boundary in bdry:
        mass = gmsh.model.occ.getMass(boundary[0], boundary[1])
        cmass = gmsh.model.occ.getCenterOfMass(boundary[0], boundary[1])
        if not (np.allclose(mass, xLgth) or np.allclose(mass, yLgth)):
            inclusion_facet.append(boundary[1])
        else:
            domain_bdry.append(boundary[1])
            
    
    gmsh.model.addPhysicalGroup(gdim -1, inclusion_facet, tag=1)
    gmsh.model.addPhysicalGroup(gdim -1, domain_bdry, tag=0)
    gmsh.model.occ.synchronize()
                
    # Add marker for the vacuum
    #gmsh.model.addPhysicalGroup(2, background_surfaces, tag=0)



    # Create mesh resolution that is fine around the wires and
    # iron cylinder, coarser the further away you get
    
    gmsh.model.mesh.field.add("Distance", 1)
    edges = gmsh.model.getBoundary(inclusion_surfaces, oriented=False)
    #edge2 = gmsh.model.getBoundary(material_surfaces)
    #edge3 = edges + edge2
    gmsh.model.mesh.field.setNumbers(1, "EdgesList", [e[1] for e in edges])
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "IField", 1)
    gmsh.model.mesh.field.setNumber(2, "LcMin", d/3)
    gmsh.model.mesh.field.setNumber(2, "LcMax", 1.5*d)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.9*a)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 3*a)
    gmsh.model.mesh.field.add("Min", 5)
    gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(5)
    # Generate mesh
    # For full-triangular mesh 
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.occ.synchronize()
    # For full-quad mesh
    #gmsh.option.setNumber('Mesh.RecombineAll', 1)
    #gmsh.option.setNumber('Mesh.RecombinationAlgorithm', 3)
    #gmsh.option.setNumber('Mesh.ElementOrder', 1)

    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.mesh.generate(gdim)   
    gmsh.model.mesh.optimize("Netgen")
    gmsh.model.mesh.removeDuplicateNodes() 
    gmsh.model.mesh.removeDuplicateElements()

    if plot:
            mesh_plotter(gmsh)

    from dolfinx.io import gmsh as gmshio
    from mpi4py import MPI
    model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    domn, ct, ft = gmshio.model_to_mesh(gmsh.model,mesh_comm,  model_rank, gdim=2)

    gmsh.finalize()

    if a !=  b:
        mesh_id = "mesh/micromagnet_mesh_ellipse_.xdmf"
    else:
        mesh_id = "mesh/micromagnet_mesh_circle_.xdmf"

    import dolfinx
    points = domn.geometry
    with dolfinx.io.XDMFFile(MPI.COMM_WORLD, mesh_id, "w") as xdmf:
        xdmf.write_mesh(domn)
        xdmf.write_meshtags(ct, points)
        xdmf.write_meshtags(ft, points)
    
        
    return mesh_id


def generate_mesh_inclusion3D(xLgth, yLgth, zLgth,
                              d, a, b, c, angle=np.pi/4, plot=False):

    gmsh.initialize()

    gdim = 3  # Geometric dimension of the mesh
    ri = 0.001 # Initial radius


    # Define geometry for the matrix
    matrix = gmsh.model.occ.addBox(-xLgth/2, -yLgth/2, -zLgth/2, 
                               xLgth, yLgth, zLgth)
    gmsh.model.occ.synchronize()
    
    # Define the inclusion
    inclusion = gmsh.model.occ.addSphere(0, 0, 0, ri)

    gmsh.model.occ.dilate([(gdim, inclusion)], 0, 0, 0, a/ri, b/ri, c/ri)
    gmsh.model.occ.synchronize()
    gmsh.model.occ.rotate([(gdim, inclusion)], 0, 0, 0, 0, 0, 1., angle)

    gmsh.model.occ.synchronize()

  
    # Resolve all boundaries of the matrix and inclusion
    all_surfaces = []
    all_surfaces.extend([(gdim, inclusion)])
    
    mat_domain = gmsh.model.occ.fragment([(gdim, matrix)], all_surfaces)
    whole_domain = mat_domain
    gmsh.model.occ.synchronize()

    #Define physical group
    inclusion_tag = 1
    background_surfaces = []
    inclusion_surfaces = []
    matrix_surfaces = []


    for dom in whole_domain[0]:
        com = gmsh.model.occ.getCenterOfMass(dom[0], dom[1])
        mass = gmsh.model.occ.getMass(dom[0], dom[1])
        area_d = 4*np.pi*a*b*c/3
        
  
            # Identify the inclusions by their mass
        if np.less(mass, area_d):
            gmsh.model.addPhysicalGroup(dom[0], [dom[1]], inclusion_tag)
            inclusion_tag +=1
            inclusion_surfaces.append(dom)
 
        # Identify the matrix by its mass
        elif np.greater(mass, area_d):
        #    background_surfaces.append(dom[1])
            gmsh.model.addPhysicalGroup(dom[0], [dom[1]], tag=0)
            matrix_surfaces.append(dom)
            
    gmsh.model.occ.synchronize()
    
    gmsh.model.mesh.field.add("Distance", 1)
    edges = gmsh.model.getBoundary(inclusion_surfaces, oriented=False)
    #edge2 = gmsh.model.getBoundary(material_surfaces)
    #edge3 = edges + edge2
    gmsh.model.mesh.field.setNumbers(1, "FacesList", [e[1] for e in edges])
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "IField", 1)
    gmsh.model.mesh.field.setNumber(2, "LcMin", d/3)
    gmsh.model.mesh.field.setNumber(2, "LcMax", d)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.9*a)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 3*a)
    gmsh.model.mesh.field.add("Min", 3)
    gmsh.model.mesh.field.setNumbers(3, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(3)
    # Generate mesh
    # For full-triangular mesh 
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 4)
    
    # For full-quad mesh
    #gmsh.option.setNumber('Mesh.RecombineAll', 1)
    #gmsh.option.setNumber('Mesh.RecombinationAlgorithm', 3)
    #gmsh.option.setNumber('Mesh.ElementOrder', 1)

    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.mesh.generate(gdim)   
    gmsh.model.mesh.optimize("Netgen")
    gmsh.model.mesh.removeDuplicateNodes() 
    gmsh.model.mesh.removeDuplicateElements()

    if plot:
            mesh_plotter(gmsh)

    from dolfinx.io.gmshio import model_to_mesh
    from mpi4py import MPI
    model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    domn, ct, x = model_to_mesh(gmsh.model,mesh_comm,  model_rank, gdim=gdim)
    gmsh.finalize()

    mesh_id = "mesh/micromagnet3D_mesh.xdmf"

    import dolfinx
    points = domn.geometry
    with dolfinx.io.XDMFFile(MPI.COMM_WORLD, mesh_id, "w") as xdmf:
        xdmf.write_mesh(domn)
        xdmf.write_meshtags(ct, domn.geometry)
    
        
    return mesh_id


def generate_mesh_holes(xLgth, yLgth, Space, d, r1, r2, meshid, plot=False):

    gmsh.initialize()

    
    gdim = 2  # Geometric dimension of the mesh
    from dolfinx.io.gmshio import model_to_mesh
    from mpi4py import MPI
    model_rank = 0
    mesh_comm = MPI.COMM_WORLD


    # Define geometry for iron cylinder
    p1 = gmsh.model.occ.addPoint(xLgth/2, yLgth/2, 0)
    p2 = gmsh.model.occ.addPoint(-xLgth/2, yLgth/2, 0)
    p3 = gmsh.model.occ.addPoint(-xLgth/2, -yLgth/2, 0)
    p4 = gmsh.model.occ.addPoint(xLgth/2, -yLgth/2, 0)
    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(p2, p3)
    l3 = gmsh.model.occ.addLine(p3, p4)
    l4 = gmsh.model.occ.addLine(p4, p1)
    rec = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
    material = gmsh.model.occ.addPlaneSurface([rec])
    gmsh.model.occ.synchronize()

    
    # Define the 5 major holes starting with the 1st quadrant
    x_H = 0; y_H = 0
    M_hole=[]  
    M_hole.append((2, gmsh.model.occ.addDisk(0, 0, 0, r1, r1)))  
    n_y = 0 
    while y_H <= yLgth/2:
        n_x = 0
        x_H = ((n_y%2) + 2*n_x)*(r1 + Space/2)
        while x_H <= xLgth/2: 
            if (x_H == 0 and y_H == 0): 
                pass
            else:
                a_m = [(np.pi/4 + i*np.pi/2 + i*(x_H==0)*np.pi/2) for i in range(2*(x_H>0) + 2)]
                Hle = [(2, gmsh.model.occ.addDisk(x_H*np.sign(np.cos(v)),\
                        y_H*np.sign(np.sin(v)), 0, r1, r1)) for v in a_m]
#
                [M_hole.append(k) for k in Hle]
            n_x += 1
            x_H = ((n_y%2) + 2*n_x)*(r1 + Space/2) 
        n_y += 1
        y_H = n_y*(r1 + Space/2)
        
   
    # Define the 5 major holes starting with the 1st quadrant
    x_H = 0; y_H = 0
    m_hole=[]  
    n_y = 0 
    while y_H <= yLgth/2:
        n_x = 0
        x_H = (((n_y%2)==0) + 2*n_x)*(r1 + Space/2)
        while x_H <= xLgth/2: 
            if (x_H == 0 and y_H == 0): 
                pass
            else:
                a_m = [(np.pi/4 + i*np.pi/2 + i*(x_H==0)*np.pi/2) for i in range(2*(x_H>0) + 2)]
                Hle = [(2, gmsh.model.occ.addDisk(x_H*np.sign(np.cos(v)),\
                        y_H*np.sign(np.sin(v)), 0, r2, r2)) for v in a_m]
                
                for idx, ang in enumerate(a_m):
                    gmsh.model.occ.rotate([(Hle[idx][0], Hle[idx][1])], x_H*np.sign(np.cos(ang)),\
                        y_H*np.sign(np.sin(ang)), 0, 0, 0, 1., np.pi/2)
     
                [m_hole.append(k) for k in Hle]
            n_x += 1
            x_H = (((n_y%2)==0) + 2*n_x)*(r1 + Space/2) 
        n_y += 1
        y_H = n_y*(r1 + Space/2)
        
    

    gmsh.model.occ.synchronize()
    # Resolve all boundaries of the different wires in the background domain
    all_surfaces = []
    all_surfaces.extend(M_hole)
    all_surfaces.extend(m_hole)
    #all_surfaces.extend(chg_i)
    mat_domain = gmsh.model.occ.cut([(2, material)], all_surfaces)
     
    gmsh.model.occ.synchronize()

    # For full-quad mesh
    #for c in gmsh.model.getEntities(1):
    #    gmsh.model.mesh.setTransfiniteCurve(c[1], 33)
    #for s in gmsh.model.getEntities(2):
    #    gmsh.model.mesh.setRecombine(s[0], s[1])
        #gmsh.model.mesh.setSmoothing(s[0], s[1], 100)


    elastomer_marker = 1
    if mesh_comm.rank == model_rank:
        surf = gmsh.model.getEntities(dim=gdim)
        assert(len(surf) == 1)
        gmsh.model.addPhysicalGroup(surf[0][0], [surf[0][1]], elastomer_marker)
        gmsh.model.setPhysicalName(surf[0][0], elastomer_marker, "Elastomer")
        gmsh.model.occ.synchronize()

    inner_bd_marker = 5; x_bd_ml = 1; x_bd_mr = 2; y_bd_mt = 3; y_bd_mb = 4
    inner_bd = []; x_bd_l = []; x_bd_r = []; y_bd_t = []; y_bd_b = []
    if mesh_comm.rank == model_rank:
        boundaries = gmsh.model.getBoundary(surf, oriented=False)
        for boundary in boundaries:
            mass = gmsh.model.occ.getMass(boundary[0], boundary[1])
            cmass = gmsh.model.occ.getCenterOfMass(boundary[0], boundary[1])
            if np.allclose(mass, 2*np.pi*r1):
                inner_bd.append(boundary[1])
            elif np.allclose(mass, 2*np.pi*r2):
                inner_bd.append(boundary[1])
            elif np.allclose(cmass, [-xLgth/2, 0, 0]):
                x_bd_l.append(boundary[1]) 
            elif np.allclose(cmass, [xLgth/2, 0, 0]):
                x_bd_r.append(boundary[1])
            elif np.allclose(cmass, [0, -yLgth/2, 0]):
                y_bd_b.append(boundary[1]) 
            elif np.allclose(cmass, [0, yLgth/2, 0]):
                y_bd_t.append(boundary[1])

          
        gmsh.model.addPhysicalGroup(1, x_bd_l, x_bd_ml)
        gmsh.model.addPhysicalGroup(1, x_bd_r, x_bd_mr)
        
        gmsh.model.addPhysicalGroup(1, y_bd_t, y_bd_mt)
        gmsh.model.addPhysicalGroup(1, y_bd_b, y_bd_mb)
          
        gmsh.model.addPhysicalGroup(1, inner_bd, inner_bd_marker)
        gmsh.model.occ.synchronize()
        

        
    
    gmsh.option.setNumber('Mesh.MeshSizeMin', d)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 4 * d)
    gmsh.option.setNumber('General.Terminal', 0)

    gmsh.option.setNumber("Mesh.Algorithm", 7)
    


    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.mesh.generate(gdim)   
    gmsh.model.mesh.removeDuplicateNodes() 
    gmsh.model.mesh.removeDuplicateElements()
    
    if plot:
        mesh_plotter(gmsh)

    domn, _, ft = model_to_mesh(gmsh.model,mesh_comm,  model_rank, gdim=2)
    gmsh.finalize()

    
    import dolfinx
    with dolfinx.io.XDMFFile(mesh_comm, meshid, "w") as xdmf:
        xdmf.write_mesh(domn)
        xdmf.write_meshtags(ft, domn.geometry) 
        
        
def real_Element(domain, comm, Subspace, loc):
    import dolfinx
    
    # Find the master dof (x = 0.0) and owning process ##### Translation
    nodes = dolfinx.mesh.locate_entities_boundary(domain, 0, \
        lambda x: np.isclose(x[0], loc[0]) & np.isclose(x[1], loc[1]))
    master_dof = dolfinx.fem.locate_dofs_topological(Subspace, 0, nodes, remote=False)
    master_owner = np.full_like(master_dof, comm.rank)
    master_dof = np.concatenate(comm.allgather(master_dof), dtype=np.int64)
    master_owner = np.concatenate(comm.allgather(master_owner), dtype=np.int32)

    # Check for correct number of master dofs (should be 1)
    if master_dof.size != 1:
        raise RuntimeError(f"Expected to find 1 master dof; found {master_dof.size}")
    else:
        print(f"[{comm.rank}]: master_dof = {master_dof}")
        print(f"[{comm.rank}]: master_owner = {master_owner}")

    # Find the slave nodes
    nodes_ = dolfinx.mesh.locate_entities(domain, 0, \
        lambda x: ~(np.isclose(x[0], loc[0]) & np.isclose(x[1], loc[1])))
    slave_dof = dolfinx.fem.locate_dofs_topological(Subspace, 0, nodes_, remote=False)
    print(f"[{comm.rank}]: slave_dof = {slave_dof}")

    # Create arrays for use by mpc.add_constraint
    from petsc4py.PETSc import ScalarType
    masters = np.broadcast_to(master_dof, slave_dof.shape)
    owners = np.broadcast_to(master_owner, slave_dof.shape)
    coeffs = np.ones_like(masters, dtype=ScalarType)
    offsets = np.arange(masters.size + 1, dtype=np.int32)

    return slave_dof, masters, coeffs, owners, offsets

def Energy_Mechanics(u, mu, penalty):
    import ufl
    
    # Spatial dimension
    d = len(u)

    # Identity tensor
    I = ufl.variable(ufl.Identity(d))

    # Deformation gradient
    F = ufl.variable(I + ufl.grad(u))

    # Right Cauchy-Green tensor
    C = ufl.variable(F.T * F)

    # Invariants of deformation tensors
    Ic = ufl.variable(ufl.tr(C))
    IIc = ufl.variable(1/2*ufl.tr(C)*ufl.tr(C) + 1/2*ufl.tr(C*C))
    J  = ufl.variable(ufl.det(F))


    # Elasticity parameters
    # Stored strain energy density (compressible neo-Hookean model)
    W_mech = (mu / 2) * (Ic - 2.) 
    # Stress
    
    # Incompressible Neo-Hookean
    Engy = W_mech + (penalty/2)*(J-1.0)*(J-1.0) 
    
    return Engy, F

    

    


