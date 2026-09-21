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
