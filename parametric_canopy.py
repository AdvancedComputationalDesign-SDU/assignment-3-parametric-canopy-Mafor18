"""
Assignment 3: Parametric Structural Canopy

Author: Martin Forná

Description:
This script generates a parametric structural canopy using depth maps and recursive geometry generation.
It explores different functions to control both the depth map and the fractal geometry to generate a
structural system composed of:
- A shell/gridshell
- A set of vertical supports

The script also combines different strategies for surface tessellation, quadratic, triangular and hexagonal elements to achieve 
a non-uniform tessellation of the input surface.

Note: This script is intended to be used within Grasshopper's Python scripting component.
"""

# Import necessary libraries
import Rhino
import Rhino.Geometry as rg
import rhinoscriptsyntax as rs
import math
import random
from scipy.spatial import Voronoi
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib import colormaps
import ghpythonlib.treehelpers as th

# Recursion_params has not been added as a grasshopper parameter yet as it seems to only transfer as a list instead of a dictionary
recursion_params = {
 'max_depth': 3,
 'vertical_angle': 20,
 'length': 4,
 'length_reduction': 0.5,
 'branches': 3,
 'angle_variation': 15
}

def ensure_surface(geometry):
    """
    Ensures the input geometry is a Rhino.Geometry.Surface. as passing it directly seems to reference the geometry 
    as a GUID instead of a surface

    Parameters:
    - geometry: The input geometry, potentially a Brep, GUID or other type.

    Returns:
    - A Rhino.Geometry.Surface object, if convertible.
    """
    # Coerce geometry to a Rhino object
    geometry = rs.coercegeometry(geometry)

    # Check if the geometry needs altering, skip if yes.
    if isinstance(geometry, rg.Surface):
        return geometry

    # If it's a Brep, try to extract the underlying surface and only continue if the Berp is a single face
    if isinstance(geometry, rg.Brep):
        if geometry.Faces.Count == 1: 
            return geometry.Faces[0].UnderlyingSurface()
        else:
            raise ValueError("The input Brep has multiple faces. Expected a single surface.")

    # Raise an error if not a valid surface
    raise ValueError("Input must be a Surface or a single-face Brep.")

# Attempt to transform given surface into a usable surface
try:
    base_surface = ensure_surface(base_surface)
    print("Successfully converted to a Surface!")
except ValueError as e:
    print(f"Error: {e}")


# Makes sure the depth divisions are passed as integers
depth_divisions = int(depth_divisions)

def generate_depth_map(surface, control_value, variation_type=""):
    """
    Modifies the input surface based on a control function to create a depth map.

    Parameters:
    - surface: The base surface (rg.Surface)
    - control_value: A numerical value controlling the depth variation (float)
    - variation_type: The type of depth variation ("Sine", "Ripple" or "Roof") (string)
    - depth_divisions: Number of divisions along U and V directions (int)

    Returns:
    - modified_surface: The surface after applying the depth map
    """
    
    # Set the domain the u (x) and v (y) directions
    u_domain = surface.Domain(0)
    v_domain = surface.Domain(1)  
    
    # Create an empty list to store modified points
    modified_points = []

    # Create a point grid with refinement depending on depth_division    
    for i in range(depth_divisions):
        u_param = u_domain[0] + (u_domain[1] - u_domain[0]) * i / (depth_divisions - 1)
        for j in range(depth_divisions):
            v_param = v_domain[0] + (v_domain[1] - v_domain[0]) * j / (depth_divisions - 1)
            
            # Evaluate the point on the surface at (u_param, v_param)
            point = surface.PointAt(u_param, v_param)
            
            # Depending on selected variation_type modify the surface in different ways 
            # Add control_value to change the intensity of the effect
            if variation_type == "Sine":
                # Sine and cosine effect
                depth_variation = (
                    math.sin(u_param * math.pi * 2) * math.cos(v_param * math.pi * 2) * control_value
                )

            elif variation_type == "Ripple":
                # Ripple effect: sin(10 * (x^2 + y^2)) / 10
                r_squared = point.X ** 2 + point.Y ** 2
                depth_variation = math.sin(10 * (r_squared) * control_value / 10)

            elif variation_type == "Roof":
                # Roof effect: create a peak in the center with linear slopes
                u_center = (u_domain[0] + u_domain[1]) / 2
                v_center = (v_domain[0] + v_domain[1]) / 2
                u_distance = abs(u_param - u_center) / (u_domain[1] - u_domain[0])
                v_distance = abs(v_param - v_center) / (v_domain[1] - v_domain[0])
                depth_variation = (1 - max(u_distance, v_distance)) * control_value

            else:
                raise ValueError("Invalid variation type. Choose 'Sine' or 'Ripple' or 'Roof'.")
            
            # Modify the Z-coordinate of the point
            modified_point = rg.Point3d(point.X, point.Y, point.Z + depth_variation)
            modified_points.append(modified_point)
    
    # Now creates a nurbs surface from the modified points
    # Pass the modified points as a flat list and provide grid dimensions
    modified_surface = rg.NurbsSurface.CreateThroughPoints(
        modified_points, depth_divisions, depth_divisions, 3, 3, False, False
    )

    return modified_surface

def generate_surface_points(surface):
    """
    Generates points along a given surface.

    Parameters:
    - surface: The base surface (rg.Surface)

    Returns:
    - surface_points: points along the input surface as rg.Point3d
    """
    
    # Set the domain the u (x) and v (y) directions
    u_domain = surface.Domain(0)
    v_domain = surface.Domain(1)  
    
    # Create an empty list to store modified points
    modified_points = []

    # Create a point grid with refinement depending on depth_division    
    for i in range(depth_divisions):
        u_param = u_domain[0] + (u_domain[1] - u_domain[0]) * i / (depth_divisions - 1)
        for j in range(depth_divisions):
            v_param = v_domain[0] + (v_domain[1] - v_domain[0]) * j / (depth_divisions - 1)
            
            # Evaluate the point on the surface at (u_param, v_param)
            point = surface.PointAt(u_param, v_param)
            
            # Modify the Z-coordinate of the point
            modified_point = rg.Point3d(point.X, point.Y, point.Z)
            modified_points.append(modified_point)
          
    # Convert the list to a data tree
    surface_points = th.list_to_tree(modified_points)

    return surface_points

# Divisions for the colormaps. Increase number to increase image resolution.
u_divisions = v_divisions = 200

def generate_depth_map_image(surface, u_divisions, v_divisions, colormap=''):
    """
    Generates a 2D depth map with a specified colormap based on the Z-values of the surface.

    Parameters:
    - surface: The input surface (rg.Surface)
    - u_divisions: Number of divisions in the U direction for sampling (int)
    - v_divisions: Number of divisions in the V direction for sampling (int)
    - colormap: Name of the colormap to apply (string)

    Returns:
    - colored_depth_map_image: A Pillow Image object with the depth map in the specified colormap.
    """
    # Get the domain of the surface.
    u_domain = surface.Domain(0)
    v_domain = surface.Domain(1)

    # Empty list to store z_values.
    z_values = []

    # Sample the surface at regular intervals in the U and V directions.
    for i in range(u_divisions):
        u_param = u_domain[0] + (u_domain[1] - u_domain[0]) * i / (u_divisions - 1)
        row_z_values = []
        for j in range(v_divisions):
            v_param = v_domain[0] + (v_domain[1] - v_domain[0]) * j / (v_divisions - 1)
            point = surface.PointAt(u_param, v_param)
            row_z_values.append(point.Z)
        z_values.append(row_z_values)

    # Convert the values into a numpy array.
    z_values = np.array(z_values)

    # Normalize the Z-values in range of 0-1.
    z_min = np.min(z_values)
    z_max = np.max(z_values)
    normalized_z = (z_values - z_min) / (z_max - z_min)

    # Apply the colormap
    colormap_function = colormaps[colormap]
    colored_data = colormap_function(normalized_z)
    
    # Convert the colormap to desired type
    colored_image = (colored_data[:, :, :3] * 255).astype(np.uint8)
    colored_depth_map_image = Image.fromarray(colored_image)

    return colored_depth_map_image

# Ensure mesh_rinfement is passed as an integer instead of float.
# Increase to mesh refinement parameter will create a finer mesh.
mesh_refinement = int(mesh_refinement)

def tessellate_surface(surface, strategy=''):
    """
    Tessellates the input surface using the specified strategy.

    Parameters:
    - surface: The surface to tessellate (rg.Surface)
    - strategy: The tessellation method ('quad', 'triangular', 'voronoi')

    Returns:
    - tessellated_mesh: A mesh representing the tessellated surface
    """
    # Sample the surface domain
    u_domain = surface.Domain(0)
    v_domain = surface.Domain(1)

    u_divisions = v_divisions = mesh_refinement 
    

    # Quadratic tessellation
    if strategy == 'quad':
        mesh = rg.Mesh()
        for i in range(u_divisions - 1):
            for j in range(v_divisions - 1):
                # Get corner points of each quad
                u0 = u_domain[0] + (u_domain[1] - u_domain[0]) * i / (u_divisions - 1)
                v0 = v_domain[0] + (v_domain[1] - v_domain[0]) * j / (v_divisions - 1)

                u1 = u_domain[0] + (u_domain[1] - u_domain[0]) * (i + 1) / (u_divisions - 1)
                v1 = v_domain[0] + (v_domain[1] - v_domain[0]) * (j + 1) / (v_divisions - 1)

                pt0 = surface.PointAt(u0, v0)
                pt1 = surface.PointAt(u1, v0)
                pt2 = surface.PointAt(u1, v1)
                pt3 = surface.PointAt(u0, v1)

                # Add a vertex for each point and join them together for create one quad cellmesh.
                v0_idx = mesh.Vertices.Add(pt0)
                v1_idx = mesh.Vertices.Add(pt1)
                v2_idx = mesh.Vertices.Add(pt2)
                v3_idx = mesh.Vertices.Add(pt3)

                mesh.Faces.AddFace(v0_idx, v1_idx, v2_idx, v3_idx)
                
        
        mesh.Normals.ComputeNormals()
        mesh.Compact()
        return mesh

    # Triangular tessellation
    elif strategy == 'triangular':
        mesh = rg.Mesh()
        for i in range(u_divisions):
            for j in range(v_divisions):
                # Generate points in a grid
                u = u_domain[0] + (u_domain[1] - u_domain[0]) * i / (u_divisions - 1)
                v = v_domain[0] + (v_domain[1] - v_domain[0]) * j / (v_divisions - 1)
                pt = surface.PointAt(u, v)
                mesh.Vertices.Add(pt)

        # Connect points to form triangles
        for i in range(u_divisions - 1):
            for j in range(v_divisions - 1):
                v0 = i * v_divisions + j
                v1 = (i + 1) * v_divisions + j
                v2 = v0 + 1
                v3 = v1 + 1

                mesh.Faces.AddFace(v0, v1, v2)
                mesh.Faces.AddFace(v2, v1, v3)

        mesh.Normals.ComputeNormals()
        mesh.Compact()
        return mesh

    # Hexagonal tessellation
    elif strategy == 'hexagonal':
        mesh = rg.Mesh()
        
        # Compute hexagon size based on mesh refinement
        u_domain = surface.Domain(0)
        v_domain = surface.Domain(1)
        u_size = (u_domain[1] - u_domain[0]) / mesh_refinement
        v_size = (v_domain[1] - v_domain[0]) / mesh_refinement
        hex_radius = min(u_size, v_size) / 2

        # Generate hexagons
        hex_height = math.sqrt(3) * hex_radius
        row_count = int((v_domain[1] - v_domain[0]) / hex_height) + 1
        col_count = int((u_domain[1] - u_domain[0]) / (1.5 * hex_radius)) + 1

        for row in range(row_count):
            for col in range(col_count):
                # Compute center of the hexagon
                center_u = u_domain[0] + col * 1.5 * hex_radius
                center_v = v_domain[0] + row * hex_height
                if col % 2 != 0: 
                    center_v += hex_height / 2
                
                if center_u > u_domain[1] or center_v > v_domain[1]:
                    continue

                center = surface.PointAt(center_u, center_v)
                
                # Generate hexagon vertices
                hex_vertices = []
                for k in range(6):
                    angle = math.pi / 3 * k
                    vertex_u = center_u + hex_radius * math.cos(angle)
                    vertex_v = center_v + hex_radius * math.sin(angle)
                    vertex = surface.PointAt(vertex_u, vertex_v)
                    hex_vertices.append(rg.Point3d(vertex.X, vertex.Y, vertex.Z))
                
                # Add the hexagon as triangles
                center_idx = mesh.Vertices.Add(center)
                for k in range(6):
                    v0_idx = mesh.Vertices.Add(hex_vertices[k])
                    v1_idx = mesh.Vertices.Add(hex_vertices[(k + 1) % 6])
                    mesh.Faces.AddFace(center_idx, v0_idx, v1_idx)
        
        mesh.Normals.ComputeNormals()
        mesh.Compact()
        return mesh


    else:
        raise ValueError(f"Unknown tessellation strategy: {strategy}")

def extract_mesh_edges(mesh):
    """
    Extracts edges from the tessellated mesh and converts them into Rhino lines.

    Parameters:
    - mesh: A Rhino.Geometry.Mesh object.

    Returns:
    - edge_lines: A list of Rhino.Geometry.LineCurve objects representing the edges.
    """
    edge_lines = []

    # Loop through each face of the mesh
    for face in mesh.Faces:
        # Get the corner vertices of the face then close the loop at the end
        if face.IsQuad:
            # For quad faces
            v0, v1, v2, v3 = face.A, face.B, face.C, face.D
            vertices = [v0, v1, v2, v3, v0]  
        elif face.IsTriangle:
            # For triangular faces
            v0, v1, v2 = face.A, face.B, face.C
            vertices = [v0, v1, v2, v0]  

        # Convert vertices to points and create lines between them
        for i in range(len(vertices) - 1):
            pt_start = mesh.Vertices[vertices[i]]
            pt_end = mesh.Vertices[vertices[i + 1]]
            line = rg.LineCurve(rg.Point3d(pt_start), rg.Point3d(pt_end))
            edge_lines.append(line)

    return edge_lines


def generate_recursive_supports(start_point, params, depth=0):
    """
    Generates recursive geometry (e.g., fractal patterns) for vertical supports.

    Parameters:
    - start_point: The starting point for recursion (rg.Point3d)
    - params: A dictionary containing parameters for recursion control
    - depth: The current recursion depth

    Returns:
    - curves: A list of generated Rhino.Geometry.LineCurve objects representing the supports
    """
    # Stops recursion if depth increases higher than max_depth to prevent infinite recursion
    if depth >= params['max_depth']:
        return []
    
    # Calculate the end point for the main branch
    length = params['length'] * (params['length_reduction'] ** depth)
    end_point = rg.Point3d(
        start_point.X,
        start_point.Y,
        start_point.Z + length
    )
    
    # Create the main vertical line
    line = rg.LineCurve(start_point, end_point)
    curves = [line]

    # Create branches
    for i in range(params['branches']):
        # Evenly distribute branches in XY-plane
        base_angle = 360 / params['branches'] * i  
        angle_variation = random.uniform(-params['angle_variation'], params['angle_variation'])
        branch_angle_xy = math.radians(base_angle + angle_variation)

        # Adjust angle between Z-axis and XY-plane
        z_to_xy_radians = math.radians(params['vertical_angle'])
        branch_vector = rg.Vector3d(
            math.cos(branch_angle_xy) * math.sin(z_to_xy_radians),  # X
            math.sin(branch_angle_xy) * math.sin(z_to_xy_radians),  # Y
            math.cos(z_to_xy_radians)                               # Z
        )
        branch_vector.Unitize()
        branch_vector *= length

        # Calculate branch end point
        branch_end_point = rg.Point3d(
            end_point.X + branch_vector.X,
            end_point.Y + branch_vector.Y,
            end_point.Z + branch_vector.Z
        )
        
        # Create the branch line
        branch_line = rg.LineCurve(end_point, branch_end_point)
        curves.append(branch_line)

        # Recursively generate sub-branches
        curves.extend(generate_recursive_supports(branch_end_point, params, depth + 1))

    return curves


if base_surface and depth_map_control and tessellation_strategy and recursion_params and support_points:
    # Generate modified surface with depth map
    modified_surface = generate_depth_map(base_surface, depth_map_control,variation_type)
    surface_points = generate_surface_points(modified_surface)
    
    # Generate and show the color map for the surface
    colored_depth_map_image = generate_depth_map_image(modified_surface, u_divisions, v_divisions, colormap = 'cool')
    # colored_depth_map_image.show()
    # # Save the image if desired
    # colored_depth_map_image.save("C:/Users/Martin/Documents/GitHub/assignment-3-parametric-canopy-Mafor18/images/roof_canopy_depth_map.png")

    # Tessellate the modified surface
    canopy_mesh = tessellate_surface(modified_surface, tessellation_strategy)
    
    # Extract and visualize edges of the tessellated mesh
    canopy_edges = extract_mesh_edges(canopy_mesh)

    # # Generate vertical supports
    supports = []
    for pt in support_points:
        curves = generate_recursive_supports(pt, recursion_params)
        supports.extend(curves)

    
else:
    # Handle cases where inputs are not provided
    print("Skipped something because of missing input")
    pass