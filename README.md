[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/t-2I_Jw_)
# Assignment 3: Parametric Structural Canopy

![Example Canopy](images/canopy.jpg)

## Objective

In this assignment, you will explore the creation of complex architectural structures using parametric design principles, building upon your knowledge of recursion, loops, and geometric transformations. Your task is to create a Python script that generates a parametric structural canopy within **Grasshopper** using the **GhPython** component.

Your script should:

- Implement functions that manipulate a base surface to create depth variations using mathematical functions.
- Generate vertical supports using recursive geometry to mimic branching structural elements.
- Combine different strategies for surface tessellation to achieve non-uniform patterns.
- Allow customization of parameters such as angles, length scaling factors, recursion depth, and tessellation strategies.
- Produce at least **three distinctly different outputs** for each of the following categories by varying parameters or methods:
  - (a) **Base Surface Shape**
  - (b) **Surface Tessellation Pattern**
  - (c) **Vertical Supports**

---

## Table of Contents

- [Tasks](#tasks)
- [Getting Started](#getting-started)
- [Repository Structure](#repository-structure)
- [Instructions](#instructions)
- [Submission Guidelines](#submission-guidelines)
- [Evaluation Criteria](#evaluation-criteria)
- [Resources](#resources)
- [Contact](#contact)

---

## Tasks

1. **Implement Parametric Canopy Script**

   - Create a Python script (`parametric_canopy.py`) that generates a parametric structural canopy using Grasshopper's Python component.
   - Use **Rhino.Geometry** and **RhinoScriptSyntax** to manipulate geometric objects within Rhino.
   - Define the parameters of your canopy explicitly in your script:
     - **Base Surface**: The initial surface for the canopy.
     - **Depth Map Control**: Parameters controlling the depth variation of the canopy.
     - **Recursion Parameters**: Parameters controlling the recursive generation of supports.
     - **Tessellation Strategy**: The method used for surface tessellation.
     - **Support Points**: Points where vertical supports will be generated.
   - Allow customization of these parameters to enable exploration of different designs.
   - **Optional:** Incorporate randomness to achieve greater visual complexity.

2. **Generate Distinct Outputs**

   - Use your script to produce at least **three distinctly different outputs** for each category by varying parameters or methods:
     - (a) **Base Surface Shape**
     - (b) **Surface Tessellation Pattern**
     - (c) **Vertical Supports**
   - Save each output as an image in the `images/` folder with descriptive filenames.

3. **Visualization**

   - Use Grasshopper's visualization tools to display your canopy designs.
   - Enhance your visualizations with custom color schemes, materials, or other visual effects.
   - Ensure that your visualizations effectively communicate the design intent and variations.

4. **Explore Advanced Concepts**

   - **Extend to Complex Tessellations**:
     - Experiment with advanced tessellation methods like Voronoi diagrams or Delaunay triangulation.
     - Handle non-planar tessellation challenges.
   - **Modify Geometric Rules**:
     - Incorporate randomness or use different geometric transformations.
     - Explore biomimicry by mimicking natural forms in your supports or canopy.

5. **Documentation**

   - Write detailed pseudo-code explaining your functions in `DOCUMENTATION.md`.
   - Provide technical explanations of your algorithms and design decisions.
   - Include a short report discussing the process and the mathematical principles behind your canopy designs.
   - Insert images of your generated canopies into the documentation.

## Getting Started

### Prerequisites

- **Rhinoceros 3D** with **Grasshopper** plugin
- **Python scripting in Grasshopper** (GhPython component)
- Familiarity with **Rhino.Geometry** and **RhinoScriptSyntax** modules
- **Python 3.x** (if using Rhino 8's Python 3 support)

### Software Installation

- Ensure that you have Rhino 8 installed, which includes Grasshopper.
- Install any additional plugins or libraries if required.

## Repository Structure

```
Assignment3/
├── README.md               # Assignment handout (this file)
├── DOCUMENTATION.md        # Documentation with pseudo-code and explanations
├── parametric_canopy.py    # Python script for parametric canopy generation
├── parametric_canopy.gh    # Grasshopper definition file
├── images/
│   └── (Output Images)
└── src/
    └── (Additional scripts or components)
```

---

## Instructions

1. **Clone the Repository**

   - Accept the assignment on GitHub Classroom and clone the repository to your local machine using GitHub Desktop or the command line.

2. **Set Up Your Environment**

   - Open Rhino and launch Grasshopper.
   - Familiarize yourself with the `parametric_canopy.py` script and the repository structure.

3. **Implement the Parametric Canopy Script**

   - Open the GhPython component in Grasshopper.
   - Copy the template code from `parametric_canopy.py` into the GhPython editor.
   - **Define the Parameters** explicitly in your script:
     - **Base Surface**: Use planes, curves, or surfaces as starting points.
     - **Depth Map Control**: Implement mathematical functions (e.g., sine, cosine) to manipulate the surface.
     - **Recursion Parameters**: Control the branching patterns of supports (e.g., angle, length reduction, depth).
     - **Tessellation Strategy**: Choose between quadrilateral grids, triangular meshes, Voronoi patterns, etc.
     - **Support Points**: Define where the supports connect to the canopy.
   - **Implement the Functions**:
     - Write a function to generate the depth map on the base surface.
     - Write a function to tessellate the modified surface based on the chosen strategy.
     - Write a recursive function to generate the vertical supports.
   - **Visualization within Grasshopper**:
     - Use Grasshopper components to visualize your geometry.
     - Ensure that the outputs are properly connected and displayed.

4. **Customize Parameters**

   - Allow users to customize parameters through Grasshopper sliders or input fields.
   - Provide examples of different parameter sets in your script or as preset configurations.
   - Encourage exploration by making it easy to adjust parameters.

5. **Generate and Save Outputs**

   - Run your script with different parameter sets to generate at least three distinct variations for each category.
   - Save each output image in the `images/` folder with descriptive filenames (e.g., `canopy_variation1.png`).
   - Arrange the images in a grid format for easy comparison.

6. **Explore Advanced Concepts (Optional for Extra Points)**

   - **Implement Advanced Tessellation**:
     - Use Grasshopper plugins or scripts to generate Voronoi or Delaunay tessellations.
   - **Incorporate Randomness**:
     - Introduce randomness in angles, lengths, or tessellation patterns to create more organic designs.
   - **Extend to Complex Forms**:
     - Experiment with double-curved surfaces or free-form geometries.

7. **Documentation**

   - Write detailed pseudo-code explaining your functions in `DOCUMENTATION.md`.
   - Provide technical explanations of your algorithms and design decisions.
   - Include images and explanations of each of your design variations.
   - Discuss any challenges faced and how you overcame them.

8. **Version Control with Git**

   - Make regular commits with meaningful messages.
   - Push your commits to the GitHub repository.

9. **Review and Finalize**

   - Ensure your code is well-commented and organized.
   - Verify that all generated images are saved correctly.
   - Double-check your documentation for clarity and completeness.

---

## Submission Guidelines

- **What to Submit**:
  - Your `parametric_canopy.py` script and any additional scripts in the `src/` directory.
  - Your Grasshopper definition file `parametric_canopy.gh`.
  - Completed `DOCUMENTATION.md` with pseudo-code, explanations, and images.
  - Generated images saved in the `images/` folder.
  - Any additional documentation or resources in the `docs/` folder.

- **Submission Checklist**:
  - [ ] Code runs without errors in Grasshopper.
  - [ ] Code is well-commented and follows best practices.
  - [ ] `DOCUMENTATION.md` is complete and thorough.
  - [ ] At least three variations for each category are saved and referenced in your documentation.
  - [ ] All changes are committed with meaningful messages.
  - [ ] All commits are pushed to your GitHub repository.

---

## Evaluation Criteria

- **Implementation and Functionality**
  - Correct implementation of the parametric canopy script in Grasshopper Python.
  - Effective use of depth map control and recursive geometry functions.
  - Successful combination of surface tessellation strategies.
  - Generation of the structural system (shell/gridshell and vertical supports).

- **Design Exploration**
  - Variety and creativity in design variations.
  - Quality and clarity of visualizations.

- **Technical Understanding**
  - Demonstrated understanding of parametric modeling, recursion, and tessellation.
  - Clear explanations in code comments and documentation.

- **Code Quality**
  - Clean, readable, and well-organized code.
  - Use of meaningful variable names and proper formatting.

- **Documentation**
  - Detailed pseudo-code and technical explanations.
  - Inclusion and discussion of generated design variations.
  - Discussion of challenges and solutions.

- **Use of Git and Version Control**
  - Regular commits with meaningful messages.
  - Proper repository structure and organization.

---

## Resources

- **Python in Rhino and Grasshopper**
  - [Rhino.Python Guides](https://developer.rhino3d.com/guides/rhinopython/)
  - [RhinoScriptSyntax](https://developer.rhino3d.com/api/RhinoScriptSyntax/)

- **Surface Tessellation Techniques**
  - [Voronoi Diagrams](https://en.wikipedia.org/wiki/Voronoi_diagram)
  - [Delaunay Triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation)
  - [Diagrid Structures](https://en.wikipedia.org/wiki/Diagrid)

---

## Contact

If you have any questions or need assistance, please reach out to the instructor via email or the course forum.

---