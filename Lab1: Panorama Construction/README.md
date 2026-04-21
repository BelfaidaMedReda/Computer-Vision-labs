# Lab 1 : Panorama Construction

## Overview
This project demonstrates how to build a **panorama from two images** using the **Imagine++** library. The program allows the user to manually select point correspondences between two views, computes the homography relating them, and merges the images into a unified coordinate frame.

---

## Features

### 1. Interactive Point Selection
Two windows display the input images. The user selects matching points between the two views:

- **Left click**: registers a point in the active window.
- **Right click**: ends the point selection process.

Corresponding points are stored in synchronized vectors for later processing.

---

### 2. Homography Estimation
With at least four point correspondences, the program computes a projective transformation:

- Builds the linear system \( A h = B \) from the projective mapping constraints.
- Solves the system using `linSolve`.
- Reconstructs the **3×3 homography matrix**.
- Performs a brief consistency check by evaluating residuals for each correspondence.

---

### 3. Panorama Construction
Once the homography is computed:

- The corners of the first image are projected into the coordinate system of the second.
- A new output image is created large enough to contain both warped images.
- For each pixel in the panorama:
  - Check whether it lies inside image 2.
  - Compute its pre-image in image 1 using \( H^{-1} \).
  - Blend contributions when both images overlap (simple averaging).

The final panorama is displayed in a separate window.

---

## Dependencies
- **Imagine++** library
- Standard C++ libraries (iostream, vector, etc.)

---

## Building the Project

### Using `make`
Create a simple Makefile based on the provided Imagine++ examples, then compile:

```bash
make
```

### Manual Compilation Example
```bash
g++ panorama.cpp -o panorama -std=c++11 \
    -I/path/to/Imagine++/include \
    -L/path/to/Imagine++/lib -lImagineGraphics -lImagineImages -lImagineLinAlg
```
Replace the paths with those of your Imagine++ installation.

---

## Running the Program
```bash
./panorama imageA.jpg imageB.jpg
```
If no arguments are provided, the program loads default images (`image0006.jpg` and `image0007.jpg`).

### Usage Steps
1. Select matching points in both windows using left clicks.
2. Right-click to validate the selection.
3. The homography is computed, and the resulting panorama is displayed.

---

## Tips for Better Results
- Choose points that are easy to identify (corners, intersections, distinct textures).
- Avoid selecting collinear or extremely close points.
- Use more than four correspondences for improved robustness.

---

## Limitations
- No outlier filtering or robustness checks.
- Pixel lookup uses direct access without bilinear interpolation.
- Blending is limited to simple averaging.
- No exposure compensation or seam optimization.

---

## Credits
This project was developed as part of a practical session on projective geometry, homographies, and image stitching.

