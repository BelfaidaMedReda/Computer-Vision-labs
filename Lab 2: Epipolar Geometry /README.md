# Lab 2 : Fundamental Matrix Estimation & Epipolar Geometry

## Overview

This project estimates the **fundamental matrix** between two views of the same scene using:

- Automatic **SIFT feature detection and matching**
- The **normalized 8-point algorithm** (with Hartley normalization)
- A robust **RANSAC** scheme to reject outliers
- An interactive visualization of **epipolar geometry**

The program takes two images as input, finds point correspondences automatically, estimates the fundamental matrix \(F\), filters inliers, and then lets the user click in one image to display the corresponding epipolar lines in the other.

---

## Features

### 1. Automatic SIFT Feature Matching

The function `algoSIFT`:

- Detects SIFT keypoints in both images using `SIFTDetector` from Imagine++.
- Draws the detected features on each image for visualization.
- Computes pairwise descriptor distances between features in image 1 and image 2.
- Keeps all pairs whose squared descriptor distance is below a fixed threshold:

  ```cpp
  const double MAX_DISTANCE = 100.0 * 100.0;
  ```

- Stores matches as `Match {x1, y1, x2, y2}` where:
  - $x_1, y_1$ is the point in the first image
  - $x_2, y_2$ is the corresponding point in the second image

This produces an initial (possibly noisy) set of candidate correspondences.

---

### 2. Hartley Normalization of Points

To improve the numerical stability of the 8-point algorithm, the code implements **Hartley normalization**:

- Computes the mean of the coordinates in each image.
- Computes the average distance of points to their centroid.
- Derives a scale factor so that the average distance becomes $\sqrt{2}$.
- Builds normalization matrices $T$ and $T'$ for the first and second image.

Each match $x_1, y_1, x_2, y_2$ is then transformed:

- $\tilde{x} = T \cdot (x_1, y_1, 1)^T$
- $\tilde{x}' = T' \cdot (x_2, y_2, 1)^T$

The function `applyHartleyNormalization`:

- Normalizes all matches **in-place**.
- Returns the pair of normalization matrices $(T, T')$ so they can later be used to denormalize the fundamental matrix.

---

### 3. Fundamental Matrix Estimation (Normalized 8-Point Algorithm)

The function `estimateF` implements the **normalized 8-point algorithm**:

1. **Normalization**  
   Calls `applyHartleyNormalization` to normalize the points and retrieve matrices $T$ and $T'$.

2. **Building the linear system**  
   For each correspondence, it builds one row of the matrix $A$ such that:

   
   $$A \cdot f = 0$$
   

   with $f$ being the 9-vector of the entries of $F$. Each row is of the form:

   ```cpp
   [x1*x2, x1*y2, x1,
    y1*x2, y1*y2, y1,
    x2,    y2,    1]
   ```

3. **SVD to find the least-squares solution**  
   Performs SVD on $A$ and takes the singular vector corresponding to the smallest singular value as the solution $f$. This vector is reshaped into a $3 \times 3$ matrix $F$.

4. **Enforcing rank-2 constraint**  
   - Performs SVD on $F$: $F = U \Sigma V^T$.
   - Sets the smallest singular value in $\Sigma$ to zero.
   - Reconstructs the rank-2 fundamental matrix:

    $$ F' = U \Sigma' V^T$$

5. **Denormalization**  
   Since the estimation was done on normalized points, the fundamental matrix is denormalized:

   $$ F = T^T \cdot F' \cdot T' $$

The resulting matrix $F$ satisfies the epipolar constraint $x'^T F x = 0$ for corresponding points.

---

### 4. Robust Estimation with RANSAC

The function `computeF` wraps the 8-point algorithm in a **RANSAC** loop to handle outliers:

- Uses all SIFT matches as the initial set of observations.
- At each iteration:
  1. Randomly selects 8 matches to form a minimal sample.
  2. Calls `estimateF` to compute a candidate fundamental matrix $F$.
  3. For every match:
     - Computes the epipolar line in the second image: $ l' = F^T x$.
     - Computes the distance of the point $x'$ to the line $l'$:

       $$ d = \frac{|x'^T l'|}{\|l'\|}$$

     - Classifies the correspondence as **inlier** if `d <= distMax` with:

       ```cpp
       const float distMax = 1.5f; // pixel threshold
       ```

  4. Keeps the model with the largest inlier set.

- The number of RANSAC iterations is **adaptively updated** based on the current inlier ratio $ \hat{w} $, using the theoretical formula:

  $$ N_\text{necessary} = \frac{\log(\text{BETA})}{\log(1 - \hat{w}^8)} $$

  where `BETA` is the failure probability (`0.01` in the code).

- Once the best inlier set is found:
  - `matches` is **filtered in-place** to keep only the inliers.
  - `estimateF` is called again on this refined set to obtain the final fundamental matrix.

---

### 5. Visualization of Matches and Inliers

After estimating $F$:

- The program displays both images side-by-side.
- For each remaining match (RANSAC inliers):
  - Draws a small colored circle at the correspondence position in both images.
  - Uses the same random color for both points of a match for easy identification.
- Displays the total number of inliers vs. original matches.

This gives an immediate visual check of the quality of the robust estimation.

---

### 6. Interactive Epipolar Geometry

The function `displayEpipolar` lets the user explore the epipolar geometry:

- The user clicks with the **left mouse button** in either image:
  - A small red circle marks the clicked point.
  - If the click is on the **left image**:
    - The corresponding epipolar line is drawn in the **right image**.
  - If the click is on the **right image**:
    - The corresponding epipolar line is drawn in the **left image**.
- Lines are computed using:
  - $l' = F^T x$ for a point $x$ in the first image.
  - $l = F x'$ for a point $x'$ in the second image.
- The user can **right-click** to exit the epipolar exploration.

This interactive tool provides a concrete visualization of the epipolar constraint.

---

## Program Workflow

1. **Initialization**
   - Seed the random number generator with `time(0)`.
   - Load two images from command-line arguments if provided, otherwise from default paths:
     ```cpp
     string s1 = argc>1? argv[1]: srcPath("im1.jpg");
     string s2 = argc>2? argv[2]: srcPath("im2.jpg");
     ```

2. **Display**
   - Open a window wide enough to display both images side-by-side.
   - Draw the images in the window.

3. **SIFT Matching**
   - Call `algoSIFT` to detect features and build initial matches.
   - Display the number of matches and wait for a click to continue.

4. **Fundamental Matrix Estimation with RANSAC**
   - Call `computeF` to:
     - Filter matches to inliers.
     - Estimate the final fundamental matrix \(F\).
   - Print \(F\) to the console.

5. **Inlier Visualization**
   - Redisplay the images.
   - Draw inlier matches with random colors.
   - Display the inlier count and wait for a click.

6. **Epipolar Line Exploration**
   - Redisplay the images.
   - Call `displayEpipolar` to interactively inspect epipolar lines.
   - Exit on right-click, then close the graphics window.

---

## Dependencies

- **Imagine++**:
  - Graphics: image display, drawing primitives, mouse input.
  - Linear algebra: `Matrix`, `Vector`, `FMatrix`, `FVector`, SVD.
  - SIFT feature detection.
- **Standard C++**:
  - `<vector>`, `<cmath>`, `<cstdlib>`, `<ctime>`, `<random>`, `<iostream>`.

---

## Building the Project

### Using `make`

Create a `Makefile` similar to other Imagine++ examples:

```make
all: fundamental

fundamental: fundamental.cpp
	$(CXX) fundamental.cpp -o fundamental -std=c++11 	    -I/path/to/Imagine++/include 	    -L/path/to/Imagine++/lib 	    -lImagineGraphics -lImagineImages -lImagineLinAlg
```

Adjust the paths to match your Imagine++ installation.

Then compile:

```bash
make
```

### Manual Compilation Example

```bash
g++ fundamental.cpp -o fundamental -std=c++11     -I/path/to/Imagine++/include     -L/path/to/Imagine++/lib     -lImagineGraphics -lImagineImages -lImagineLinAlg
```

---

## Running the Program

```bash
./fundamental im1.jpg im2.jpg
```

If no image paths are provided, the program tries to load default images `im1.jpg` and `im2.jpg` from the Imagine++ resource path.

---

## Usage Steps

1. **Start the program** with two views of the same static scene.
2. **Observe SIFT features and initial matches** (first window, then click to continue).
3. **Wait for RANSAC** to estimate the fundamental matrix and filter inliers.
4. **Inspect the filtered matches** shown with random colors.
5. **Click to enter epipolar exploration mode**.
6. **Left-click** points in either image to see the corresponding epipolar lines in the other.
7. **Right-click** to exit and close the program.

---

## Limitations & Possible Improvements

- Descriptor matching uses a **simple distance threshold**, no ratio test (e.g., Lowe’s ratio).
- RANSAC uses **uniform random sampling**; no guided sampling or more advanced robust estimators (LMedS, MLESAC, etc.).
- The epipolar distance is computed in a **single image** using a symmetric-like distance, but alternative robust error measures could be added.
- No explicit handling of degenerate cases (e.g., planar scenes, too few inliers).
- No saving of the estimated \(F\) or visualization results to disk.

Possible extensions:

- Implement Lowe’s ratio test for more reliable initial matches.
- Draw **epipoles** and epipolar lines for multiple random points.
- Add an option to **save the fundamental matrix** to a file.
- Visualize **Sampson distance** or other robust geometric errors.

---

## Credits

This project was developed as part of a practical session on **multiple-view geometry**, focusing on:

- Automatic **feature detection and matching**
- **Fundamental matrix** estimation
- **RANSAC**-based robust fitting
- Interactive visualization of **epipolar geometry**