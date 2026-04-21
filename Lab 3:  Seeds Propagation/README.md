# Lab 3: Seeds Propagation & Stereo Matching

## Overview

This project implements a **seed-based stereo matching algorithm** to compute dense disparity maps from a pair of stereo images. The approach uses:

- **NCC (Normalized Cross-Correlation)** as a similarity metric for matching image patches
- **Seed selection** to identify reliable correspondences with high correlation
- **Priority queue-based propagation** to expand the disparity map from seeds to neighboring regions
- **Flexible disparity range** to handle different baseline configurations
- **3D reconstruction** visualization of the resulting stereo geometry

The algorithm works in three stages: computing a dense disparity map, selecting high-confidence seeds, and propagating those seeds to neighboring pixels using a greedy approach guided by correlation quality.

---

## Features

### 1. Normalized Cross-Correlation (NCC) Matching

The function `ccorrel` computes the **normalized cross-correlation** between two image patches:

- Takes two images and two patch centers: $(i_1, j_1)$ in image 1 and $(i_2, j_2)$ in image 2.
- Extracts patches of size $(2 \cdot \text{win} + 1) \times (2 \cdot \text{win} + 1)$ where:
  ```cpp
  static const int win = (9-1)/2;  // win = 4, so patches are 9x9
  ```
- Computes mean-centered correlation:
  
  ```math
    \text{NCC} = \frac{\sum (p_1 - m_1)(p_2 - m_2)}{\sqrt{\sum(p_1 - m_1)^2 \cdot \sum(p_2 - m_2)^2 + \varepsilon}}
  ```

  where:
  - $m_1, m_2$ are the patch means
  - $\varepsilon = 0.1$ (constant to avoid division by zero on uniform patches)

- Returns a normalized similarity score in $[-1, 1]$, where 1 indicates perfect correspondence.

---

### 2. Seed Selection

The function `find_seeds` identifies **reliable match seeds** that will be used as starting points for propagation:

- Scans every pixel in image 1 (excluding border pixels within `win` of the image edge).
- For each pixel, searches over the disparity range $[\text{dmin}, \text{dmax}]$.
- For each candidate disparity $d$, computes the NCC with the corresponding location in image 2: $(i_1 + d, j_1)$.
- Records the disparity $d$ that achieves the **best (maximum) NCC** at that location.
- If the best NCC exceeds the seed threshold:
  ```cpp
  static const float nccSeed = 0.95f;
  ```
  then:
  - Marks the pixel as a **seed** in the `seeds` image.
  - Stores the disparity in the `disp` map.
  - Adds the seed to a **priority queue** `Q`, ordered by NCC (highest first).

This produces a sparse but highly reliable set of disparity estimates.

**Progress indication:** The function prints percentage progress as it scans through the rows.

---

### 3. Priority Queue-Based Propagation

The function `propagate` expands the disparity solution from seeds to neighboring pixels:

- **Data structure:** Uses a `priority_queue<Seed>`, which orders seeds by NCC in descending order (highest confidence first).
  
  ```cpp
  bool operator<(const Seed& s1, const Seed& s2) {
      return (s1.ncc < s2.ncc);  // Max-heap behavior
  }
  ```

- **4-neighborhood:** Each seed can propagate to its 4 neighbors:
  ```cpp
  static const int dx[] = {+1,  0, -1,  0};
  static const int dy[] = { 0, -1,  0, +1};
  ```

- **Propagation step:** For each seed popped from the queue:
  1. Checks the 4 neighboring pixels.
  2. For each unseed neighbor at $(x, y)$:
     - Searches disparities in a **local window** around the parent's disparity:
       ```
       shifts: [-1, 0, +1]
       ```
       This leverages the assumption that disparity varies smoothly.
     - Computes NCC for each candidate disparity.
     - Selects the best disparity and its NCC.
  3. Marks the neighbor as seeded, stores its disparity, and pushes it onto the queue.

- **Greedy property:** By processing seeds in order of *decreasing* NCC, the algorithm prioritizes propagating from higher-confidence regions, leading to a more reliable dense map.

---

### 4. Disparity Map Visualization

The function `displayDisp` converts the disparity map into a color image for visualization:

- **Invalid disparities:** Pixels with disparity outside the range $[\text{dmin}, \text{dmax}]$ are colored **cyan**.
- **Valid disparities:** Mapped to a **grayscale gradient**:
  
  ```math 
    \text{gray} = 255 \cdot \frac{d - \text{dmin}}{\text{dmax} - \text{dmin}}
  ```
  where darker values represent disparities closer to `dmin` and lighter values closer to `dmax`.

- **Three outputs are saved:**
  1. `0dense.png`: Dense disparity (all pixels matched, ignoring NCC threshold).
  2. `1seeds.png`: Only seed pixels (high-confidence matches).
  3. `2final.png`: Final propagated disparity map (seeds + propagated neighbors).

The cyan border indicates regions where stereo matching could not be computed (e.g., occlusions or image boundaries).

---

### 5. 3D Reconstruction Visualization

The function `show3D` reconstructs and displays the scene in 3D space using the disparity map:

- **Prerequisites:** Requires Imagine++ built with **OpenGL support**.
- **Camera calibration:** Uses intrinsic parameters from the Middlebury stereo datasets:
  ```cpp
  const float f = 3740;        // Focal length in pixels
  const float B = 0.160;       // Baseline in meters
  const float d0 = -200;       // Disparity offset (image cropping)
  const float zoom = 2;        // Original images were 2x larger
  ```

- **3D reconstruction:** Converts each valid disparity to a 3D point using the camera model:
  
  ```math
    z = \frac{B \cdot f}{d_\text{scaled}}
  ```
  
  where $d_\text{scaled} = \text{zoom} \cdot d + d_0$.
  
  The $(x, y)$ coordinates are backprojected using the inverse intrinsic matrix $K^{-1}$.

- **Output:** A colored point cloud where each point is colored by the corresponding pixel in image 1.
- **Interaction:** The window opens in 3D space; use **shift-click** to animate and explore the reconstruction.

---

## Algorithm Parameters

The behavior of the algorithm is controlled by several constants:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `nccSeed` | 0.95 | NCC threshold for seed selection |
| `win` | 4 | Radius of correlation patch (9×9 patches) |
| `EPS` | 0.1 | Constant to avoid division by zero in NCC |
| `dmin` | -30 | Minimum disparity to search |
| `dmax` | -7 | Maximum disparity to search |

**Notes:**
- Negative disparities indicate rightward shifts (image 2 is shifted right relative to image 1).
- The disparity range $[\text{dmin}, \text{dmax}]$ can be adjusted via command-line arguments.
- A **tighter seed threshold** (closer to 1.0) produces fewer but more reliable seeds.
- **Larger patch windows** (`win`) are more robust to texture variations but slower and less precise.

---

## Program Workflow

1. **Initialization**
   - Load two stereo images from command-line arguments or default paths:
     ```cpp
     string im1 = argc > 1 ? argv[1] : srcPath("im1.jpg");
     string im2 = argc > 1 ? argv[2] : srcPath("im2.jpg");
     int dmin = argc > 3 ? stoi(argv[3]) : -30;
     int dmax = argc > 3 ? stoi(argv[4]) : -7;
     ```
   - Create images to store the disparity map, seed markers, and a priority queue.

2. **Dense Disparity (Baseline)**
   - Call `find_seeds` with `nccSeed = -1.0f` to match *every* pixel without filtering.
   - Save the result as `0dense.png`.
   - This provides a reference for the quality of seed-based approaches.

3. **High-Confidence Seed Selection**
   - Call `find_seeds` with `nccSeed = 0.95f` to select only reliable seeds.
   - Save the sparse seed map as `1seeds.png`.
   - Populate the priority queue with these seeds.

4. **Propagation**
   - Call `propagate` to expand the disparity map from seeds to all reachable regions.
   - The priority queue ensures propagation is guided by match confidence.
   - Save the final dense map as `2final.png`.

5. **3D Reconstruction (Optional)**
   - Call `show3D` to visualize the stereo pair as a 3D point cloud.
   - If OpenGL is available, the user can interact with the 3D model.
   - If OpenGL is unavailable, a message is printed to the console.

6. **Cleanup**
   - Close the graphics window and exit.

---

## Dependencies

- **Imagine++**:
  - Graphics: image I/O, window management, 2D/3D display.
  - Linear algebra: 3×3 matrices (`FMatrix`), vector operations.
  - 3D visualization (optional): OpenGL integration for 3D point cloud rendering.
- **Standard C++**:
  - `<queue>`: priority queue for seed management.
  - `<string>`, `<iostream>`: file and console I/O.
  - `<cmath>`: numerical operations (sqrt, pow).

---

## Building the Project

### Using CMake

A `CMakeLists.txt` is provided to build the project using CMake:

```bash
mkdir cmake-build-debug
cd cmake-build-debug
cmake ..
make
```

### Manual Compilation

Adjust paths to match your Imagine++ installation:

```bash
g++ Seeds.cpp -o Seeds -std=c++11 \
    -I/path/to/Imagine++/include \
    -L/path/to/Imagine++/lib \
    -lImagineGraphics -lImagineImages -lImagineLinAlg
```

---

## Running the Program

### Default Images

```bash
./Seeds
```

Loads default images `im1.jpg` and `im2.jpg` from the Imagine++ resource path, with disparity range $[-30, -7]$.

### Custom Images and Disparity Range

```bash
./Seeds im1.jpg im2.jpg dmin dmax
```

Example:
```bash
./Seeds left.jpg right.jpg -60 -10
```

### Output Files

The program saves three PNG images to the Imagine++ resource directory:
- `0dense.png`: Dense disparity without seed filtering.
- `1seeds.png`: Sparse high-confidence seeds only.
- `2final.png`: Final propagated disparity map.

---

## Usage Steps

1. **Prepare stereo images:**
   - Rectified stereo pairs work best (epipolar lines are horizontal).
   - Images should be roughly the same size.

2. **Run the program:**
   ```bash
   ./Seeds im1.jpg im2.jpg dmin dmax
   ```

3. **Observe the results:**
   - **Window 0 & 1:** Input images (image 1 and image 2).
   - **Window 2:** Dense disparity map (all pixels with high or low NCC).
   - **Window 3:** Seed map (only high-confidence NCC ≥ 0.95).
   - **Window 4:** Final propagated disparity map.
   - The cyan border indicates regions where matching failed.

4. **Inspect 3D reconstruction (if OpenGL available):**
   - A 3D window displays the reconstructed scene as a point cloud.
   - Use **shift-click** to animate/rotate the view.
   - Close the 3D window to exit.

5. **Check saved outputs:**
   - View `0dense.png`, `1seeds.png`, and `2final.png` to compare stages.

---

## Example Results

The three output images show the progression of the algorithm:

1. **Dense map (0dense.png):**
   - Computed without any NCC filtering.
   - Shows all matched disparities, including unreliable ones.
   - Likely contains noise and outliers.

2. **Seeds (1seeds.png):**
   - Only pixels with NCC ≥ 0.95 are marked.
   - Creates a sparse but very reliable set of seed points.
   - Cyan areas indicate regions not selected as seeds.

3. **Final map (2final.png):**
   - Propagation from seeds fills in neighboring regions.
   - Produces a denser map than seeds alone.
   - Cyan borders remain where propagation could not reach (occlusions, boundaries).

---

## Limitations & Possible Improvements

**Current limitations:**

- **No occlusion handling:** Occluded regions are colored cyan and left empty.
- **No left-right consistency check:** No validation that left-to-right and right-to-left disparities agree.
- **Smooth surfaces only:** The local disparity window (±1 pixel) assumes smooth depth variation; steep edges may be poorly estimated.
- **Fixed parameters:** `nccSeed`, `win`, and disparity range are hardcoded; no adaptive tuning.
- **No sub-pixel refinement:** Disparity estimates are integer-valued.
- **OpenGL dependency:** 3D visualization requires OpenGL; headless systems cannot visualize 3D points.

**Possible extensions:**

- **Left-right consistency check:** Ensure disparities are consistent in both directions; fill holes using the opposite direction.
- **Adaptive seed threshold:** Automatically select `nccSeed` based on image statistics.
- **Edge-preserving propagation:** Respect image boundaries and avoid propagating across high-contrast regions.
- **Sub-pixel disparity:** Refine disparity estimates using parabolic fitting or other sub-pixel methods.
- **Multi-scale approach:** Process at coarse resolution first, then refine at finer scales.
- **Occlusion estimation:** Detect and mark occluded regions explicitly.
- **Save intermediate data:** Export disparity maps and 3D point clouds in standard formats (e.g., PLY, VTK).

---

## Credits

This project demonstrates a practical implementation of **seed propagation** for stereo matching, a classical approach in computer vision that balances speed and accuracy. The algorithm is inspired by methods used in real-time stereo systems and provides insight into the trade-off between robustness (seed selection) and density (propagation).

Key concepts:
- **Normalized Cross-Correlation (NCC)** for similarity metric
- **Priority-queue-based propagation** for greedy dense estimation
- **Stereo geometry** and disparity-to-depth conversion
- **3D reconstruction** from stereo pairs

This implementation is part of a practical session on **multi-view 3D reconstruction**.