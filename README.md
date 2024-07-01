<h1 align="center"> EnPrO: Enhancing Precision Through
Optimization in Image-Guided Spine Surgical
Procedures</h1>

<p  align="center">  
  
 Advancements in intraoperative visualization have driven the adoption of fluoroscopy-based imaging in `Image-Guided Surgery Systems (IGSS)` for spine procedures. Precise localization and identification of fiducials are crucial for IGSS accuracy. This study enhances IGSS precision by optimizing the mapping between `2D fluoroscopic images` and `patient anatomy`, minimizing `reprojection error (RPE)` by selecting `optimal fiducial points`. The Least-Squares Method is used to compute weights for the fiducials, excluding those that significantly contribute to RPE, thus reducing error and the number of fiducials required. This optimization is beneficial in scenarios where certain fiducials are not detected due to occlusion or low contrast. By selecting the most effective fiducials, IGSS accuracy is improved, ensuring better surgical outcomes.
</p>

<h3 > <i>Index Terms</i> </h3> 

 :diamond_shape_with_a_dot_inside:Image Guided Spine Surgery
  :diamond_shape_with_a_dot_inside: Fiducials
  :diamond_shape_with_a_dot_inside:Reprojection Error

</div>


## <div align="center">Getting Started</div>

<details>
  <summary><i>Fiducial Extraction</i></summary>

  - Calibration drums have over 80 fiducials. Detection and identification involve image processing steps like filtering and connected component analysis.
</details>
<details open>
<summary><i>Study of Fiducials</i></summary>
  
  - Analyzed fiducials' contribution to RPE, comparing errors before and after removing individual fiducials.

</details>

<details open>
<summary><i>Optimization for Fiducial Weights</i></summary>
  
  - Using fiducial coordinates and RPE to determine optimal fiducials, employing the Least Squares method for real-time implementation.

</details>
<details open>
  <summary><i>Camera Calibration</i></summary>
  
- The fiducials chosen based on their weights are used for calibration. The reprojection
error obtained when all 17 fiducials are used in calibration is compared with the
errors obtained when the top ten fiducials and the bottom ten fiducials, along with the
significant fiducials, are used for calibration. The fiducials that give the minimum error
are the optimal fiducials.
</details>

## <div align="center">Methodology</div>

<p align="center">
  <img src="Block Diagram.png">
</p>
<div align = "center">
  
  :small_orange_diamond: Fig 4:Block diagram of proposed work - EnPrO
</div>


## <div align="center">Pre-requisites</div>

Before installing and running the project, ensure you have the following prerequisites:

 :grey_exclamation: Download and install CVX matlab toolbox from the [CVX](https://cvxr.com/cvx/)
 
  :grey_exclamation: Download and install MATLAB.

  
## <div align="center">Installation</div>
:arrow_right:Clone the Repository
```bash
git clone https://github.com/Deepa-Ramki/Fiducial-Optimization.git
```

:arrow_right:Navigate to the Project Directory
```bash
cd Fiducial-Optimization
```
:arrow_right:Install Dependencies mentioned in Pre-requisites

## <div align="center">Environments</div>

<div align="center">
  <a href="https://www.mathworks.com/products/matlab.html">
    <img src="https://upload.wikimedia.org/wikipedia/commons/2/21/Matlab_Logo.png" width="10%" /></a>
</div>

<div align="center">
  
**This repository contains the codebase. Please note that sample data is not included yet but will be made available soon.**
</div>
