# Anisotropically Scaled Iterative Closest Point Algorithm

### The Project
Implementation of the Anisotropically Scaled Iterative Closest Point Algorithm done in C++ using the Eigen library.
The algorithm solves for the transformation between two sets of points. This includes rotation, scaling, and translation.
And unlike the normal ICP algorithm it will solve for transforamtions that have and inherent anisotropic scaling relationship.

The algorithm used is as described by the Chen et al. paper which is mostly the same as the normal ICP solution but instead
uses Mahalanobis distance as a metric for finding point correspondence instead of Eculidean distance.
It also uses an iterative Orthogonal Procrustes Analysis solution for finding the transformations from the corresponding points as given by Dosse and Berge
instead of the more conventional closed form solutions as it requires an algorithm that solves for anisotropic scaling.

### The Algorithm
Given two sets of Points X and Y which can be of different sizes(generally more in Y than X)
find the transformation that minimises the difference between the transformed X and Y.

The ASICP algorithm is as follows:
1. Translate points to be centred at the origin
2. Find the Correspondence of points in X and Y using Mahalanobis distance
3. Find the Transformation between the points selected in X and Y using ASOPA
4. Calculate error between the the transfomred X points and Y

### Dependencies
* Eigen Library (tested with 3.3.7)
* C++ compiler (tested with g++ 8 and 9 and clang++ 6)

### Usage
1. ``$ make``
2. ``$ ./asicp``

### References
* Elvis C.S. Chen, A. Jonathan McLeod, John S.H. Baxter, Terry M. Peters,
  "Registration of 3D shapes under anisotropic scaling"
  International Journal of Computer Assisted Radiology and Surgery
  June 2015, Volume 10, Issue 6, pp 867–878
  
* Mohammed Bennani Dosse and Jos Ten Berge (2010),
  "Anisotropic Orthogonal Procrustes Analysis"
  Journal of Classification 27:111-128
