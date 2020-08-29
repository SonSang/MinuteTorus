# MinuteTorus

<p align="center">
  <img src="Image/main.png">
</p>

**C++ library that supports basic math operations related to torus :  Computing Osculating Toroidal Patch, Computing Binormal Lines between Toroidal Patches, Computing Gaussmaps of Toroidal patches, etc.**

Please refer to following papers to learn about details of this library :

[1] Sang-Hyun Son, Seung-Hyun Yoon, Myung-Soo Kim, Gershon Elber, *Efficient Minimum Distance Computation for Solids of Revolution*, Eurographics & Eurovis (Computer Graphics Forum 2020)

# Circle
Circles play important roles in geometric computations related to tori : they offer efficient algorithm to find local extremums of distance between different tori. Therefore, this library offers algorithms to find **binormals - line that passes through two geometric entities orthogonally at the same time -** between two circles in 3D space. Also, it offers an algorithm to find **minimum distance** between two circles in 3D space, which is much faster than finding every binormal, as minimum distance is just the shortest binormal.

# Torus
