# Neural-Oscillations-Cpp
The c++ code for neural oscillation project.

## Way to run the reduced and coarse grained models

* step1: Run **model_full.cpp** and generate **membrane_potential_sample.txt**. 
  * Set a large value for *terminate_time* is recommended to get a stable probability.
  * You will get other two txt files as by-product.
* step2: Run **P_generation.m** in matlab and generate **Probability.txt**.
* step3: Now you can run **model_reduced_network.cpp** and **model_coarse_grained.cpp**.

