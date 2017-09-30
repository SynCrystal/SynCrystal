# SynCrystal

1. OVERVIEW

SynCrystal is a MATLAB toolbox to analyze atomic crystal images. For a given atomic crystal image, it contains several tools to identify grain boundary, crystal orientation, elastic deformation. A few examples of synthetic and real atomic crystal images are provided to illustrate how to use these tools.

2. INTRODUCTION

SSTmethod: A collection of MATLAB and MEX routines which applies 2D synchrosqueezed transforms proposed in [1][2] to analyze crystal images. Detailed description of this method is in [4]. An instruction for parameter tuning in the synchrosqueezed transform is in [3]. 

VarSSTmethod: A collection of MATLAB and MEX routines which applies a variational model to optimize the crystal analysis results provided by 2D synchrosqueezed transforms. Detailed description of this method is in [5]. 


[1] H. Yang and L. Ying. Synchrosqueezed wave packet transform for 2d mode decompo- sition. SIAM Journal on Imaging Sciences, 6(4):1979–2009, 2013.

[2] H. Yang and L. Ying. Synchrosqueezed curvelet transform for two-dimensional mode decomposition. SIAM Journal on Mathematical Analysis, 46(3):2052–2083, 2014.

[3] H. Yang. Statistical Analysis of Synchrosqueezed Transforms. Applied and Computational Harmonic Analysis, 2017.

[4] H. Yang, J. Lu, and L. Ying. Crystal image analysis using 2d synchrosqueezed transforms. arXiv:1402.1262 [math.NA], SIAM Multiscale Modeling and Simulation, 2015.

[5] J. Lu, B. Wirth and H. Yang. Compbining 2d synchrosqueezed wave packet transforms with optimization for crystal image analysis. Journal of the Mechanics and Physics of Solids, 2016.

[6] J. Lu and H. Yang. Phase Space Sketching for Crystal Image Analysis based on Synchrosqueezed Transforms. arXiv:1703.10877 [cond-mat.mtrl-sci]

3. INSTALLING SYNLAB

Run the file SetPath.m first. It will automatically add all the MATLAB codes to your MATLAB path and compile all MEX files. After this, you can run all demo codes to see how to use this tool box.

4. COPY RIGHT

SynCrystal is copyright reserved. For further information, please contact 
Contact information
Jianfeng Lu at jianfeng@math.duke.edu
Benedikt Wirth at Benedikt.Wirth@uni-muenster.de
Haizhao Yang at matyh@nus.edu.sg
Lexing Ying at lexing@math.stanford.edu