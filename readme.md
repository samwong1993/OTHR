#Code for Over-the-horizon Radar  
Published by Sen Huang  
Email: hsen@se.cuhk.edu.hk  
IDE: MATLAB R2019a  
The convergence of the algorithm has been proved. Thanks for help of Prof. Anthony So.

##Function  
* main.m: You can test our algorithm.
* MonteCarlo.m: Monte Carlo Simulation for the real-world data/sythetic data.
* ABC.m: Return value of A/B/C.
* beta_bound.m: Return the EDR of the sensors by golden section search.
* CRLB_tdoaOTHR.m: Calculate the Crame-Rao Lower Bound.
* LGLTtoXYZ.m/XYZtoLGLT.m: Transformation between Latitude/Longitude and Euclidean coordinate system.
* solve_eq.m: Recover the feasibility of beta for corresponding x.
* UpBound_beta.m: Return the upper bound of beta.
* proj_ball.m: Projection on a ball.
* PD.m/graPD.m/hessPD.m/GraF.m: 
    * 1. P/D value. 
    * 2. Gradient of P/D on beta. 
    * 3. Hessian of P/D on beta.
    * 4. Gradient of objective function on beta.  
* generator.m: Generate sythetic simulation data.  
* raypath.m: Plot the ray path.  
* Initialization.m: Our initialization strategy.  
* hildreth.m: Projection onto polyhedron.
* GPGD.m: Our proposed method.
* pltMCCRLB.m: The Monte Carlo simulation can be plotted directly by pltMCCRLB('realdata_log.txt','*k-').
* plotobj.m: Plot the figure 2 in the paper.
* Set plt =1 to visualize the simulation process (figure 4).
* perturbation.m: Plot figure 7.
* plotCRLB.m: An demo for pltMCCRLB.m (figure 5, figure 6).  
##Base Line
* ICASSP.m: T. Wang's method (delta is 0.1).
* CD.m:Coordinate descent  
##Real-World Data  
Here we give three real-world settings. beta0 is the true flying angle and each case has five sensors.  
Real-world data 1:  
beta0 = [0.133647283541992,0.100329606166839,0.335342045407703,0.413523978325060,0.0846681067915992];  
XYZ = zeros(M,3);  
%Xi Ning  
[x0 y0 z0] = LGLTtoXYZ(101.74,36.56,R);  
emitter = [x0 y0 z0]';  
%Hang Zhou  
[x0 y0 z0] = LGLTtoXYZ(120.19,30.26,R);  
XYZ(1,:) = [x0 y0 z0];  
%Wen Chang  
[x0 y0 z0] = LGLTtoXYZ(110.72,19.61,R);  
XYZ(2,:) = [x0 y0 z0];  
%Zheng Zhou  
[x0 y0 z0] = LGLTtoXYZ(113.65,34.76,R);  
XYZ(3,:) = [x0 y0 z0];  
%Tai Yuan  
[x0 y0 z0] = LGLTtoXYZ(112.53,37.87,R);  
XYZ(4,:) = [x0 y0 z0];  
%Seoul  
[x0 y0 z0] = LGLTtoXYZ(126.58,37.33,R);  
XYZ(5,:) = [x0 y0 z0];  
Real-world data 2:  
beta0 = [0.114957231412252,0.449398124172348,0.277420425918117,0.0168095219080640,0.103488345084960];  
XYZ = zeros(M,3);  
%Hong Kong  
[x0 y0 z0] = LGLTtoXYZ(114.16,22.28,R);  
emitter = [x0 y0 z0]';  
%Bei Jing   
[x0 y0 z0] = LGLTtoXYZ(116.41,39.90,R);  
XYZ(1,:) = [x0 y0 z0];  
%Wu Han  
[x0 y0 z0] = LGLTtoXYZ(114.31,30.59,R);  
XYZ(2,:) = [x0 y0 z0];  
%Shang Hai  
[x0 y0 z0] = LGLTtoXYZ(121.47,31.23,R);  
XYZ(3,:) = [x0 y0 z0];  
%Tokyo  
[x0 y0 z0] = LGLTtoXYZ(139.69,35.69,R);  
XYZ(4,:) = [x0 y0 z0];  
%Seoul  
[x0 y0 z0] = LGLTtoXYZ(126.58,37.33,R);  
XYZ(5,:) = [x0 y0 z0];  
Real-world data 3:  
beta0 = [0.349025176895742,0.279243732945529,0.174632834143026,0.183395025597251,0.285385790728298];  
or  
beta0 = [0.674175558279844,0.678354140357644,0.679512211766613,0.679496774102502,0.678155318189844];  
XYZ = zeros(M,3);  
[x0 y0 z0] = LGLTtoXYZ(116.24,39.55,R);  
emitter = [x0 y0 z0]';  
[x0 y0 z0] = LGLTtoXYZ(128.72,40.55,R);  
XYZ(1,:) = [x0 y0 z0];  
[x0 y0 z0] = LGLTtoXYZ(130.42,38.68,R);  
XYZ(2,:) = [x0 y0 z0];  
[x0 y0 z0] = LGLTtoXYZ(132.94,33.82,R);  
XYZ(3,:) = [x0 y0 z0];  
[x0 y0 z0] = LGLTtoXYZ(130.90,31.84,R);  
XYZ(4,:) = [x0 y0 z0];  
[x0 y0 z0] = LGLTtoXYZ(129.06,35.63,R);  
XYZ(5,:) = [x0 y0 z0];  
