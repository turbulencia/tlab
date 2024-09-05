# Run Tlab

This is a summary of the workflow. The details of each tool and the corresponding input data should be read in the corresponding source file.

## Preprocessing Tools

| binary    | inputs                            | outputs |
| --------- | --------------------------------- | ------- |
|inigrid.x   | tlab.ini                         | grid |
|inirand.x   | tlab.ini, grid                      |  [flow,scal].rand.? |        
|iniflow.x |  tlab.ini, grid [,flow.rand.?]     | flow.ics.?
|iniscal.x |   tlab.ini, grid [,scal.rand.?]     |    scal.ics.?

## Simulation Tools

| binary    | inputs                            | outputs |
| --------- | --------------------------------- | ------- |
|dns.x      |tlab.ini, grid, flow.*.?, scal.*.? |   flow.*.?, scal.*.? |

## Postprocessing Tools

| binary    | inputs                            | outputs |
| --------- | --------------------------------- | ------- |
|averages.x | tlab.ini, grid, flow.*.?, scal.*.?| avg*
|pdfs.x     | tlab.ini, grid, flow.*.?, scal.*.?| pdf*
|spectra.x  | tlab.ini, grid, flow.*.?, scal.*.?| xsp*, zsp*
|visuals.x  | tlab.ini, grid, flow.*.?, scal.*.?| *variable files* |

## List of examples

### 2D cases

* Case01. Shear layer with broadband ICs. Uniform grid. Kelvin-Helmholtz.  
* Case02. Same as Case01, but compressible
* Case03. Same as Case01, but with stretched grid.  
* Case04. Same as Case03, but compressible.  
* Case05. Same as Case03, but 2 scalars with different Schmidt * numbers.  
* Case06. Stably stratified density interface with discrete ICs: oscillating inversion.  
* Case07. Unstable  stratified density interface: Rayleigh-Taylor.  
* Case08. Stably stratified shear layer.  
* Case09. Heated plate.  
* Case10. Same as Case09, but two scalars with different BCs.  
* Case11. Half Channel flow.  
* Case12. Same as 11, implicit solver.  
* Case13. Channel flow.  
* Case14. Cloud-top mixing layer, airwater equilibrium compressible formulation.  
* Case15. Cloud-top mixing layer; airwater equilibrium incompressible formulation.  
* Case16. Cloud-top mixing layer; airwaterlinear, evaporation only. (Case06 plus buoyancy reversal.)  
* Case17. Cloud-top mixing layer; airwaterlinear, radiation only.  
* Case18. Cloud-top mixing layer; airwaterlinear, radiation case with evaporative cooling.  
* Case19. Cloud-top mixing layer; airwaterlinear, radiation case with evaporative cooling and settling.  
* Case20. Subsiding shell; airwater. Gravity vector along the horizontal.  
* Case21. Subsiding shell; airwater. Broadband perturbation.  
* Case22. Convective boundary layer with quadratic chemistry in three passive scalars.  
* Case23. Rayleigh-Benard convection with Dirichlet boundary conditions. 
* Case24. Airwater equilibrium incompressible formulation of stratocumulus-topped boundary layer.  
* Case25. Anelastic formulation of CBL.  
* Case26. Airvapor anelastic formulation of CBL.  
* Case27. Airwater equilibrium anelastic formulation of stratocumulus-topped boundary layer.   
* Case28. Same as Case27, but adding subsidence and sedimentation.  
* Case29. Same as Case28, but with dimensions.  
* Case30. Same as Case29, but with different radiation model.  

### 2D cases: Lagrangian routines

* Case31. Case01 (shear layer), saving only the particles.  
* Case32. Case01 (shear layer), saving trajectories as well (but fewer particles).  
* Case33. Case32, but inertia particles instead of tracers.  
* Case34. Case17 (cloud-top), solving liquid equation w/ & w/o diffusion.  
* Case35. Same as 34, but with a stratified bottom interface.  

### Gravity waves

* Case37. Wave maker.

### 3D cases

* Case41. Neutral Ekman layer without sponge at the top.  
* Case42. Neutral Ekman layer with sponge at the top.  
* Case44. Stable Ekman layer with sponge at the top.  
* Case45. Same as 21, implicit solver.  
* Case46. Same as 22, implicit solver.  
* Case47. Same as 23, implicit solver.  
* Case48. Neutral Ekman layer with sponge at the top and interactive BC at the bottom.  

### 1D cases

* Case50. 1D perturbed laminar Ekman layer, implicit solver.  

### 3D cases: Channel flow

* Case61. Channel flow with constant streamwise pressure gradient (Re_tau=180, rotation term for turbulent transition!).  
* Case62. Same as 61, with horizontal pressure staggering and a compact vertical pressure filter (without rotation term).  
* Case63. Same as 62, with IBM (streamwise aligned bars on lower and upper boundary).  

<!-- make checkrl/checkdb runs the check.sh bash-script inside each directory -->
