function [gradPlantHat,uk_h,yk_h,xk_h] = CDAGradient(dxk,uk,parPlant,h,H)
% Estimate plant gradients via Central Difference approach (CDA)
%
% Inputs:
%    dxk = state values (guess)
%    uk = inputs (central value for the CDA)
%    parPlant = plant parameters
%    h = perturbation step
%    H = output map function (H*zk = yk)
%
% Outputs:
%    gradPlantHat - estimated gradient (measurements - defined by H)
%    uk_h = required input excitation
%    yk_h = measurements related to these inputs
%    xk_h = plant states related to these inputs

% Other m-files required: PlantModel
% Subfunctions: none
% MAT-files required: none

%perturbing the inputs
uh = uk + h;
u_h = uk - h;

%use the perturbed input to probe the plant
[~,~,~,~,xh,yh,~] = PlantModel(dxk,uh,parPlant);
[~,~,~,~,x_h,y_h,~] = PlantModel(dxk,u_h,parPlant);

%measure plant values at the probe points
xk_h = [H*xh, H*x_h];
yk_h = [H*yh, H*y_h];
uk_h = [uh,u_h];

%gradients estimates
gradPlantHat = (H*yh - H*y_h)/(2*h);

end

