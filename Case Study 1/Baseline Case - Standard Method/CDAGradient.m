function [gradPlantHat,uk_h,yk_h,xk_h] = CDAGradient(dxk,uk,parPlant,h,H)
% Estimate plant gradients via Forward Difference approach (FDA)
%
% Inputs:
%    yValuePlant = measured value of the initial point
%    dxk = initial differential states
%    uk = inputs
%    parPlant = plant parameters
%    h = method parameters
%    H = output map function (H*zk = yk)
%
% Outputs:
%    gradPlantHat - estimated gradient (measurements - defined by H)
%    uk_h = required input incitation
%    yk_h = measurements related to the inputs

% Other m-files required: PlantModel
% Subfunctions: none
% MAT-files required: none
%
% Author: Jose Otavio Matias
% email: jose.otavio@usp.br 
% March 2018; Last revision: 13-Mar-2018

%perturbation on the inputs
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

