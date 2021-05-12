function [gradPlantHat,uk_h1,uk_h2,yk_h1,yk_h2] = CDAGradient(yk,xk,uk,fk,parPlant,h,H,count)
% Estimate plant gradients via Central Difference approach (CDA)
%
% Inputs:
%    yValuePlant = current measurements
%    xk = current plant state
%    uk = current inputs
%    fk = current feed
%    parPlant = plant parameters
%    h = method parameters
%    H = output map function (H*zk = yk)
%
% Outputs:
%    gradPlantHat - estimated gradient (measurements - defined by H)
%    uk_h1,uk_h2 = required input incitation'
%    yk_h1,yk_h2 = measurements related to the inputs

% Other m-files required: PlantModel
% Subfunctions: none
% MAT-files required: none
%
% Author: Jose Otavio Matias
% email: jose.otavio@usp.br 
% October 2017; Last revision: 11-Oct-2019

%perturbation on the inputs
uk_h1 = uk + [1;0]*h;
uk_h2 = uk + [0;1]*h;
uk_h_1 = uk - [1;0]*h;
uk_h_2 = uk - [0;1]*h;

%use the perturbed input to probe the plant
[~,xm_h1,~,~] = PlantModel(xk,uk_h1,fk,parPlant,count);
[~,xm_h2,~,~] = PlantModel(xk,uk_h2,fk,parPlant,count);
[~,xm_h_1,~,~] = PlantModel(xk,uk_h_1,fk,parPlant,count);
[~,xm_h_2,~,~] = PlantModel(xk,uk_h_2,fk,parPlant,count);

%measure plant values at the probe points
yk_h1 = H*xm_h1;
yk_h2 = H*xm_h2;
yk_h_1 = H*xm_h_1;
yk_h_2 = H*xm_h_2;

%gradients estimates
gradPlantHat = [(yk_h1 - yk_h_1)/(2*h),(yk_h2 - yk_h_2)/(2*h)];

end

