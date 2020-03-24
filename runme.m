% CMH: Coordinates Manifold Harmonics fo Functional Remeshing (3DOR 2019)
% CODE for remeshing using CMH pipeline
% AUTHORS:  Simone   Melzi , simone.melzi@univr.it
%           Riccardo Marin , riccardo.marin_01@univr.it
%           Pietro   Musoni , pietro.musoni@gmail.com
%           Filippo  Bardon , bardon.filippo@gmail.com
%           Marco    Tarini , marco.tarini@unimi.it
%           Umberto  Castellani , umberto.castellani@univr.it
% 
%              
% -------------------------------------------------------------------------
% The output of this code is the remeshed shape using CMH pipeline
% --------------------------------------+----------------------------------
% select input meshes as struct shapes, with the following fields:
% X:         matrix of dimension N x 3 that are the X,Y,Z coordinates of 
%            vertices in R^3
% T:         matrix of M x 3 matrix where each row is a triangle in
%            the mesh.

% landmarks: a vector 2x5 containing the indices of the
%            landmarks for the two input shapes
%
% In this code we make use of code from:
% "Functional Maps: A Flexible Representation of Maps Between Shapes" 
% "gptoolbox: Geometry Processing Toolbox"

clear all;  close all; 
addpath(genpath(pwd));

pathToRemesh = './shapes/TposeManBad';
pathConnectivity = './shapes/TposeWomanFat';

landmarks1 = [1581, 851, 741, 2362, 2782]'; 
landmarks2 = [5638, 2336, 1066, 8266, 8599]'; 

landmarks = [landmarks1 landmarks2];

timeTot = [];


%% parameters
visualize = 1; % put at 1 to visualize the plots
% the number of eigenfunctuions used (excluding CMH)
kM1 = 50; kM2 = 50; 
% 1 to use the relative descriptor, 0 to not use it
desc.WKS = 1; desc.HKS = 0; desc.LAND = 1; 
% subsampling of the descriptor function used. With more functions is more 
% accurate but slower
desc.skipSize = 10; 
hands_radius = 0.13;
a_arap=0.2;
delta_t=0.03;
%%

[Mesh1, L1] = getShape(pathToRemesh, landmarks1);
[Mesh2, L2] = getShape(pathConnectivity, landmarks2);

% scaling the two input meshes with the area of the SMPL model
[M1, M2] = scaling_shapes(Mesh1, Mesh2);

% computation of the CMH functions, the returning matrix is 
% #vertices x #eigenfunctions+3 
fprintf('Computing CMH for M1 and M2...\n');

timeCMH = tic;
M1_CMH = compute_CMH(M1,L1.evecs(:,1:kM1));
M2_CMH = compute_CMH(M2,L2.evecs(:,1:kM2));
timeTot = [timeTot, toc(timeCMH)];
fprintf('done in %f seconds\n', timeTot(end));

% compute the map using CMH, using the chosen descriptors 
fprintf('Estimating a functional map in the CMH basis...\n');
timeMap = tic;
[CMH_C, CMH_p2p] = compute_map_CMH(M1_CMH, M2_CMH, landmarks, kM1+3, ...
    kM2+3, L1, L2, M1, M2, desc);
timeTot = [timeTot, toc(timeMap)];
fprintf('done in %f seconds\n', timeTot(end));

% this plot shows the two input shapes with the 3 CMH functions
if visualize
    plotCMH(M1, M2, M1_CMH, M2_CMH, kM1, kM2);
end
% computes the map using MH (without CMH functions), using the chosen 
% descriptors 
[MH_C, MH_p2p] = compute_map_CMH(M1_CMH, M2_CMH, landmarks, kM1, kM2, ...
    L1, L2, M1, M2, desc);

if visualize
    plot_comparison_mh_cmh(M1, M2, MH_p2p, CMH_p2p);
end
% calculates the remeshing pipeline using the maps evaluated before 
%timeRem = tic;
[new, new2, timeTot] = perform_functional_remeshing(M1, M2, landmarks, ...
    CMH_p2p, hands_radius, a_arap, delta_t, timeTot);

% timeTot = [timeTot, toc(timeRem)];
% fprintf('      done %f \n', timeTot(end));

fprintf('Total execution time: %f seconds\n', sum(timeTot));
% plot of the two input shapes and the final result
plot_results(M1, M2, new2);
