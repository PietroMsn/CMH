% Runs the MATLAB CMH pipeline on the example shapes (CPD module off) and
% saves the intermediate results in tests/reference/reference_cmh.mat.
% The Python tests (tests/test_against_matlab.py) compare against them.
%
% Run from any folder with:
%   matlab -batch "run('<path-to-repo>/matlab/export_reference.m')"

clear all; close all;
root = fileparts(fileparts(mfilename('fullpath')));
cd(root);
addpath(genpath(root));

landmarks1 = [1581, 851, 741, 2362, 2782]';
landmarks2 = [5638, 2336, 1066, 8266, 8599]';
landmarks = [landmarks1 landmarks2];

kM1 = 50; kM2 = 50;
desc.WKS = 1; desc.HKS = 0; desc.LAND = 1;
desc.skipSize = 10;
hands_radius = 0.13;
a_arap = 0.2;
delta_t = 0.03;

[Mesh1, L1] = getShape('./shapes/TposeManBad', landmarks1);
[Mesh2, L2] = getShape('./shapes/TposeWomanFat', landmarks2);
[M1, M2] = scaling_shapes(Mesh1, Mesh2);

M1_CMH = compute_CMH(M1, L1.evecs(:,1:kM1));
M2_CMH = compute_CMH(M2, L2.evecs(:,1:kM2));

[CMH_C, CMH_p2p, CMH_C0] = compute_map_CMH(M1_CMH, M2_CMH, landmarks, ...
    kM1+3, kM2+3, L1, L2, M1, M2, desc);
[MH_C, MH_p2p, MH_C0] = compute_map_CMH(M1_CMH, M2_CMH, landmarks, ...
    kM1, kM2, L1, L2, M1, M2, desc);

% a single ARAP gradient evaluation on the initial reconstruction
arap_U = M1.VERT(CMH_p2p, :);
arap_G = arap_gradient(M2.VERT, M2.TRIV, arap_U);

[new, new2] = perform_functional_remeshing(M1, M2, landmarks, CMH_p2p, ...
    hands_radius, a_arap, delta_t, [], 0);

M1_VERT = M1.VERT; M1_TRIV = M1.TRIV; M1_A = full(diag(M1.A));
M2_VERT = M2.VERT; M2_TRIV = M2.TRIV; M2_A = full(diag(M2.A));
L1_evecs = L1.evecs; L1_evals = L1.evals; L1_A = full(diag(L1.A));
L2_evecs = L2.evecs; L2_evals = L2.evals; L2_A = full(diag(L2.A));
new_VERT = new.VERT;
new2_VERT = new2.VERT;

if ~exist('./tests/reference', 'dir')
    mkdir('./tests/reference');
end
save('./tests/reference/reference_cmh.mat', ...
    'M1_VERT', 'M1_TRIV', 'M1_A', 'M2_VERT', 'M2_TRIV', 'M2_A', ...
    'L1_evecs', 'L1_evals', 'L1_A', 'L2_evecs', 'L2_evals', 'L2_A', ...
    'M1_CMH', 'M2_CMH', 'CMH_C', 'CMH_C0', 'CMH_p2p', 'MH_C', 'MH_C0', ...
    'MH_p2p', 'arap_U', 'arap_G', 'new_VERT', 'new2_VERT', '-v7');
fprintf('Reference saved in tests/reference/reference_cmh.mat\n');
