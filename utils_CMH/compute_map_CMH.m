function [F_lb2, pF_lb2] = compute_map_CMH(lmh1,lmh2, landmarks,numEigsSrc,numEigsTar, L1, L2, DATO, MAGLIA, desc, D1, D2)

Phi1 = lmh1(:,1:numEigsSrc); 
Phi2 = lmh2(:,1:numEigsTar); 


%% Descriptors
fct_src = [];
fct_tar = [];
%fprintf('     Computing the descriptors...\n');%tic;

%compute the descriptors based on Wave Kernel Signature
if desc.WKS
    fct_src = [fct_src, waveKernelSignature(L1.evecs, L1.evals, L1.A, 200)];
    fct_tar = [fct_tar, waveKernelSignature(L2.evecs, L2.evals, L2.A, 200)];
end

% compute the descriptors using Wave Kernel Map around the landmarks
if desc.LAND
    fct_src = [fct_src, waveKernelMap(L1.evecs, L1.evals, L1.A, 200, landmarks(:,1))];
    fct_tar = [fct_tar, waveKernelMap(L2.evecs, L2.evals, L2.A, 200, landmarks(:,2))];
end

% compute the descriptors based on Heat Kernel Signature
if desc.HKS
    nt = 100;        % number of time steps
    min_t = 0.03;   % minimum time step
    max_t = 0.50;   % maximum time step
    
    log_ts = linspace(log(min_t), log(max_t), nt);
    ts = exp(log_ts);
    fct_src = [fct_src, log(hks(L1, ts))];
    fct_tar = [fct_tar, log(hks(L2, ts))];
end


% subsampling of the descriptor functions, too functions could slow down
% the process
fct_src = fct_src(:,1:desc.skipSize:end);
fct_tar = fct_tar(:,1:desc.skipSize:end);




if exist('D1')
    fct_src = [fct_src, waveKernelMapR(L1.evecs, L1.evals, L1.A, 20, D1(:,2:end))];
    fct_tar = [fct_tar, waveKernelMapR(L2.evecs, L2.evals, L2.A, 20, D2(:,2:end))];
end


%fprintf('done computing descriptors (%d on source and %d on target)\n',size(fct_src,2),size(fct_tar,2)); %toc;

assert(size(fct_src,2)==size(fct_tar,2));

% Normalization
no = sqrt(diag(fct_src'*L1.A*fct_src))';
fct_src = fct_src ./ repmat(no, [size(L1.evecs,1),1]);
fct_tar = fct_tar ./ repmat(no, [size(L2.evecs,1),1]);

%fprintf('Pre-computing the multiplication operators...');%tic;
%% Multiplication Operators
numFct = size(fct_src,2);
OpSrc = cell(numFct,1);
OpTar = cell(numFct,1);
for i = 1:numFct
    OpSrc{i} = Phi1'*L1.A*(repmat(fct_src(:,i), [1,numEigsSrc]).*Phi1);
    OpTar{i} = Phi2'*L2.A*(repmat(fct_tar(:,i), [1,numEigsTar]).*Phi2);
end

Fct_src = Phi1'*L1.A*fct_src;
Fct_tar = Phi2'*L2.A*fct_tar;
%fprintf('done\n');%toc;

%% Fmap Computation
%fprintf('Optimizing the functional map...\n');%tic;
Dlb = (repmat(L1.evals(1:numEigsSrc), [1,numEigsTar]) - repmat(L2.evals(1:numEigsTar)', [numEigsSrc,1])).^2;
Dlb = Dlb/norm(Dlb, 'fro')^2;
constFct = sign(Phi1(1,1)*Phi2(1,1))*[sqrt(sum(sum(L2.A))/sum(sum(L1.A))); zeros(numEigsTar-1,1)];

a = 1e-1; % Descriptors preservation
b = 1;    % Commutativity with descriptors
c = 1e-3; % Commutativity with Laplacian 
funObj = @(F) deal( a*sum(sum((reshape(F, [numEigsTar,numEigsSrc])*Fct_src - Fct_tar).^2))/2 + b*sum(cell2mat(cellfun(@(X,Y) sum(sum((X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y).^2)), OpTar', OpSrc', 'UniformOutput', false)), 2)/2 + c*sum( (F.^2 .* Dlb(:))/2 ),...
            a*vec((reshape(F, [numEigsTar,numEigsSrc])*Fct_src - Fct_tar)*Fct_src') + b*sum(cell2mat(cellfun(@(X,Y) vec(X'*(X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y) - (X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y)*Y'), OpTar', OpSrc', 'UniformOutput', false)), 2) + c*F.*Dlb(:));
funProj = @(F) [constFct; F(numEigsTar+1:end)];

F_lb = zeros(numEigsTar*numEigsSrc, 1); F_lb(1) = constFct(1);

% Compute the optional functional map using a quasi-Newton method.
options.verbose = 0;
F_lb = reshape(minConf_PQN(funObj, F_lb, funProj, options), [numEigsTar,numEigsSrc]);
%fprintf('done fmap optimization.\n');%toc;


%%
%fprintf('ICP refinement...');%tic;
[F_lb2, ~] = icp_refine(Phi1, Phi2, F_lb, 5);
%fprintf('done\n');%toc;

%% Evaluation
% Compute the p2p map

%fprintf('Mapping the vertices of M1 on M2...');%tic;
% fmap after ICP 
pF_lb2 = knnsearch((F_lb2*Phi1')', Phi2);
%fprintf('done\n');%toc;

