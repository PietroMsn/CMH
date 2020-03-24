function [p2p_n,matches,P1,P2] = LocalPat(M1,M2,landmarks,p2p,dist_thresh)


P1 = MESH.extract_patch(M1, calc_dist_map(M1, landmarks(1)) < dist_thresh);
P2 = MESH.extract_patch(M2, calc_dist_map(M2, landmarks(2)) < dist_thresh);
n_angles = 16;
angles = linspace(0, 2*pi*(1-1/n_angles), n_angles);
N_orig = P2;
solutions = cell(1,n_angles);
parfor iter=1:n_angles
    
    %fprintf('outer iter %d/%d\n', iter, n_angles)
    
    angle = angles(iter);
    R = [cos(angle) -sin(angle) 0 ; sin(angle) cos(angle) 0 ; 0 0 1];
    
    N_deformable = N_orig;
    N_deformable.VERT = N_orig.VERT*R;

    
    N_deformable = align_icp_cpd_hand(P1, N_deformable, 3);
    matches = knnsearch(P1.X, N_deformable.VERT);
    
    solutions{iter}.N = N_deformable;
    solutions{iter}.matches = matches;
   
end
best_sol = 0;
best_distortion = Inf;
for i=1:length(solutions)
    
    matches = solutions{i}.matches;
    
    N_pts = fps_euclidean(N_orig.VERT, 100, 1);
    M_pts = matches(N_pts);
    
    DM = calc_dist_matrix(P1, M_pts);
    DN = calc_dist_matrix(N_orig, N_pts);
    distortion = sum(sum(abs(DM - DN)));
    
    if distortion<best_distortion
        best_distortion = distortion;
        best_sol = i;
    end
    %fprintf('(%d) Distortion: %.4f\n', i, distortion)
end

%fprintf('Best solution: %.4f\n', best_distortion)
matches = solutions{best_sol}.matches;
Dato_idx = N_orig.fullshape_idx;
smpl_idx = P1.fullshape_idx(matches);

p2p_n=p2p;
p2p_n(Dato_idx)=smpl_idx;



end