dist_thresh = 0.20;
P1 = MESH.extract_patch(M1, calc_dist_map(M1, landmarks(2,1)) < dist_thresh);
P2 = MESH.extract_patch(M2, calc_dist_map(M2, landmarks(2,2)) < dist_thresh);

n_angles = 16;
angles = linspace(0, 2*pi*(1-1/n_angles), n_angles);

N_orig = P2;
solutions = cell(1,n_angles);

parfor iter=1:n_angles
    
    fprintf('outer iter %d/%d\n', iter, n_angles)
    
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
    
    fprintf('(%d) Distortion: %.4f\n', i, distortion)
    
end

fprintf('Best solution: %.4f\n', best_distortion)
matches = solutions{best_sol}.matches;
Dato_idx = N_orig.fullshape_idx;
smpl_idx = P1.fullshape_idx(matches);

p2p_n=p2p;
p2p_n(Dato_idx)=smpl_idx;

MESH.plot_dense_matches(M1,M2,p2p_n);


P1_B=MESH.func_basis(P1,120);
P2_B=MESH.func_basis(P2,120);
Pi=sparse([1:P2.nv], matches,1,P2.nv,P1.nv);
C2= P2_B'*P2.Ae*Pi*P1_B;


[Cm2, p2pm, w_m] = fmap_to_precise(P2, P1, C2, P2_B, P1_B);
 new = full(Cm2 * P1.VERT);
new_mirela=new;
 
P3=MESH('New',new,P2.TRIV);
a_arap=0.5;
delta_t=0.01;

for i=1:200
[G1,E] = arap_gradient(P2.X,P2.T,new);
new=new-delta_t*(a_arap*G1);
subplot(121);
trisurf(P3.TRIV,P3.X(:,1),P3.X(:,2),P3.X(:,3),zeros(P3.nv,1),'SpecularStrength',0.15);
axis equal; axis off, view([0 90]); light; lighting phong; 
subplot(122);
trisurf(P3.TRIV,new(:,1),new(:,2),new(:,3),zeros(P3.nv,1),'SpecularStrength',0.15);
axis equal; axis off, view([0 90]); light; lighting phong; 
pause(0.01);
end
P_Final_Hand=MESH('pp',new,P2.TRIV);
new_M=new;


P3=MESH('New',P1.VERT(matches,:),P2.TRIV);
 new = P1.VERT(matches,:);

a_arap=0.5;
delta_t=0.01;

for i=1:200
[G1,E] = arap_gradient(P2.X,P2.T,new);
new=new-delta_t*(a_arap*G1);
subplot(121);
trisurf(P3.TRIV,P3.X(:,1),P3.X(:,2),P3.X(:,3),zeros(P3.nv,1),'SpecularStrength',0.15);
axis equal; axis off, view([0 90]); light; lighting phong; 
subplot(122);
trisurf(P3.TRIV,new(:,1),new(:,2),new(:,3),zeros(P3.nv,1),'SpecularStrength',0.15);
axis equal; axis off, view([0 90]); light; lighting phong; 
pause(0.01);
end
P_Final_Hand_CPD=MESH('ppc',new,P2.TRIV);


 
 
 subplot(321);
trisurf(P1.TRIV,P1.X(:,1),P1.X(:,2),P1.X(:,3),zeros(P1.nv,1),'SpecularStrength',0.15);
axis equal; axis off, view([0 90]); light; lighting phong; 
subplot(322);
trisurf(P2.TRIV,P2.X(:,1),P2.X(:,2),P2.X(:,3),zeros(P2.nv,1),'SpecularStrength',0.15);
axis equal; axis off, view([0 90]); light; lighting phong; 
subplot(323);
trisurf(P2.TRIV,P1.VERT(matches,1),P1.VERT(matches,2),P1.VERT(matches,3),zeros(P2.nv,1),'SpecularStrength',0.15);
axis equal; axis off, view([0 90]); light; lighting phong; 
title('CPD')

subplot(324);
trisurf(P2.TRIV,new_mirela(:,1),new_mirela(:,2),new_mirela(:,3),zeros(P2.nv,1),'SpecularStrength',0.15);
axis equal; axis off, view([0 90]); light; lighting phong; 
title('CPD+Mirela')

subplot(325);
trisurf(P2.TRIV,P_Final_Hand_CPD.X(:,1),new(:,2),new(:,3),zeros(P2.nv,1),'SpecularStrength',0.15);
axis equal; axis off, view([0 90]); light; lighting phong; 
title('(CPD)+ARAP')

subplot(326);
trisurf(P2.TRIV,new_M(:,1),new_M(:,2),new_M(:,3),zeros(P2.nv,1),'SpecularStrength',0.15);
axis equal; axis off, view([0 90]); light; lighting phong; 
title('(CPD+Mirela)+ARAP')


X1= M3.X;
L = X1(P2.fullshape_idx,:);
X1(P2.fullshape_idx,:)=P_Final_Hand_CPD.X;

trisurf(M3.T,X1(:,1),X1(:,2),X1(:,3),zeros(Tar.nv,1),'SpecularStrength',0.15);
axis equal; axis off, view([0 90]); light; lighting phong; 


