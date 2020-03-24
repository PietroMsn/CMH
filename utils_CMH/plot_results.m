function plot_results(M1, M2, new2)

%     figure;
%     subplot(131)
%     trisurf(M1.TRIV, M1.VERT(:,1), M1.VERT(:,2), M1.VERT(:,3));
%     view([0, 90]); axis off; axis equal;
%     subplot(132)
%     trisurf(M2.TRIV, M2.VERT(:,1), M2.VERT(:,2), M2.VERT(:,3));
%     view([0, 90]); axis off; axis equal;
%     subplot(133)
%     trisurf(new2.TRIV, new2.VERT(:,1), new2.VERT(:,2), new2.VERT(:,3));
%     view([0, 90]); axis off; axis equal;
%     
    
    yrotation = -30;

    to_remesh.VERT = M1.VERT;
    to_remesh.TRIV = M1.TRIV;
    to_remesh.n = size(M1.VERT, 1);
    to_remesh.m = size(M1.TRIV, 1);
    to_remesh = mesh_with_consistent_fields(to_remesh);
    
    to_transfer.VERT = M2.VERT;
    to_transfer.TRIV = M2.TRIV;
    to_transfer.n = size(M2.VERT, 1);
    to_transfer.m = size(M2.TRIV, 1);
    to_transfer = mesh_with_consistent_fields(to_transfer);
    
    res.VERT = new2.VERT;
    res.TRIV = new2.TRIV;
    res.n = size(new2.VERT, 1);
    res.m = size(new2.TRIV, 1);
    res = mesh_with_consistent_fields(res);

    %h = figure;
    %subplot(131);
    figure;
    
    render_mesh_with_edges(to_remesh,...
        'RotationOps',{[-90,0,30],[0,0,0]},... % [-90,0,30] for others
        'CameraPos',[-40,3],...
        'FaceAlpha',1);
    title('Geometry');
    %subplot(132);
    figure;
    
    render_mesh_with_edges(to_transfer,...
        'RotationOps',{[-90,0,30],[0,0,0]},... % [-90,0,30] for others
        'CameraPos',[-40,3],...
        'FaceAlpha',1);
    title('Connectivity');
    %subplot(133);
    figure;
    
    render_mesh_with_edges(res,...
        'RotationOps',{[-90,0,30],[0,0,0]},... % [-90,0,30] for others
        'CameraPos',[-40,3],...
        'FaceAlpha',1);
    title('Geometry + connectivity');

    %set(h,'Position', get(0, 'Screensize')/0.7);
 

end