function plot_comparison_mh_cmh(M1, M2, MH_p2p, CMH_p2p)

    mhRes.VERT = M1.VERT(MH_p2p,:);
    mhRes.TRIV = M2.TRIV;
    
    cmhRes.VERT = M1.VERT(CMH_p2p,:);
    cmhRes.TRIV = M2.TRIV;


    geom.VERT = M1.VERT;
    geom.TRIV = M1.TRIV;
    geom.n = size(M1.VERT, 1);
    geom.m = size(M1.TRIV, 1);
    geom = mesh_with_consistent_fields(geom);
    
    mhRes.VERT = M1.VERT(MH_p2p,:);
    mhRes.TRIV = M2.TRIV;
    mhRes.n = size(M1.VERT, 1);
    mhRes.m = size(M2.TRIV, 1);
    mhRes = mesh_with_consistent_fields(mhRes);
    
    cmhRes.VERT = M1.VERT(CMH_p2p,:);
    cmhRes.TRIV = M2.TRIV;
    cmhRes.n = size(M1.VERT, 1);
    cmhRes.m = size(M2.TRIV, 1);
    cmhRes = mesh_with_consistent_fields(cmhRes);
    
    figure;
    subplot(131)
    render_mesh_with_edges(geom,...
        'RotationOps',{[-90,0,30],[0,0,0]},... % [-90,0,30] for others
        'CameraPos',[-40,3],...
        'FaceAlpha',1);
    title('Geometry');
    
    subplot(132)
    render_mesh_with_edges(mhRes,...
        'RotationOps',{[-90,0,30],[0,0,0]},... % [-90,0,30] for others
        'CameraPos',[-40,3],...
        'FaceAlpha',1);
    title('MH reconstruction k=50');
    
    subplot(133)
    render_mesh_with_edges(cmhRes,...
        'RotationOps',{[-90,0,30],[0,0,0]},... % [-90,0,30] for others
        'CameraPos',[-40,3],...
        'FaceAlpha',1);
    title('CMH reconstruction k=53');
    

end