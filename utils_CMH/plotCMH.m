function plotCMH(M1, M2, M1_CMH, M2_CMH, kM1, kM2)


    figure;
    subplot(231); colormap(bluewhitered);
    trisurf(M1.TRIV, M1.VERT(:,1), M1.VERT(:,2), M1.VERT(:,3), M1_CMH(:,kM1+1),'SpecularStrength',0.15); 
    view([0, 90]); axis equal; axis off, view([0 90]); light; lighting phong; shading interp;
    subplot(232);
    trisurf(M1.TRIV, M1.VERT(:,1), M1.VERT(:,2), M1.VERT(:,3), M1_CMH(:,kM1+2),'SpecularStrength',0.15); 
    view([0, 90]); axis off; axis equal; light; lighting phong; shading interp;
    subplot(233);
    trisurf(M1.TRIV, M1.VERT(:,1), M1.VERT(:,2), M1.VERT(:,3), M1_CMH(:,kM1+3),'SpecularStrength',0.15); 
    view([0, 90]); axis off; axis equal; light; lighting phong; shading interp;
    subplot(234);
    trisurf(M2.TRIV, M2.VERT(:,1), M2.VERT(:,2), M2.VERT(:,3), M2_CMH(:,kM2+1),'SpecularStrength',0.15); 
    view([0, 90]); axis off; axis equal; light; lighting phong; shading interp;
    subplot(235);
    trisurf(M2.TRIV, M2.VERT(:,1), M2.VERT(:,2), M2.VERT(:,3), M2_CMH(:,kM2+2),'SpecularStrength',0.15); 
    view([0, 90]); axis off; axis equal; light; lighting phong; shading interp;
    subplot(236);
    trisurf(M2.TRIV, M2.VERT(:,1), M2.VERT(:,2), M2.VERT(:,3), M2_CMH(:,kM2+3),'SpecularStrength',0.15); 
    view([0, 90]); axis off; axis equal; light; lighting phong; shading interp;

end