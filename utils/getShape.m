function [M, L] = getShape(path, land)

    shape = load_shape(path);
    
    M = shape; M.nv = size(M.VERT,1);
    M.path = path;
    M.A = calc_mass_matrix(M);
    M.lm = land;
    
    S.surface.X = M.VERT(:,1); S.surface.Y = M.VERT(:,2); S.surface.Z = M.VERT(:,3);
    S.surface.TRIV = M.TRIV;
    S.surface.nv = size(M.VERT,1);
    

    L = lb_basis_surface(S, 200);
    S = MESH('Src',M.VERT,M.TRIV);
   M.angview = S.angview; M.X = M.VERT; M.T = M.TRIV;
    
    M.Ae = L.A;
end