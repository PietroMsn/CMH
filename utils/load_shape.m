function [M] = load_shape(path)

    if exist([path,'.off'])
        M = load_off([path,'.off']);
    elseif exist([path,'.ply'])
        M = load_ply([path,'.ply']);
    elseif exist([path,'.obj'])
        M = load_obj([path,'.obj']);
        if size(M.TRIV,2) == 1
            clear M;
            [VERT,TRIV] = readOBJ([path,'.obj']);
            M.VERT = VERT; M.TRIV = TRIV;
            M.n = size(M.VERT,1); M.m = size(M.TRIV,1);
        end
        M.path = path;
    elseif exist([path,'.mat'])
        tmp = load([path,'.mat']);
        if isfield(tmp, 'shape')
           S = tmp.shape;
        elseif isfield(tmp, 'surface')
           S = tmp.surface;
        elseif  isfield(tmp, 'Src')

               S.X = tmp.Src.X(:,1);
               S.Y = tmp.Src.X(:,2);
               S.Z = tmp.Src.X(:,3);
               S.TRIV = tmp.Src.T;
           if isfield(tmp, 'Area')
              M.Area = tmp.Area; 
           end
           if isfield(tmp, 'landmarks1')
               M.landmark = tmp.landmarks1;       
           end
        else
            disp('ERROR: check the file .mat')
        end
        M.VERT = [S.X,S.Y,S.Z]; M.TRIV = S.TRIV;  
        M.n = size(S.X,1); M.m = size(S.TRIV,1);
        
        M.path = path;
    elseif exist([path,'.pcd'])
        M = pcread([path,'.pcd']);

    else
        disp('ERROR: no data available check path or name')
    end

    