% Written by Omri Azencot 
classdef MESH < handle
    
    properties (Access='public')
        
        name
        
        vertices, X, VERT
        triangles, T, TRIV
        area, Ae
        
        nv,n              % #vertices
        nf,m              % #faces
        
        va              % vertex areas
        vta             % vertex total area
        ta              % face areas
        percva
        
        VA, IVA         % (lumped) vertices mass matrix and inverse
        TA, ITA         % (lumped) triangles mass matrix and inverse
        
        Nf              % normal-per-face
                
        G               % gradient of vertex functions
        D               % divergence of face fields
        
        E1, E2, E3      % triangle edges
        
        mel             % mean edge length
        
        vdeg            % vertex degrees
        rv              % regular vertices
        E2V, T2E, E2T, T2T
        AR              % Face Quality
        degtri          % Degenerate Triangles
        angview=[0 90]
        fullshape_idx=[]
    end
    
    
    methods
        function [ mesh ] = MESH( meshname, v, t)
            if nargin < 1; meshname = 'sphere_s3'; end
            mesh.name = meshname;

            if nargin < 3   
                [v, t] = MESH_IO.roff([meshname '.off']);
            end
            
            mesh.vertices = v; mesh.X=full(mesh.vertices); mesh.VERT=v;
            mesh.triangles = double(t); mesh.T=full(mesh.triangles); mesh.TRIV=t;
            
            % mesh data structures
                        
            mesh.nv = size(mesh.vertices,1);mesh.n=mesh.nv;
            mesh.nf = size(mesh.triangles,1);mesh.m=mesh.nf;
            
            
            % simplices' areas and mass matrices
            mesh.ta = face_areas( mesh );
            mesh.va = vertex_areas( mesh );
            mesh.vta = sum(mesh.va);
            mesh.percva=mesh.va./mesh.vta;
            [mesh.VA,mesh.IVA,mesh.TA,mesh.ITA] = mass_matrices( mesh );
            mesh.area=mesh.VA; mesh.Ae=mesh.VA;
            [mesh.E2V, mesh.T2E, mesh.E2T, mesh.T2T] = conn(mesh);
            
            a = unique(mesh.E2V);
            out = [a,histc(mesh.E2V(:),a)];
            mesh.vdeg = out(:,2);
            rv = sum(mesh.vdeg==6);
            % simplices' normals
            mesh.Nf = face_normals( mesh );
            
            % geometric and differential operators
            [mesh.E1,mesh.E2,mesh.E3] = face_edges( mesh );
                        
            mesh.G = grad( mesh );
            mesh.D = div( mesh );
            
            mesh.AR= quality(mesh);
            mesh.degtri= degenerate_triangles(mesh);
        end
        
        function AR = quality(mesh)
            a = MESH.normv(mesh.E1); b = MESH.normv(mesh.E2); c=MESH.normv(mesh.E3);
            S = (a + b + c)/2;
            AR = (a.*b.*c) ./ (8 .* (S-a).*(S-b) .* (S-c));
            
        end
        
        function [E2V, T2E, E2T, T2T] = conn(mesh) 
           [E2V, T2E, E2T, T2T]=connectivity(mesh.triangles);
        end
        
        function [ c] = degenerate_triangles( mesh )
        %DEGENERATE TRIANGLES
                c=0;
                edge = [MESH.normv(mesh.E1), MESH.normv(mesh.E2), MESH.normv(mesh.E3)];
                edge = sort(edge,2);

                % check if triangle is degenerate
                flag = (edge(:,1) + edge(:,2)) <= edge(:,3);
                c = nnz(flag);
        end
        
        function [ ta ] = face_areas( mesh )
            X = mesh.vertices;
            T = mesh.triangles;
            
            P1 = X(T(:,1),:) - X(T(:,2),:);
            P2 = X(T(:,1),:) - X(T(:,3),:);
            
            ta = mesh.normv( cross( P1, P2 ) ) / 2;
        end
        
        function [ va ] = vertex_areas( mesh )
            va = full( sum( mass_matrix(mesh), 2 ));
            
        end
        
        function [ M ] = mass_matrix( mesh )
            T = double( mesh.triangles ); 
            
            I = [T(:,1);T(:,2);T(:,3)];
            J = [T(:,2);T(:,3);T(:,1)];
            Mij = 1/12*[mesh.ta; mesh.ta; mesh.ta];
            Mji = Mij;
            Mii = 1/6*[mesh.ta; mesh.ta; mesh.ta];
            In = [I;J;I];
            Jn = [J;I;I];
            Mn = [Mij;Mji;Mii];
            M = sparse(In,Jn,Mn,mesh.nv,mesh.nv);            
        end

        function [ VA, IVA, TA, ITA] = mass_matrices( mesh )
            sv = mesh.nv; sf = mesh.nf; 
            VA = spdiags(mesh.va,0,sv,sv);
            IVA = spdiags(1./mesh.va,0,sv,sv);
            TA = spdiags([mesh.ta; mesh.ta; mesh.ta],0,3*sf,3*sf);
            ITA = spdiags(1./[mesh.ta; mesh.ta; mesh.ta],0,3*sf,3*sf);
        end
        
        function [ Nf ] = face_normals( mesh )          
            X = mesh.vertices;
            T = mesh.triangles;
            
            P1 = X(T(:,1),:) - X(T(:,2),:);
            P2 = X(T(:,1),:) - X(T(:,3),:);
            
            Nf = cross( P1, P2 );
            Nf = MESH.normalize_vf( Nf );
        end
        
        function [ E1, E2, E3 ] = face_edges( mesh )
            X = mesh.vertices;
            T = mesh.triangles;
            
            E1 = X(T(:,3),:) - X(T(:,2),:);
            E2 = X(T(:,1),:) - X(T(:,3),:);
            E3 = X(T(:,2),:) - X(T(:,1),:);
            
            
            
            E = [E1; E2; E3];
            
            mesh.mel = mean( MESH.normv( E ) );
        end
        
        
        function [ G ] = grad( mesh )
            % G corresponds to eq. (3.9) in Polygon mesh processing book
            I = repmat(1:mesh.nf,3,1);
            II = [I(:); I(:)+mesh.nf; I(:)+2*mesh.nf];

            J = double( mesh.triangles' );
            JJ = [J(:); J(:); J(:)];

            RE1 = mesh.rotate_vf( mesh.E1 );
            RE2 = mesh.rotate_vf( mesh.E2 );
            RE3 = mesh.rotate_vf( mesh.E3 );

            S = [RE1(:) RE2(:) RE3(:)]'; SS = S(:);
            
            G = sparse(II,JJ,SS,3*mesh.nf,mesh.nv);
            G = .5 * mesh.ITA * G;
        end
        
        function [ D ] = div( mesh )
            % D corresponds to eq. (3.12) in Polygon mesh processing book
            D = - mesh.IVA * mesh.G' * mesh.TA;
        end

        function [ rvf ] = rotate_vf( mesh, vf )
            vf = reshape(vf,mesh.nf,3);
            rvf = cross( mesh.Nf, vf );
        end
        function plot(mesh,f)
             if not(exist('f'))
                 f=mesh.vertices(:,1);
             end
            trisurf(mesh.triangles,mesh.vertices(:,1),mesh.vertices(:,2),mesh.vertices(:,3),f,'EdgeColor','None');
            axis off equal; light; lighting phong;
            view(mesh.angview);
        end
       


    end

         
    methods (Static)
        
        function [ nf ] = normalize_f( f )
            nf = (f - min(f)) / (max(f)-min(f));
        end
        
        function [ snv ] = snormv( v )
            snv = sum(v.^2,2);
        end
        
        function [ nv ] = normv( v )
            nv = sqrt(sum(v.^2,2));
        end
        
        function [ nnv ] = normalize_vf( v )
            vn = MESH.normv(v); I = vn > 0;
            nnv = v;
            nnv(I,:) = v(I,:) ./ repmat(vn(I),1,3);
        end

        function [ B, BI, D ] = func_basis( M, k )
            L = - M.D * M.G;
            A = M.VA;
            W = A*L;
            tic; fprintf('\tEIGS of Laplace--Beltrami operator: ');
            [ev,el] = eigs(W,A,k,'SM');
            fprintf('%f sec\n',toc);
            
            [~,ii] = sort(diag(el)); el = el(ii,ii); ev = ev(:,ii);
            B = ev; BI = ev'*A; D = el;
        end
        
        function [ d ] = distance_p2s( points,M )
            d = p2s_dist(points', M.vertices', M.triangles', knnsearch(M.vertices, points)');
        end
         
        function [ e, metro ] = Metro( M1,M2 )
            S_1=M1.vertices;
            e = MESH.distance_p2s(S_1,M2);
            metro = 1/size(S_1,1) * sum (e);
        end   
        
        function A = calc_adj_matrix(M, sym)

            % asymmetric
            A = sparse(...
                [M.T(:,1); M.T(:,2); M.T(:,3)], ...
                [M.T(:,2); M.T(:,3); M.T(:,1)], ...
                ones(3 * M.nf, 1), ...
                M.nv, M.nv, 3 * M.nf);

            if nargin==2 && sym
                A = A+A';
                A = double(A~=0);
            end

        end


        
       function [lut edges] = Adj2Lut(adjacencyMtx)
                nVertices = length(adjacencyMtx);

                [r c] = find(tril(adjacencyMtx)); %<- tril = assuming non-directed graph
                nEdges = length(r);

                lut = cell(nVertices,1);

                % % pre-allocate?
                % nNeighbors   = accumarray([r;c],ones(1,nEdges*2));
                % maxNeighbors = max(nNeighbors);

                for iEdge = 1:nEdges
                    currR = r(iEdge);
                    currC = c(iEdge);
                    lut{currR} = [lut{currR} currC];
                    lut{currC} = [lut{currC} currR];
                end

                if nargout > 1
                    edges = [r c];
                end

        end
            
            
       function onering = calc_onering(M)
        % find 1-ring neighbors for each shape vertex

        if isfield(M, 'adj')
            adj = M.adj;
        else
            adj = MESH.calc_adj_matrix(M);
        end

        onering = MESH.Adj2Lut(adj + adj');


        end


        function mask = grow(M , seed, n_edges)
            if n_edges<0
                error('Cannot grow a negative number of steps.')
            end

            onering = MESH.calc_onering(M);

            dilated = false(M.nv,n_edges);
            mask = seed;

            if n_edges==0
                return;
            end

            for k=1:n_edges
                for i=1:M.nv
                    if ~mask(i) && ~any(dilated(i,1:k-1))
                        if any(mask(onering{i})) || any(any(dilated(onering{i},1:k-1),2))
                            dilated(i,k) = 1;
                        end
                    end
                end
                mask = mask | dilated(:,k);
            end
        end
        
        function P = extract_patch(M, tris_to_keep)
            if size(tris_to_keep,1)==M.nv
                tris_to_keep = tris_to_keep(M.T(:,1)) | tris_to_keep(M.T(:,2)) | tris_to_keep(M.T(:,3));
            end

            P = {};

            P.T = M.T(tris_to_keep,:);
            P.nf = size(P.T,1);

            verts_to_keep = unique(P.T);
            P.V = M.X(verts_to_keep,:);
            P.nv = size(P.V,1);

            P.fullshape_idx = verts_to_keep;

            old_to_new = zeros(M.nv,1);
            old_to_new(verts_to_keep) = 1:length(verts_to_keep);

            P.T = reshape(old_to_new(P.T(:)), P.nf, 3);
            
            P = MESH('patch',P.V,P.T);
            P.fullshape_idx = verts_to_keep;
            P.angview= M.angview;
        end


        function bd = calc_boundary_edges(M1)

            triangles = M1.T;
            [c,d,~,~] = MESH.get_boundary(triangles);
            bd = zeros(length(d),2);

            for i=1:length(c)
                t = triangles(c(i),:);
                v = true(1,3);
                v(d(i)) = false;
                v = t(v);
                bd(i,1) = v(1);
                bd(i,2) = v(2);
            end
        end
        
function [ c,d,I,v ] = get_boundary( tri )
        %GET_BOUNDARY determines the boundary edges of a triangular mesh 
        %   [c,d] = get_boundary(tri) takes as input a list tri of consistently oriented triangles
        %   returns the indices c of the triangles the boundary edges belong to and
        %   the (local) indices d (in {1,2,3}) of the vertices opposing the boundary
        %   edge
        %   One gets the global indices of those vertices via F(sub2ind(size(F),c,d))
        %   Via 
        %   d1 = mod(d+1,3); d1(d1==0) = 3;
        %   d2 = mod(d+2,3); d2(d2==0) = 3;
        %   one gets the local indices of the boundary vertices.


        % Transpose matrix if neccessary
        if size(tri,1)<size(tri,2)
            tri=tri';
        end
        m = size(tri,1);

        %% Check for opposing halfedges

        % Matrix of directed edges
        I = [tri(:,1) tri(:,2);
             tri(:,2) tri(:,3);
             tri(:,3) tri(:,1)];

        b = not(ismember(I(:,[2 1]),I,'rows'));
        b = find(b);



        % Triangle indices
        c = mod(b,m);
        c(c==0) = m;

        % vertex opposing boundary edge
        d = floor((b-1)/m);
        d(d==0)=3;


        % % Directed boundary edges
        I=I(b,:);

        % Boundary vertices
        [~,~,v] = find(I);
        v = unique(v);

end
function  plot_dense_matches(M1, M2, matches)

        num_sub=size(matches,2);

        % M is the full shape
        % N is the partial shape
        M1.n = size(M1.X,1);
        M2.n = size(M2.X,1);

        colors = create_colormap(M1,M1);
        figure,
        subplot(1,num_sub+1,1);
        colormap(colors),plot_scalar_map(M1, 1:M1.n),axis off
        shading flat
        for i=1:num_sub
        freezeColors
        %N2 = M1;
        valid = matches(:,i)>0;
        N2.X = M1.X(matches(valid,i),:);
        % N2.VERT(~valid,:) = N2.VERT(~valid,:).*0;
        colors=[];
        colors(valid,:) = create_colormap(N2,M1);
        colors(~valid,:) = colors(~valid,:).*nan;
        subplot(1,num_sub+1,i+1);
        colormap(colors),plot_scalar_map(M2, 1:M2.n),axis off
        shading flat
        end
end
    end


        
end