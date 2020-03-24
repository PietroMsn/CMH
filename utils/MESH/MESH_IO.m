% Written by Omri Azencot and Danielle Ezuz
classdef MESH_IO
    
    properties
    end
    
    methods (Static)
        
        %%%%%%%%
        % Read %
        %%%%%%%%

        function [X,T] = roff(filename)
            
            fid = fopen(filename,'r');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end
            
            str = fgets(fid);   % -1 if eof
            
            if strcmp(str(1:4), 'COFF')
                [X,T,~] = readCoff(filename,4); % assume 4 color channels
                return;
            end
            
            if ~strcmp(str(1:3), 'OFF')
                error('The file is not a valid OFF one.');
            end
           
            str = fgets(fid);
            sizes = sscanf(str, '%d %d', 2);
            while length(sizes) ~= 2
                str = fgets(fid);
                sizes = sscanf(str, '%d %d', 2);
            end
            nv = sizes(1);
            nf = sizes(2);
            
            % Read vertices
            [X,cnt] = fscanf(fid,'%lf %lf %lf\n', [3,nv]);
            if cnt~=3*nv
                warning('Problem in reading vertices.');
            end
            X = X';
            
            [T,cnt] = fscanf(fid,'3 %ld %ld %ld\n', [3,inf]);
            T = T'+1;
            
            fclose(fid);
        end
        
        % generate texture coordinates per vertex by projecting to the plane of 
        % col1, col2 coordinates (signed) and stretch so that the texture appears 
        % mult_const times.
        function vt = generate_tex_coords(v, col1, col2, mult_const)

            vt = [sign(col1)*v(:, abs(col1)), sign(col2)*v(:, abs(col2))];
            vt = bsxfun(@minus, vt, min(vt));

            max_vt = max(vt(:));
            vt = mult_const * vt / max_vt;
        end
            
        % Output M to an obj file 
        % uv - texture coordinates.
        % iname - filename for texture image
        % fname - output file name
        function wobj(M, uv, iname, fname)

            options.object_texture = uv;
            options.nm_file = iname;

            [pathstr,name] = fileparts(fname); 
            write_obj2(pathstr, [name '.obj'], M.vertices, M.triangles, options);
        end

    end
end

