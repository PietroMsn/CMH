function [new2, times] = arap_dato(new, M, N, delta_t, a_arap, times)

    new2.VERT=new.VERT; new2.TRIV=new.TRIV;


    [ distances, surface_points ] = point2trimesh('Faces',M.TRIV,'Vertices',M.VERT,...
        'QueryPoints',new2.VERT,'Algorithm','parallel_vectorized_subfunctions' );


    tarapdato = tic;
    for i=1:399
        if(mod(i,50)==0)
            disp(['iter:',num2str(i)]);
                          [ distances, surface_points ] = point2trimesh('Faces',M.TRIV,'Vertices',M.VERT,...
                                        'QueryPoints',new2.VERT,'Algorithm','parallel_vectorized_subfunctions' );
        end
        [G1,E] = arap_gradient(N.VERT,N.TRIV,new2.VERT);
        G2=new2.VERT-surface_points;


        new2.VERT=new2.VERT-delta_t*(a_arap*G1+G2);    
    end

    times = [times, toc(tarapdato)];
    %fprintf('      done %f \n', times(end));
    new2MESH=MESH('n2',new2.VERT,new2.TRIV);



end