function [new, new2, times] = perform_functional_remeshing(M1, M2, landmarks, p2p, hands_radius, a_arap, delta_t, times)

    S1 = MESH('Src',M1.VERT,M1.TRIV);
    S2 = MESH('Tar',M2.VERT,M2.TRIV);
    
    
    % run CPD on the local patches around the landmarks relative to the
    % hands
    fprintf('Applying Coherent Point Drift on local patches...\n')
    tCPD = tic;
    [p2pCPD,matches,P1,P2] = LocalPat(M1,M2,[landmarks(2,1),landmarks(2,2)] ...
        ,p2p,hands_radius); 
    [p2pCPD,matches,P1,P2] = LocalPat(M1,M2,[landmarks(5,1),landmarks(5,2)] ...
        ,p2pCPD,hands_radius); 
    new=M1.VERT(p2pCPD,:);
    M3=MESH('New',new,M2.TRIV);
    
    times = [times, toc(tCPD)];
    fprintf('done in %f seconds\n', times(end));
    


    % computes Arap between the input mesh which gives the connectivity and
    % the reconstructed mesh
    fprintf('processing arap for the mesh consistency...\n')

    timeArap = tic;
    for i=1:200
        [G1,E] = arap_gradient(M2.VERT,M2.TRIV,new);
        new=new-delta_t*(a_arap*G1);
    end
    new=MESH('New',new,M2.TRIV);
    
    times = [times, toc(timeArap)];
    fprintf('done in %f seconds\n', times(end));


    % Computes Arap+dato step
    fprintf(['Minimizing the ARAP energy to fit the original data...\n']);
    timeArapDato = tic;
    
    [new2, times] = arap_dato(new, M1, M2, delta_t, a_arap, times);
    
    fprintf('done in %f seconds\n', times(end));


end