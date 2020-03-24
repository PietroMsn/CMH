function [M1, M2] = scaling_shapes(M1, M2)
 
    load('./utils/scaling_weights');
    
    factor = full(sqrt(sum(scaling_w(:))/sum(M1.A(:))));
    M1.VERT = factor.*M1.VERT;

    factor = full(sqrt(sum(scaling_w(:))/sum(M2.A(:))));
    M2.VERT = factor.*M2.VERT;
    
end