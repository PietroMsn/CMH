function CMH = compute_CMH(shape,basis)

error = shape.VERT - full(basis*basis'*shape.Ae*shape.VERT);

CMH = basis;

for i =1:3
    Mcmh = error(:,i) - CMH*CMH'*shape.Ae*error(:,i);
    Mcmh = Mcmh./sqrt(Mcmh'*shape.Ae*Mcmh);
    CMH = [CMH,Mcmh]; clear Mcmh;
end


