%% LAGRANGE SIMPLICIAL
k = 2;
d = 2;
ls = LagrangeElement.create('SIMPLICIAL',k,d);
ls.plotShapeFunctions();

%% LAGRANGE TENSOR PRODUCT
k = 2;
d = 2;
lt = LagrangeElement.create('TENSOR PRODUCT',k,d);
lt.plotShapeFunctions();

%% CROUZEIX-RAVIART
d = 2;
cr = CrouzeixRaviartElement.create(d);
cr.plotShapeFunctions();

%% RAVIART-THOMAS
d = 2;
rt = RaviartThomasElement.create(d);
rt.plotShapeFunctions();

%% NÉDÉLEC
d = 2;
nd = NedelecElement.create(d);
nd.plotShapeFunctions();