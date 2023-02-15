classdef CMatrixComputer_ElasticMicro < handle

    properties (Access = private)
        dim
        quadrature
        mesh
        geometry
    end

    methods (Access = public)

        function obj = CMatrixComputer_ElasticMicro(cParams)
           obj.init(cParams); 
           obj.createQuadrature();
           obj.createGeometry();
        end
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.dim = cParams.dim;
            obj.quadrature = cParams.quadrature;
            obj.mesh = cParams.mesh;
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.quadrature = q;
        end

        function createGeometry(obj)
            q = obj.quadrature;
            int = obj.mesh.interpolation;
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end

        function Bmat = computeBMatrix(obj, igaus)
            s.dim      = obj.dim;
            s.geometry = obj.geometry;
            Bcomp      = BMatrixComputer(s);
            Bmat       = Bcomp.compute(igaus);
        end
    
    end
end