classdef RHSintegrator_FunElastic < handle

    properties (Access = private)
        mesh
        neumann
        fun
    end
    
    methods (Access = public)

        function obj = RHSintegrator_FunElastic(cParams)
            obj.init(cParams);
        end

        function Fext = compute(obj)
            Fext = obj.computePunctualFext();
        end
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.fun     = cParams.fun;
            obj.mesh    = cParams.mesh;
            obj.neumann = cParams.neumann;
        end

        function F = computePunctualFext(obj)
            %Compute Global Puntual Forces (Not well-posed in FEM)
            nNods = size(obj.mesh.coord,1);
            nDimf = obj.fun.ndimf;
            Fp = zeros(nNods, nDimf);
            coor = obj.mesh.coord';
            isFNod = obj.neumann.domainFun(coor)';
            Fp(isFNod,:) = obj.neumann.value;
            F = reshape(Fp', [nNods*nDimf, 1]);
        end

    end
    
end