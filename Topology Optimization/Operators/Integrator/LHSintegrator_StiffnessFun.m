classdef LHSintegrator_StiffnessFun < handle %LHSintegrator

    properties (Access = private)
        fun
        mesh
        material
        quadrature
    end

    methods (Access = public)

        function obj = LHSintegrator_StiffnessFun(cParams)
            obj.initFun(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            dNdx  = obj.fun.computeCartesianDerivatives(obj.quadrature,obj.mesh);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = size(obj.mesh.connec,1);
            nNodE = size(dNdx,2);
            nDofE = nNodE*obj.fun.ndimf;
            lhs = zeros(nDofE,nDofE,nElem);
            Bcomp = obj.createBComputer(dNdx);
            for igaus = 1:nGaus
                Bmat = Bcomp.compute(igaus);
                dV(1,1,:) = dVolu(igaus,:)';
                Bt  = permute(Bmat,[2 1 3]);
                BtB = pagemtimes(Bt, Bmat);
                lhs = lhs + bsxfun(@times, BtB, dV);
            end
        end

    end

    methods (Access = private)

        function initFun(obj, cParams)
            obj.fun      = cParams.fun;
            obj.mesh     = cParams.mesh;
        end
        
       function createQuadrature(obj)
           quad = Quadrature.set(obj.mesh.type);
           quad.computeQuadrature(obj.fun.order);
           obj.quadrature = quad;
       end

        function Bcomp = createBComputer(obj, dNdx)
            s.fun  = obj.fun;
            s.dNdx = dNdx;
            Bcomp = BMatrixComputerFun(s);
        end

        function LHS = assembleMatrix(obj, lhs)
            s.connec = obj.mesh.connec; % !!!
            s.fun    = obj.fun; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs);
        end
    end

end