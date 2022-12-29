classdef MatrixBuilder < handle
    % Class to build matrix for direct solution without reductions. 

    properties (Access = private)
        K
        bc
        RHS
        sK
    end

    methods (Access = public)
        function obj = MatrixBuilder(s)
            obj.init(s);
        end

        function genLHS = createLHS(obj)
            Ct = obj.createConstraintMatrix;
            genLHS = obj.createGeneralMatrix(Ct);
        end

        function genRHS = createRHS(obj)
            genRHS = obj.createGeneralVector();
        end
    end

    methods (Access = private)
        function init(obj, cParams)
            obj.K = cParams.K;
            obj.bc = cParams.bc;
            obj.RHS = cParams.RHS;
            obj.sK = size(cParams.K, 1);
        end

        function generalM = createGeneralMatrix(obj, Ct)
            % Method to create the general matrix
            C = Ct';
            sC1 = size(C, 1);
            sC2 = size(C, 2);
            generalM = zeros(obj.sK+sC2);
            generalM(1:obj.sK, 1:obj.sK) = obj.K;
            generalM(obj.sK+1:end, 1:sC1) = C';
            generalM(1:sC1, obj.sK+1:end) = C;
        end

        function fullRHS = createGeneralVector(obj)
            constDisp = obj.bc.dirichlet_values;
            sRHS = size(obj.RHS,1);
            sDirichlet = size(constDisp, 1);
            fullRHS = zeros(sRHS+sDirichlet, 1);
            fullRHS(1:sRHS) = obj.RHS;
            fullRHS(sRHS+1:end) = constDisp;
        end

        function Ct = createConstraintMatrix(obj)
            dirichletDOFs = obj.bc.dirichlet;
            nDOFr = size(dirichletDOFs, 1);
            Ct = zeros(nDOFr, obj.sK);
            for i = 1:nDOFr
                 Ct(i,dirichletDOFs(i)) = 1; 
            end     
        end
    end
end