classdef MatrixBuilder < handle
    % Class to build matrix for direct solution without reductions. 

    properties (Access = private)
        K
        bc
        RHS
        sK
        nConstraints
    end

    methods (Access = public)
        function createBuilder(obj, s)
            obj.init(s);
        end

        function defLHS = createLHS(obj)
            Ct     = obj.createConstraintMatrix;
            defLHS = obj.createGeneralMatrix(Ct);
        end

        function defRHS = createRHS(obj)
            defRHS = obj.createGeneralVector();
        end
    end

    methods (Access = private)
        function init(obj, cParams)
            obj.K   = cParams.LHS;
            obj.bc  = cParams.bc;
            obj.RHS = cParams.RHS;
            obj.sK  = size(cParams.LHS, 1);
        end

        function fullLHS = createGeneralMatrix(obj, Ct)
            % Method to create the general matrix
            C = Ct';
%             sC1 = size(C, 1);
%             sC2 = size(C, 2);
%             fullLHS = zeros(obj.sK+sC2);
%             fullLHS(1:obj.sK, 1:obj.sK) = obj.K;
%             fullLHS(obj.sK+1:end, 1:sC1) = C';
%             fullLHS(1:sC1, obj.sK+1:end) = C;
            nC = obj.nConstraints;
            Z  = zeros(nC);
            Km = obj.K;
            fullLHS = [Km C; C' Z];
        end

        function fullRHS = createGeneralVector(obj)
            uD = obj.bc.dirichlet_values;
%             sRHS = size(obj.RHS,1);
%             sDirichlet = size(constDisp, 1);
%             fullRHS = zeros(sRHS+sDirichlet, 1);
%             fullRHS(1:sRHS) = obj.RHS;
%             fullRHS(sRHS+1:end) = constDisp;
            fullRHS = [obj.RHS; uD];
        end

        function Ct = createConstraintMatrix(obj)
            bC = obj.bc;
            dirichletDOFs    = bC.dirichlet;
            obj.nConstraints = size(dirichletDOFs, 1);
            Ct               = zeros(obj.nConstraints, obj.sK);
            for i = 1:obj.nConstraints
                 Ct(i,dirichletDOFs(i)) = 1; 
            end     
        end
    end
end