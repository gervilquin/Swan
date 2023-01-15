classdef MatrixBuilder < handle
    % Class to build matrix for direct solution without reductions. 

    properties (Access = private)
        LHS
        bc
        RHS
        sizeK
        nConstraints
        solver
    end

    methods (Access = public)
        function createBuilder(obj, s)
            obj.init(s);
        end

        function [u, R] = solveSystem(obj)
            %Meter aqui el defRHS i el defLHS y hacerlos privados.
            defRHS = obj.createRHS();
            defLHS = obj.createLHS();
            sol = obj.solver.solve(defLHS, defRHS);
            u   = sol(1:obj.sizeK, 1);
            R   = -sol(obj.sizeK+1:end, 1);
        end

    end

    methods (Access = private)
        function init(obj, cParams)
            obj.LHS   = cParams.LHS;
            obj.bc  = cParams.bc;
            obj.RHS = cParams.RHS;
            obj.sizeK  = size(cParams.LHS, 1);
            obj.solver = cParams.solver;
        end

        function defLHS = createLHS(obj)
            Ct     = obj.createConstraintMatrix();
            defLHS = obj.createGeneralMatrix(Ct);
        end

        function defRHS = createRHS(obj)
            defRHS = obj.createGeneralVector();
        end

        function fullLHS = createGeneralMatrix(obj, Ct)
            C = Ct';
            nC = obj.nConstraints;
            Z  = zeros(nC);
            Km = obj.LHS;
            fullLHS = [Km C; C' Z];
        end

        function fullRHS = createGeneralVector(obj)
            uD = obj.bc.dirichlet_values;
            fullRHS = [obj.RHS; uD];
        end

        function Ct = createConstraintMatrix(obj)
            dirichletDOFs    = obj.bc.dirichlet;
            obj.nConstraints = size(dirichletDOFs, 1);
            Ct               = zeros(obj.nConstraints, obj.sizeK);
            for i = 1:obj.nConstraints
                 Ct(i,dirichletDOFs(i)) = 1; 
            end     
        end
    end
end