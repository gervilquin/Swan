classdef MatrixBuilder < handle
    % Class to build matrix for direct solution without reductions. 

    properties (Access = private)
        LHS
        bc
        RHS
        sizeK
        nConstraints
        solver
        scale
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
            R = computeReactions(obj, sol);
        end

    end

    methods (Access = private)
        function init(obj, cParams)
            obj.LHS   = cParams.LHS;
            obj.bc  = cParams.bc;
            obj.RHS = cParams.RHS;
            obj.sizeK  = size(cParams.LHS, 1);
            obj.solver = cParams.solver;
            obj.scale = cParams.scale;
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
            switch obj.scale
                case 'MACRO'
                    uD = obj.bc.dirichlet_values;
                    fullRHS = [obj.RHS; uD];
                case 'MICRO'
                    nPerDofs = size(obj.bc.periodic_constrained, 1);
                    perVector = zeros(nPerDofs, 1);
                    uD = obj.bc.dirichlet_values;
                    fullRHS = [obj.RHS; perVector; uD];
            end
            
        end

        function Ct = createConstraintMatrix(obj)
            [CtDir, sizeDir] = obj.computeDirichletCond();
            switch obj.scale
                case 'MACRO'
                    obj.nConstraints = sizeDir;
                    Ct               = CtDir;
                case 'MICRO'
                    perDOFslave      = obj.bc.periodic_constrained;
                    perDOFmaster     = obj.bc.periodic_free;
                    sizePer          = size(perDOFslave, 1);
                    obj.nConstraints = sizeDir + sizePer; 
                    Ct               = zeros(obj.nConstraints, obj.sizeK);
                    for i = 1:sizePer
                        masterNode = perDOFmaster(i);
                        slaveNode = perDOFslave(i);
                        Ct(i, [masterNode slaveNode]) = [1 -1];
                    end
                    Ct(sizePer+1:end, :) = CtDir;
            end
        end

        function [CtDir, sizeDir] = computeDirichletCond(obj)
            dirDOFs    = obj.bc.dirichlet;
            sizeDir = size(dirDOFs, 1);
            CtDir = zeros(sizeDir, obj.sizeK);
            for i = 1:sizeDir
                CtDir(i,dirDOFs(i)) = 1; 
            end
        end

        function R = computeReactions(obj, sol)
            switch obj.scale
                case 'MACRO'
                    R = -sol(obj.sizeK+1:end, 1);
                case 'MICRO'
                    sizePer = size(obj.bc.periodic_constrained, 1);
                    R = -sol(obj.sizeK+sizePer+1:end, 1);
            end
        end
    end
end