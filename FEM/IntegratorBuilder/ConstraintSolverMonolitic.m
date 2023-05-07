classdef ConstraintSolverMonolitic < ConstraintSolverFactory

    properties (Access = private)
        solver
        bc
        K
        RHS
        sizeK
        vstrain
        sizePer
    end

    methods (Access = public)
        function obj = ConstraintSolverMonolitic(cParams)
            obj.init(cParams);
        end

        function [u, L] = solveSystem(obj, LHSMatrix, RHSMatrix, nConstraints)
            lhs = obj.createGeneralMatrix(LHSMatrix, nConstraints);
            sol         = obj.solver.solve(lhs, RHSMatrix);
            u = sol(1:obj.sizeK, 1);
            L = sol(obj.sizeK+1:end, 1);
         
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end

    methods (Access = private)
        function init(obj, cParams)
            obj.solver = cParams.solver;
            obj.bc     = cParams.bc; 
            obj.K      = cParams.LHS;
            obj.sizeK  = size(obj.K, 1);
            obj.RHS    = cParams.RHS;
            if isfield(cParams, 'vstrain')
                obj.vstrain = cParams.vstrain;
                perDOFslave = obj.bc.periodic_constrained;
                obj.sizePer = size(perDOFslave, 1);
            end 
        end

        function fullLHS = createGeneralMatrix(obj, Ct, nConstraints)
            C       = Ct';
            nC      = nConstraints;
            Z       = zeros(nC);
            Km      = obj.K;
            fullLHS = [Km C; C' Z];
        end

    end


end