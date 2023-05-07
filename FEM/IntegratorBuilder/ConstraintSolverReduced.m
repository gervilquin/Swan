classdef ConstraintSolverReduced < ConstraintSolverFactory

    properties 
        solver
        bc
        K
        RHS
        sizeK
        vstrain
        sizePer
        nConstraints
        solType
        solMode
        scale
    end

    methods (Access = public)
        function obj = ConstraintSolverReduced(cParams)
            obj.init(cParams);
        end

        function [u, L] = solveSystem(obj, lhs, rhs, nConstraints)
            sol         = obj.solver.solve(lhs, rhs);
            u = reducedToFullVector(obj, sol);
            if strcmp(obj.scale, 'MACRO')
                L = obj.getReactions(sol);
            else
                L = 0;
            end
        end
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.solver = cParams.solver;
            obj.bc     = cParams.bc; 
            obj.K      = cParams.LHS;
            obj.sizeK  = size(obj.K, 1);
            obj.RHS    = cParams.RHS;
            obj.scale  = cParams.scale;
            if isfield(cParams, 'vstrain')
                obj.vstrain = cParams.vstrain;
                perDOFslave = obj.bc.periodic_constrained;
                obj.sizePer = size(perDOFslave, 1);
            end 
            obj.solType = cParams.solType;
            obj.solMode = cParams.solMode;
        end
        
        function full = reducedToFullVector(obj, vec)
            switch obj.scale
                case 'MACRO'
                    full = obj.expandVectorDirichlet(vec);
                case 'MICRO'
                    full = obj.expandVectorPeriodic(vec);
            end
        end

        function b = expandVectorDirichlet(obj,bfree)
            dir = obj.bc.dirichlet;
            uD  = obj.bc.dirichlet_values;
            fr  = obj.bc.free;
            nsteps = length(bfree(1,:));
            ndof = sum(obj.bc.ndofs);
            uD = repmat(uD,1,nsteps);
            
            b = zeros(ndof,nsteps);
            b(fr,:) = bfree;
            if ~isempty(dir)
                b(dir,:) = uD;
            end
        end

        function b = expandVectorPeriodic(obj,bfree)
            vF = obj.bc.free;
            vP = obj.bc.periodic_free;
            vC = obj.bc.periodic_constrained;
            vI = setdiff(vF,vP);
            b = zeros(obj.bc.ndofs,1);
            b(vI) = bfree(1:1:size(vI,2));
            b(vP) = bfree(size(vI,2)+1:1:size(bfree,1));
            b(vC) = b(vP);
        end

        function R = getReactions(obj, sol)
            R         = zeros(obj.sizeK, 1);
            dirich    = obj.bc.dirichlet;
            dirichV   = obj.bc.dirichlet_values;
            free      = obj.bc.free;
            dl        = sol;
            R(dirich) = obj.K(dirich, free)*dl + obj.K(dirich, dirich)*dirichV;
            R         = R(dirich);
        end

    
    end

end