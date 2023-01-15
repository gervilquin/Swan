classdef ReducedBuilder < handle

    properties (Access = private)
        bc
        lhs
        rhs
        solver
    end

    methods (Access = public)

        function createBuilder(obj, cParams)
            obj.init(cParams);
        end

        function [u,R] = solveSystem(obj)
            %Meter aqui el defRHS i el defLHS y hacerlos privados.
            defRHS = obj.createRHS();
            defLHS = obj.createLHS();
            sol = obj.solver.solve(defLHS, defRHS);
            u   = obj.bc.reducedToFullVector(sol);
            R   = obj.getReactions(sol);
        end
    end

    
    methods (Access = private)

        function init(obj, cParams)
            obj.bc     = cParams.bc;
            obj.lhs    = cParams.LHS;
            obj.rhs    = cParams.RHS;
            obj.solver = cParams.solver;
        end

        function defLHS = createLHS(obj)
            defLHS = obj.bc.fullToReducedMatrix(obj.lhs);
        end

        function defRHS = createRHS(obj)
            R       = obj.computeReactions();
            fullRHS = obj.rhs + R;
            defRHS  = obj.bc.fullToReducedVector(fullRHS);
        end

        function R = computeReactions(obj)
            boundaryCond  = obj.bc;
            K             = obj.lhs;
            dirich  = boundaryCond.dirichlet;
            dirichV = boundaryCond.dirichlet_values;
            if ~isempty(dirich)
                R = -K(:,dirich)*dirichV;
            else
                R = zeros(sum(obj.dim.ndofs(:)),1);
            end
        end

        function R = getReactions(obj, sol)
            K = obj.lhs;
            sizeK = size(K, 1);
            R = zeros(sizeK, 1);
            dirich  = obj.bc.dirichlet;
            dirichV = obj.bc.dirichlet_values;
            free = obj.bc.free;
            dl = sol;
            R(dirich) = K(dirich, free)*dl + K(dirich, dirich)*dirichV;
            R = R(dirich);
        end

    end
   
end