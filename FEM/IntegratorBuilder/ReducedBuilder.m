classdef ReducedBuilder < handle

    properties (Access = private)
        bc
        lhs
        rhs
    end

    methods (Access = public)

        function createBuilder(obj, cParams)
            obj.init(cParams);
        end

        function defLHS = createLHS(obj)
            defLHS = obj.bc.fullToReducedMatrix(obj.lhs);
        end

        function defRHS = createRHS(obj)
            R = obj.computeReactions();
            fullRHS = obj.rhs + R;
            defRHS = obj.bc.fullToReducedVector(fullRHS);
        end
    end

    
    methods (Access = private)

        function init(obj, cParams)
            obj.bc = cParams.bc;
            obj.lhs = cParams.LHS;
            obj.rhs = cParams.RHS;
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

    end
   
end