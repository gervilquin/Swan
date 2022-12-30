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
            defRHS = obj.bc.fullToReducedVector(obj.rhs);
        end
    end

    
    methods (Access = private)

        function init(obj, cParams)
            obj.bc = cParams.bc;
            obj.lhs = cParams.LHS;
            obj.rhs = cParams.RHS;
        end
    end
end