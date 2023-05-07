classdef ConstraintSolverFactory < handle

    properties
        Property1
    end

    methods (Access = public, Static)
        function obj = create(cParams)
            switch cParams.solType
                case 'REDUCED'
                    obj = ConstraintSolverReduced(cParams);
                case 'MONOLITIC'
                    obj = ConstraintSolverMonolitic(cParams);
            end
        end

    end
end