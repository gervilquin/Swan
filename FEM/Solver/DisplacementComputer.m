classdef DisplacementComputer < handle

    methods (Static)

        function u = computeU(cParams, sol)
            switch cParams.type
                case 'MONOLITIC'
                    uDofs = size(cParams.LHS, 1);
                    u = sol(1:uDofs);
                case 'REDUCED'
                    boundaryCond = cParams.bc;
                    u = boundaryCond.reducedToFullVector(sol);
            end
        end
    end
end