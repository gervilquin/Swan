classdef IntegratorBuilder < handle

    methods (Static)
        function iType = create(cParams)
            switch cParams.type
                case 'MONOLITIC'
                    iType = MatrixBuilder();

                case 'REDUCED'
                    iType = ReducedBuilder();
            end
            
        end

        function u = computeU(bc, cParams, sol)
            switch cParams.type
                case 'MONOLITIC'
                    uDofs = size(cParams.LHS, 1);
                    u = sol(1:uDofs);
                case 'REDUCED'
                    u = bc.reducedToFullVector(u);

            end
        end
    end
end