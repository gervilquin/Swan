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

    end
end