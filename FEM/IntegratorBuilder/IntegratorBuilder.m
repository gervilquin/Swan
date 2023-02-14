classdef IntegratorBuilder < handle

    methods (Static)
        function iType = create(cParams)
            switch cParams.btype
                case 'MONOLITIC'
                    iType = MatrixBuilder();

                case 'REDUCED'
                    iType = ReducedBuilder();

                case 'MONOLITIC_MICRO'
                    iType = MicroBuilder();
            end
            
        end

    end
end