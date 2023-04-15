classdef IntegratorBuilder < handle

    methods (Static)
        function iType = create(cParams)
            switch cParams.btype
                case 'MONOLITIC_DIR'
                    iType = DirichletBuilder();

                case 'REDUCED'
                    iType = ReducedBuilder();

                case 'MONOLITIC_MICRO'
                    iType = MicroBuilder();
                
                case 'MONOLITIC_MICRO_CoV'
                    iType = MicroBuilderCoV();
            end
            
        end

    end
end