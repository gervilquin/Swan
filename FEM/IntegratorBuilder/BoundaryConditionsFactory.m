classdef BoundaryConditionsFactory < handle

    properties
       
    end

    methods (Access = public, Static)
        function obj = create(cParams)
            switch cParams.solMode
                case 'DISP'
                    obj = BoundaryConditionsDisp(cParams);
                case 'FLUC'
                    obj = BoundaryConditionsFluc(cParams);
            end
        end

    end
end