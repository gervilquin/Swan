classdef FlucDispChanger < handle

    properties
        
    end

    methods (Access = public)
        function obj = FlucDispChanger(cParams)
           obj.init(cParams);
        end

        function uTotal = computeTotalDisplacements(obj, u)
            coords  = obj.mesh.coord';
            nel     = size(coords, 2);
            strainM = [obj.vstrain(1) obj.vstrain(3); 
                        obj.vstrain(3) obj.vstrain(2)];
            uTotal  = zeros(obj.sizeK, 1);
            for i=1:nel
                uTotal([2*i-1 2*i], 1) = strainM*coords(:, i);
            end
        end

    end

    methods (Access = private)
    
    end
    
end