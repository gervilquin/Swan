classdef FunThermalProblem < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = FunThermalProblem(cParams)
            obj.init(cParams)
            obj.createTemperatureFunction();
            obj.createSolver();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
        end
        
        function createTemperatureFunction(obj)
        end

        function createSolver(obj)
        end

    end
    
end