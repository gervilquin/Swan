classdef BoundaryCondition < handle
    
    properties (Access = public)
        type
        domainFun
        value
    end
    
    properties (Access = private)
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = BoundaryCondition(cParams)
            obj.init(cParams)
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            % We are not specifying which DOF is being restricted. For now,
            % let's assume we restrict all of them at the same time
            obj.type      = cParams.type;
            obj.value     = cParams.value;
            obj.domainFun = cParams.domainFun;
        end
        
    end
    
end