classdef ProjectedCostDerivator < handle 
    properties (Access = public)
        derivedCost
    end
    properties (Access = private)
        structure
        mesh 
        displacement 
        
    end
    methods (Access = public)
        function obj = ProjectedCostDerivator(cParams)
            obj.inputData(cParams);
        end
        function compute(obj)
            obj.computeCost();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.mesh.conectivityMatrixMat = cParams.conectivityMatrixMat;
            obj.mesh.elementNumberX = cParams.elementNumberX;
            obj.mesh.elementNumberY = cParams.elementNumberY;
            obj.structure.elementalStiffnessMatrix = cParams.elementalStiffnessMatrix;
            obj.displacement = cParams.displacement; 
            
        end 
        function computeCost(obj) 
           obj.derivedCost  = reshape(sum((obj.displacement(obj.mesh.conectivityMatrixMat)*obj.structure.elementalStiffnessMatrix).*obj.displacement(obj.mesh.conectivityMatrixMat),2),obj.mesh.elementNumberY,obj.mesh.elementNumberX);
        end 
    end
end