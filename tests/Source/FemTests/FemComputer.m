classdef FemComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        testName
        interpolationType
        solType
        solMode
    end

    methods (Access = public)
        function obj = FemComputer(cParams)
            obj.testName          = cParams.testName;
            obj.solType = cParams.solType;
            obj.solMode = cParams.solMode;
            if isfield(cParams, 'interpolationType')
                obj.interpolationType = cParams.interpolationType;
            else
                obj.interpolationType = 'LINEAR';
            end
        end

        function compute(obj)
            a.fileName = obj.testName;
            s = FemDataContainer(a);
            s.interpolationType = obj.interpolationType;
            s = ObjectSaver.saveObj(s);
            s.solType = obj.solType;
            s.solMode = obj.solMode;
%             s.builderType = obj.builderType;
            obj.computation = FEM.create(s);
            obj.computation.solve();
        end
    end

end