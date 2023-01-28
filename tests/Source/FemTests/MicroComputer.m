classdef MicroComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        testName
        builderType
    end

    methods (Access = public)
        function obj = MicroComputer(cParams)
            obj.testName     = cParams.testName;
            obj.builderType  = cParams.builderType;
        end

        function compute(obj)
            a.fileName = obj.testName;
            s = FemDataContainer(a);
            s = ObjectSaver.saveObj(s);
            s.builderType = obj.builderType;
            femSolver = ElasticProblemMicro(s);
            femSolver.computeChomog();
            obj.computation = femSolver;
        end

    end

end
