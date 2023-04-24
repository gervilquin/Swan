classdef MicroComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        testName
        builderType
        solType
        solMode
    end

    methods (Access = public)
        function obj = MicroComputer(cParams)
            obj.testName     = cParams.testName;
%             obj.builderType  = cParams.builderType;
            obj.solType = cParams.solType;
            obj.solMode = cParams.solMode;
        end

        function compute(obj)
            a.fileName = obj.testName;
            s = FemDataContainer(a);
            s = ObjectSaver.saveObj(s);
%             s.builderType = obj.builderType;
            s.solType = obj.solType;
            s.solMode = obj.solMode;
            switch s.solMode
                case 'DISP'
                    femSolver = ElasticProblemDisp(s);
                    femSolver.computeStressHomog();
                    obj.computation = femSolver;
                case 'FLUC'
                    switch s.solType
                        case 'REDUCED'
                            femSolver = ElasticProblemMicro(s);
                            femSolver.computeChomog();
                            obj.computation = femSolver;
                        case 'MONOLITIC'
                            femSolver = ElasticProblemFluc(s);
                            femSolver.computeChomog();
                            obj.computation = femSolver;
                    end

            end
%             femSolver = ElasticProblemMicro(s);
%             femSolver.computeChomog();
%             obj.computation = femSolver;
        end

    end

end
