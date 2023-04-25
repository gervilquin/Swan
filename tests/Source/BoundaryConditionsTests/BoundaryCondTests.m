classdef BoundaryCondTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        triangle = {'test2d_triangle'}
%         micro = {'test2d_micro', 'test2d_micro_thin',...
%             'test_micro_holeinclusion'}
        micro ={'test_micro_holeinclusion'}
    end

    methods (Test, TestTags = {'Monolitic', 'Macro'})

        function testTriangleNullDisp(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = triangle;
            s.variablesToStore = {'d_u'};
            s.testSolverType   = 'MONOLITIC_DIR';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

        function testTriangleMonolitic(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = [triangle '_non_null'];
            s.variablesToStore = {'d_u'};
            s.testSolverType   = 'MONOLITIC_DIR';
            s.testResultsName  = [triangle '_non_null'];
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

    end

    methods (Test, TestTags = {'Monolitic', 'Micro'})

        function testMicroFlucReduced(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.solType   = 'REDUCED';
            s.solMode   = 'FLUC';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-4;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testMicroDispMonolitic(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.solType   = 'MONOLITIC';
            s.solMode   = 'DISP';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-4;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testMicroFlucMonolitic(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.solType   = 'MONOLITIC';
            s.solMode   = 'FLUC';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-4;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end
end