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
            s.testSolverType   = 'MONOLITIC';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

        function testTriangleMonolitic(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = [triangle '_non_null'];
            s.variablesToStore = {'d_u'};
            s.testSolverType   = 'MONOLITIC';
            s.testResultsName  = [triangle '_non_null'];
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

    end

    methods (Test, TestTags = {'Monolitic', 'Micro'})

        function testMicroReduced(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.testSolverType = 'REDUCED';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end
        
        function testMicroFluc(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.testSolverType   = 'MONOLITIC_MICRO';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testMicroChOV(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.testSolverType   = 'MONOLITIC_MICRO_CoV';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end
        
    end
end