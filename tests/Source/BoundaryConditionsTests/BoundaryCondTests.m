classdef BoundaryCondTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        triangle = {'test2d_triangle'}
    end

    methods (Test)

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
end