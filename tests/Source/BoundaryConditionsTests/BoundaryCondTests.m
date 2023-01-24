classdef BoundaryCondTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        triangle = {'test2d_triangle_non_null'}
    end

    methods (Test)
        function testTriangle(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = triangle;
            s.variablesToStore = {'d_u'};
            s.testSolverType   = 'MONOLITIC';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

    end
end