classdef BoundaryCondTestsSuite < handle & matlab.unittest.TestCase
    
    methods (Access = public)
        function obj = BoundaryCondTestsSuite()
            path = './tests/Source/FemTests/BoundaryConditionsTests/BoundaryCondTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path);
            results = suite.run;
            table(results)
        end
    end
end