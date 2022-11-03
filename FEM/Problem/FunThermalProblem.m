classdef FunThermalProblem < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private) % In
        mesh
        BC
    end
    
    properties (Access = private) % Calc
        LHS, RHS, solver
        temperature
    end
    
    methods (Access = public)
        
        function obj = FunThermalProblem(cParams)
            obj.init(cParams)
            obj.createTemperatureFunction();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeLHS();
            obj.computeRHS();
            obj.computeTemperature();
            obj.temperature.plot(obj.mesh);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.BC   = cParams.bc;
        end
        
        function createTemperatureFunction(obj)
            nDimf  = 1;
            nNodes = size(obj.mesh.coord,1);
            s.ndimf   = nDimf;
            s.connec  = obj.mesh.connec;
            s.type    = obj.mesh.type;
            s.fValues = zeros(nNodes,nDimf);
            f = P1Function(s);
            obj.temperature = f;
        end

        function createSolver(obj)
%             s.type = 'DIRECT';
            s.mesh = obj.mesh;
            s.BC   = obj.BC;
            obj.solver = FunSolver(s);
        end

        function computeLHS(obj)
            cond = 10;
            s.type  = 'FunStiffnessMatrix';
            s.mesh  = obj.mesh;
            s.fun   = obj.temperature;
            lhs = LHSintegrator.create(s);
            obj.LHS = cond*lhs.compute();
        end

        function computeRHS(obj)
            fNodal         = obj.createHeatSource();
                             obj.computeBoundaryContribution();
            s.type         = 'ShapeFunction';
            s.mesh         = obj.mesh;
            s.npnod        = size(obj.mesh.coord,1);
            s.globalConnec = obj.mesh.connec;
            rhs = RHSintegrator.create(s);
            obj.RHS = rhs.compute(fNodal);
%             nNodes = size(obj.mesh.coord,1);
%             obj.RHS = rand(nNodes,1);
        end

        function fNod = createHeatSource(obj)
            sAF.fHandle = @(x) [ x(1,:,:)*0 + 4 ];
            sAF.ndimf   = 1;
            sAF.mesh    = obj.mesh;
            xFun = AnalyticalFunction(sAF);
            
            pp1.mesh   = obj.mesh;
            pp1.connec = obj.mesh.connec;
            projP1 = Projector_toP1(pp1);
            p1fun = projP1.project(xFun);
            fNod = p1fun.fValues;
        end

        function computeBoundaryContribution(obj)
%             fBound = zeros(size(fNodal));
%             q=1e2;
%             fBound(1,:)  = q;
%             fBound(13,:) = q;
        end

        function computeTemperature(obj)
            t = obj.solver.solve(obj.LHS,obj.RHS);
            obj.temperature.fValues = t;
        end

    end
    
end