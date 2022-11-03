classdef FunElasticProblem < handle
    
    % Just a concept. Does not work.
    properties (Access = public)
        
    end

    properties (Access = private) % Calculated
        LHS, RHS, solver
        displacement
        strain, stress
    end

    properties (Access = private) % Inputs
        pdim
        mesh
        scale
        material
        BC
    end

    methods (Access = public)

        function obj = FunElasticProblem(cParams)
            obj.init(cParams);
            obj.createDisplacementFunction();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeStiffnessMatrix();
            obj.computeForces();
            obj.computeDisplacements();
            obj.computeStrain();
%             obj.computeStress();
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.pdim     = cParams.dim;
            obj.mesh     = cParams.mesh;
            obj.scale    = cParams.scale;
            obj.material = cParams.material;
            obj.BC       = cParams.bc;
        end

        function createDisplacementFunction(obj)
            strdim = regexp(obj.pdim,'\d*','Match');
            nDimf  = str2double(strdim);
            nNodes = size(obj.mesh.coord,1);
            s.ndimf   = nDimf;
            s.connec  = obj.mesh.connec;
            s.type    = obj.mesh.type;
            s.fValues = zeros(nNodes,nDimf);
            f = P1Function(s);
            obj.displacement = f;
        end

        function createSolver(obj)
%             s.type =  'DIRECT';
            s.BC   = obj.BC;
            s.mesh = obj.mesh;
            obj.solver = FunSolver(s);
        end

        function computeStiffnessMatrix(obj)
            s.type     = 'FunElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = obj.displacement;
            s.material = obj.material;
            lhs = LHSintegrator.create(s);
            obj.LHS = lhs.compute();
        end

        function computeForces(obj)
            s.type    = 'FunElastic';
            s.mesh    = obj.mesh;
            s.fun     = obj.displacement;
            s.neumann = obj.BC.neumann;
            rhs = RHSintegrator_FunElastic(s);
            obj.RHS = rhs.compute();
        end

        function u = computeDisplacements(obj)
            u = obj.solver.solve(obj.LHS,obj.RHS);
            nDimf = obj.displacement.ndimf;
            nNods = size(obj.mesh.coord,1);
            uRsh = reshape(u,[nDimf, nNods])';
            obj.displacement.fValues = uRsh;
        end

        function computeStrain(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            sig = obj.displacement.computeSymmetricGradient2(q,obj.mesh);
            obj.strain = sig;
        end

        function computeStress(obj)
            obj.stress = pagemtimes(obj.material.C, obj.strain);
        end
        

    end

end
