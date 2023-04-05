classdef MicroBuilder < handle

    properties (Access = private)
        LHS
        bc
        RHS
        sizeK
        nConstraints
        solver
        scale
        sizePer
        mesh
        vstrain
    end

    methods (Access = public)
        function createBuilder(obj, s)
            obj.init(s);
        end

        function [u, stressHomog] = solveSystem(obj)
            defLHS = obj.createLHS();
            defRHS = obj.createRHS();
            sol = obj.solver.solve(defLHS, defRHS);
            u   = sol(1:obj.sizeK, 1);
            [L, LDir, stressHomog] = obj.computeStressHomog(sol);
            uTotal = obj.computeTotalDisplacements(u);
        end

    end

    methods (Access = private)
        function init(obj, cParams)
            obj.LHS     = cParams.LHS;
            obj.bc      = cParams.bc;
            obj.RHS     = cParams.RHS;
            obj.sizeK   = size(cParams.LHS, 1);
            obj.solver  = cParams.solver;
            obj.scale   = cParams.scale;
            obj.mesh    = cParams.mesh;
            obj.vstrain = cParams.vstrain;
            perDOFslave = obj.bc.periodic_constrained;
            obj.sizePer = size(perDOFslave, 1);
        end

        function defLHS = createLHS(obj)
            Ct     = obj.createConstraintMatrix();
            defLHS = obj.createGeneralMatrix(Ct);
        end

        function defRHS = createRHS(obj)
            defRHS = obj.createGeneralVector();
        end

        function fullLHS = createGeneralMatrix(obj, Ct)
            C = Ct';
            nC = obj.nConstraints;
            Z  = zeros(nC);
            Km = obj.LHS;
            fullLHS = [Km C; C' Z];
        end

        function fullRHS = createGeneralVector(obj)
            nPerDofs = size(obj.bc.periodic_constrained, 1);
            perVector = zeros(nPerDofs, 1);
            uD = obj.bc.dirichlet_values;
            fullRHS = [obj.RHS; perVector; uD]; 
        end

        function Ct = createConstraintMatrix(obj)
            s.dirDOFs        = obj.bc.dirichlet;
            s.sizeK          = obj.sizeK;
            DirComputer      = DirichletComputer(s);
            [CtDir, sizeDir] = DirComputer.compute();
            perDOFslave      = obj.bc.periodic_constrained;
            perDOFmaster     = obj.bc.periodic_free;
            obj.nConstraints = sizeDir + obj.sizePer; 
            Ct               = zeros(obj.nConstraints, obj.sizeK);
            for i = 1:obj.sizePer
                masterNode = perDOFmaster(i);
                slaveNode = perDOFslave(i);
                Ct(i, [masterNode slaveNode]) = [1 -1];
            end
            Ct(obj.sizePer+1:end, :) = CtDir;
        end

        function [L,  LDir, stressHomog] = computeStressHomog(obj, sol)
            nEqperType = obj.sizePer/4;
            sigmaX = 0;
            sigmaY = 0;
            tauXY  = 0;
            d1  = obj.sizeK+1;
            d2  = obj.sizeK + obj.sizePer;
            L   = sol(d1:d2);
            LDir = sol(d2+1:end);
            for i = 1:nEqperType
                sigmaX = sigmaX + L(i);
            end
            for i = nEqperType+1:2*nEqperType
                tauXY = tauXY + L(i);
            end
            for i = 2*nEqperType+1:3*nEqperType
                tauXY = tauXY + L(i);
            end
            for i = 3*nEqperType+1:4*nEqperType
                sigmaY = sigmaY + L(i);
            end
            sigmaX = sigmaX + LDir(1) + LDir(5);
            sigmaY = sigmaY + LDir(2) + LDir(4);
            tauXY = (tauXY + LDir(6) + LDir(3) - LDir(7) - LDir(8))/2;
            stressHomog = [sigmaX; sigmaY; tauXY];
        end

        function uTotal = computeTotalDisplacements(obj, u)
            coords = obj.mesh.coord';
            nel = size(coords, 2);
            strainM = [obj.vstrain(1) obj.vstrain(3); 
                obj.vstrain(3) obj.vstrain(2)];
            uTotal = zeros(obj.sizeK, 1);
            for i=1:nel
                uTotal([2*i-1 2*i], 1) = strainM*coords(:, i);
            end
        end

    end
end