classdef MicroBuilder < handle

    properties (Access = private)
        LHS
        bc
        sizeK
        nConstraints
        solver
        scale
        mesh
        sizeDir
        sizePer
        vstrain
    end

    methods (Access = public)
        function createBuilder(obj, cParams)
            obj.init(cParams);
        end

        function [u, stressHomog] = solveSystem(obj)
            defLHS      = obj.createLHS();
            defRHS      = obj.createRHS();
            sol         = obj.solver.solve(defLHS, defRHS);
            u           = sol(1:obj.sizeK, 1);
            [L, stressHomog] = obj.computeStressHomog(sol);
        end

    end

    methods (Access = private)
        function init(obj, cParams)
            obj.LHS     = cParams.LHS;
            obj.bc      = cParams.bc;
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
            C       = Ct';
            nC      = obj.nConstraints;
            Z       = zeros(nC);
            Km      = obj.LHS;
            fullLHS = [Km C; C' Z];
        end
        
        function fullRHS = createGeneralVector(obj)
            nU        = size(obj.LHS, 1);
            der1Rhs   = zeros(nU, 1);
            uD        = obj.bc.dirichlet_values;
            nEqperType = obj.sizePer/4;
            perVector    = zeros(nEqperType*3, 1);
            for i = 1:nEqperType
                perVector(i, 1) = obj.vstrain(1);
            end
            % second 39
            for i = nEqperType+1:2*nEqperType
                perVector(i, 1) = obj.vstrain(3);
            end
            % third 39
            for i = 2*nEqperType+1:3*nEqperType
                perVector(i, 1) = obj.vstrain(2);
            end
            fullRHS   = [der1Rhs; perVector; uD];
        end

        function Ct = createConstraintMatrix(obj)
            s.dirDOFs        = obj.bc.dirichlet;
            s.sizeK          = obj.sizeK;
            DirComputer      = DirichletComputer(s);
            [CtDir, obj.sizeDir] = DirComputer.compute();
            %reformulation starts here
            perDOFmaster = obj.bc.periodic_free;
            perDOFslave  = obj.bc.periodic_constrained;
            nEqperType = obj.sizePer/4;
            obj.nConstraints = obj.sizeDir + nEqperType*3; 
            Ct = zeros(nEqperType*3 + obj.sizeDir, obj.sizeK);
%             Ct = zeros(nEqperType*3, obj.sizeK);
            pX = zeros(nEqperType, 1);
            pY = zeros(nEqperType, 1);
            pXY = zeros(nEqperType, 1); 
            % first 39
            for i = 1:nEqperType
                masterDOF = perDOFmaster(i);
                slaveDOF = perDOFslave(i);
                Ct(i, [masterDOF slaveDOF]) = [1 -1];
                pX(i) = i;
            end
            % second 39
            for i = nEqperType+1:2*nEqperType
                masterDOF = perDOFmaster(i);
                slaveDOF = perDOFslave(i);
                pXY(i-nEqperType) = i;
                Ct(i, [masterDOF slaveDOF]) = [1 -1];
            end
            % third 39
            for i = 2*nEqperType+1:3*nEqperType
                masterDOF = perDOFmaster(i);
                slaveDOF = perDOFslave(i);
                positionCt = i-nEqperType;
                Ct(positionCt, [masterDOF slaveDOF]) = [1 -1];
            end
            %last 39
            for i = 3*nEqperType+1:4*nEqperType
                masterDOF = perDOFmaster(i);
                slaveDOF = perDOFslave(i);
                positionCt = i-nEqperType;
                pY(i-3*nEqperType) = positionCt;
                Ct(positionCt, [masterDOF slaveDOF]) = [1 -1];
            end
            Ct(nEqperType*3+1:end, :) = CtDir;
        end

        function [L, stressHomog] = computeStressHomog(obj, sol)
            % reformulation starts here
            nEqperType = obj.sizePer/4;
            Lx  = 0;
            Ly  = 0;
            Lxy = 0;
            d1  = obj.sizeK+1;
            d2  = obj.sizeK + nEqperType*3;
            L   = sol(d1:d2);
            for i = 1:nEqperType
                Lx = Lx + L(i);
            end
            % second 39
            for i = nEqperType+1:2*nEqperType
                Lxy = Lxy + L(i);
            end
            % third 39
            for i = 2*nEqperType+1:3*nEqperType
                Ly = Ly + L(i);
            end
            stressHomog = -[Lx; Ly; Lxy];
        end
    end

end