classdef MicroBuilderCoV < handle

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
        ConditionSetter
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
            [L, LDir, stressHomog] = obj.computeStressHomog(sol);
            uTotal = obj.computeTotalDisplacements(u);
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
            s.dirDOFs        = obj.bc.dirichlet;
            s.sizeK          = obj.sizeK;
            s.vstrain        = obj.vstrain;
            s.mesh           = obj.mesh;
            obj.ConditionSetter = GlobalDofSetter(s);
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
                perVector(i, 1) = -obj.vstrain(1);
            end
            % second 39
            for i = nEqperType+1:2*nEqperType
                perVector(i, 1) = -obj.vstrain(3);
            end
            % third 39
            for i = 2*nEqperType+1:3*nEqperType
                perVector(i, 1) = -obj.vstrain(2);
            end
%             fullRHS   = [der1Rhs; perVector; uD];
            % implement new GlobalDofSetter class
            [RHSDir, RHSDirPer] = obj.ConditionSetter.setDirichletRhs();
            fullRHS   = [der1Rhs; perVector; RHSDirPer; RHSDir];

%             PerDirRhs = [-1; -1];
%             DirRhs = zeros(6,1);
%             fullRHS   = [der1Rhs; perVector; PerDirRhs; DirRhs];
        end

        function Ct = createConstraintMatrix(obj)
%             s.dirDOFs        = obj.bc.dirichlet;
%             s.sizeK          = obj.sizeK;
%             s.vstrain        = obj.vstrain;
%             s.mesh           = obj.mesh;
%             DirComputer      = DirichletComputer(s);
            % implement new GlobalDofSetter class
            [CtDir, CtPerDir] = obj.ConditionSetter.setDirichletLhs();
%             [CtDir, obj.sizeDir] = DirComputer.compute();
            %only fix nodes on left edges
%             CtDir = zeros(6, obj.sizeK);
%             CtDir(1,1) = 1;
%             CtDir(2,2) = 1;
%             CtDir(3,409) = 1;
%             CtDir(4,410) = 1;
%             CtDir(5,408) = 1;
%             CtDir(6,530) = 1;
%             CtPerDir = zeros(2, obj.sizeK);
%             CtPerDir(1,[1, 407]) = [1 -1];
%             CtPerDir(2,[409, 529]) = [1 -1];
            %reformulation starts here
            perDOFmaster = obj.bc.periodic_free;
            perDOFslave  = obj.bc.periodic_constrained;
            nEqperType = obj.sizePer/4;
%             obj.nConstraints = obj.sizeDir + nEqperType*3; 
            obj.nConstraints = size(CtDir, 1)+ size(CtPerDir, 1) + nEqperType*3; 
%             Ct = zeros(nEqperType*3 + obj.sizeDir, obj.sizeK);
            Ct = zeros(nEqperType*3+size(CtDir, 1)+ size(CtPerDir, 1), obj.sizeK);
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
%                 positionCt = 4*nEqperType+1-i;
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
            Ct(nEqperType*3+1:end, :) = [CtPerDir; CtDir];
        end

        function [L, LDir, stressHomog] = computeStressHomog(obj, sol)
            % reformulation starts here
            nEqperType = obj.sizePer/4;
            Lx  = 0;
            Ly  = 0;
            Lxy = 0;
            d1  = obj.sizeK+1;
            d2  = obj.sizeK + nEqperType*3;
            L   = sol(d1:d2);
            LDir = sol(d2+1:end);
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

            if obj.vstrain(1) == 1
                Lx = Lx + LDir(1) + LDir(2) + LDir(3) + LDir(5);
                Ly = Ly + LDir(4) + LDir(7);
            elseif obj.vstrain(2) == 1
                Ly = Ly + LDir(1) + LDir(2) + LDir(4) + LDir(6);
                Lx = Lx + LDir(3) + LDir(7);                
            else
            end
            stressHomog = [Lx; Ly; Lxy];
        end

        function fluct = computeTotalDisplacements(obj, u)
            coords = obj.mesh.coord';
            nel = size(coords, 2);
            strainM = [obj.vstrain(1) obj.vstrain(3); 
                obj.vstrain(3) obj.vstrain(2)];
            uTotal = zeros(obj.sizeK, 1);
            for i=1:nel
                uTotal([2*i-1 2*i], 1) = strainM*coords(:, i);
            end
            fluct = u - uTotal;
        end
    end

end