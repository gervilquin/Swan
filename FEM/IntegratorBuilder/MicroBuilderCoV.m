classdef MicroBuilderCoV < handle

    properties (Access = private)
        LHS
        bc
        sizeK
        nConstraints
        solver
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
            fluc = obj.computeFluctuations(u);
        end

    end

    methods (Access = private)
        function init(obj, cParams)
            obj.LHS             = cParams.LHS;
            obj.bc              = cParams.bc;
            obj.sizeK           = size(cParams.LHS, 1);
            obj.solver          = cParams.solver;
            obj.mesh            = cParams.mesh;
            obj.vstrain         = cParams.vstrain;
            perDOFslave         = obj.bc.periodic_constrained;
            obj.sizePer         = size(perDOFslave, 1);
            s.dirDOFs           = obj.bc.dirichlet;
            s.sizeK             = obj.sizeK;
            s.vstrain           = obj.vstrain;
            s.mesh              = obj.mesh;
            obj.ConditionSetter = MicroDirichletSetter(s);
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
            nU           = size(obj.LHS, 1);
            der1Rhs      = zeros(nU, 1);
            nEqperType   = obj.sizePer/4;
            perVector    = zeros(nEqperType*3, 1);
            for i = 1:nEqperType
                perVector(i, 1) = -obj.vstrain(1);
            end
            for i = nEqperType+1:2*nEqperType
                perVector(i, 1) = -obj.vstrain(3);
            end
            for i = 2*nEqperType+1:3*nEqperType
                perVector(i, 1) = -obj.vstrain(2);
            end
            [RHSDir, RHSDirPer] = obj.ConditionSetter.setDirichletRhs();
            fullRHS             = [der1Rhs; perVector; RHSDirPer; RHSDir];
        end

        function Ct = createConstraintMatrix(obj)
            [CtDir, CtPerDir] = obj.ConditionSetter.setDirichletLhs();
            perDOFmaster      = obj.bc.periodic_free;
            perDOFslave       = obj.bc.periodic_constrained;
            nEqperType        = obj.sizePer/4;
            obj.nConstraints  = size(CtDir, 1)+ size(CtPerDir, 1) ...
                + nEqperType*3; 
            Ct                = zeros(nEqperType*3+size(CtDir, 1)+ ...
                size(CtPerDir, 1), obj.sizeK);
            for i = 1:nEqperType
                masterDOF                   = perDOFmaster(i);
                slaveDOF                    = perDOFslave(i);
                Ct(i, [masterDOF slaveDOF]) = [1 -1];
            end
            for i = nEqperType+1:2*nEqperType
                masterDOF                   = perDOFmaster(i);
                slaveDOF                    = perDOFslave(i);
                Ct(i, [masterDOF slaveDOF]) = [1 -1];
            end
            for i = 2*nEqperType+1:3*nEqperType
                masterDOF                            = perDOFmaster(i);
                slaveDOF                             = perDOFslave(i);
                positionCt                           = i-nEqperType;
                Ct(positionCt, [masterDOF slaveDOF]) = [1 -1];
            end
            for i = 3*nEqperType+1:4*nEqperType
                masterDOF                            = perDOFmaster(i);
                slaveDOF                             = perDOFslave(i);
                positionCt                           = i-nEqperType;
                Ct(positionCt, [masterDOF slaveDOF]) = [1 -1];
            end
            Ct(nEqperType*3+1:end, :) = [CtPerDir; CtDir];
        end

        function [L, LDir, stressHomog] = computeStressHomog(obj, sol)
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
            for i = nEqperType+1:2*nEqperType
                Lxy = Lxy + L(i);
            end
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
                Lxy = Lxy + LDir(1) + LDir(2);
            end
            stressHomog = [Lx; Ly; Lxy];
        end

        function fluct = computeFluctuations(obj, u)
            coords  = obj.mesh.coord';
            nel     = size(coords, 2);
            strainM = [obj.vstrain(1) obj.vstrain(3)/2; 
                        obj.vstrain(3)/2 obj.vstrain(2)];
            uTotal  = zeros(obj.sizeK, 1);
            for i=1:nel
                if i == 208
                    coordP                 = coords(:, i);
                    uTotal([2*i-1 2*i], 1) = strainM*coordP;
                end
                coordP = coords(:, i);
                uTotal([2*i-1 2*i], 1) = strainM*coordP;
            end
            fluct = u - uTotal;
        end
    end

end