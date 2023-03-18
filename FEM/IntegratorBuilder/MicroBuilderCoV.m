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

        %% OLD VERSION: POINT PER POINT 3 EQUATION STYLE. 
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
            Ct(nEqperType*3+1:end, :) = CtDir;
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
            stressHomog = [Lx; Ly; Lxy];
        end

%% NEW VERSION: POINT PER POINT WITHOUT IMPLEMENTING 3RD EQ.
% 
%         function Ct = createConstraintMatrix(obj)
%             s.dirDOFs        = obj.bc.dirichlet;
%             s.sizeK          = obj.sizeK;
%             DirComputer      = DirichletComputer(s);
%             [CtDir, obj.sizeDir] = DirComputer.compute();
%             perDOFslave      = obj.bc.periodic_constrained;
%             perDOFmaster     = obj.bc.periodic_free;
%             obj.nConstraints = obj.sizeDir + obj.sizePer; 
%             Ct               = zeros(obj.nConstraints, obj.sizeK);
%             for i = 1:obj.sizePer
%                 masterNode = perDOFmaster(i);
%                 slaveNode = perDOFslave(i);
%                 Ct(i, [masterNode slaveNode]) = [1 -1];
%             end
%             Ct(obj.sizePer+1:end, :) = CtDir;
%         end
% 
% %         function [CtDir, sizeDir] = computeDirichletCond(obj)
% %             dirDOFs    = obj.bc.dirichlet;
% %             sizeDir = size(dirDOFs, 1);
% %             CtDir = zeros(sizeDir, obj.sizeK);
% %             for i = 1:sizeDir
% %                 CtDir(i,dirDOFs(i)) = 1; 
% %             end
% %         end
% 
%         function [L, stressHomog] = computeStressHomog(obj, sol)
%             nEqperType = obj.sizePer/4;
%             sigmaX = 0;
%             sigmaY = 0;
%             tauXY  = 0;
%             d1  = obj.sizeK+1;
%             d2  = obj.sizeK + obj.sizePer;
%             L   = sol(d1:d2);
%             LDir = sol(d2+1:end);
%             % first 8: xx
%             for i = 1:nEqperType
%                 sigmaX = sigmaX + L(i);
%             end
%             % second 8: xy
%             for i = nEqperType+1:2*nEqperType
%                 tauXY = tauXY + L(i);
%             end
%             % third 8: xy
%             for i = 2*nEqperType+1:3*nEqperType
%                 tauXY = tauXY + L(i);
%             end
%             % last 8: yy
%             for i = 3*nEqperType+1:4*nEqperType
%                 sigmaY = sigmaY + L(i);
%             end
%             stressHomog = -[sigmaX; sigmaY; tauXY/2];
%         end
% 
%         function fullRHS = createGeneralVector(obj)
%             nU        = size(obj.LHS, 1);
%             der1Rhs   = zeros(nU, 1);
%             uD        = obj.bc.dirichlet_values;
%             nEqperType = obj.sizePer/4;
%             perVector    = zeros(nEqperType*4, 1);
%             for i = 1:nEqperType
%                 perVector(i, 1) = obj.vstrain(1);
%             end
%             % second 39
%             for i = nEqperType+1:2*nEqperType
%                 perVector(i, 1) = obj.vstrain(3);
%             end
%             % third 39
%             for i = 2*nEqperType+1:3*nEqperType
%                 perVector(i, 1) = obj.vstrain(3);
%             end
%             % last 39
%             for i = 3*nEqperType+1:4*nEqperType
%                 perVector(i, 1) = obj.vstrain(2);
%             end
%             fullRHS   = [der1Rhs; perVector; uD];
%         end

%% NEW VERSION WITH MINIMUM CONDITIONS
%         function Ct = createConstraintMatrix(obj)
%             s.dirDOFs        = obj.bc.dirichlet;
%             s.sizeK          = obj.sizeK;
%             DirComputer      = DirichletComputer(s);
%             [CtDir, obj.sizeDir] = DirComputer.compute();
%             obj.nConstraints = obj.sizeDir + 3;
%             Ct               = zeros(obj.nConstraints, obj.sizeK);
%             nEqperType = obj.sizePer/4;
%             perDOFslave      = obj.bc.periodic_constrained;
%             perDOFmaster     = obj.bc.periodic_free;
%             for i = 1:nEqperType
%                 masterNode = perDOFmaster(i);
%                 slaveNode = perDOFslave(i);
%                 Ct(1, [masterNode slaveNode]) = [1 -1];
%             end
%             % second 8: xy
%             for i = nEqperType+1:2*nEqperType
%                 masterNode = perDOFmaster(i);
%                 slaveNode = perDOFslave(i);
%                 Ct(2, [masterNode slaveNode]) = [1 -1];                
%             end
%             % third 8: xy
%             for i = 2*nEqperType+1:3*nEqperType
%                 masterNode = perDOFmaster(i);
%                 slaveNode = perDOFslave(i);
%                 Ct(2, [masterNode slaveNode]) = [1 -1];               
%             end
%             % last 8: yy
%             for i = 3*nEqperType+1:4*nEqperType
%                 masterNode = perDOFmaster(i);
%                 slaveNode = perDOFslave(i);
%                 Ct(3, [masterNode slaveNode]) = [1 -1];               
%             end
%             Ct(4:end, :) = CtDir;
%         end
% 
%         function fullRHS = createGeneralVector(obj)
%             nPerDofs = size(obj.bc.periodic_constrained, 1);
% %             perVector = zeros(nPerDofs, 1);
%             perVector = [obj.vstrain(1); obj.vstrain(3); obj.vstrain(2)]; 
%             uD = obj.bc.dirichlet_values;
%             nU        = size(obj.LHS, 1);
%             der1Rhs   = zeros(nU, 1);
%             fullRHS = [der1Rhs; perVector; uD]; 
%         end
% 
%         function [L, stressHomog] = computeStressHomog(obj, sol)
%             d1  = obj.sizeK+1;
%             nEqperType = obj.sizePer/4;
%             L = sol(d1:end);
%             stressHomog = [L(1)*(nEqperType+1); L(3)*(nEqperType+1); L(2)*(nEqperType+1)];
%         end
    end

end