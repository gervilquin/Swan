classdef MicroBuilder < handle
% solves micro structure obtaining fluctuations as a result

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
        C
        perDirConst
    end

    methods (Access = public)
        function createBuilder(obj, s)
            obj.init(s);
        end

        function [u, stressHomog] = solveSystem(obj)
            %Meter aqui el defRHS i el defLHS y hacerlos privados.
            defLHS = obj.createLHS();
            defRHS = obj.createRHS();
            sol = obj.solver.solve(defLHS, defRHS);
            u   = sol(1:obj.sizeK, 1);
            [L, stressHomog] = obj.computeStressHomog(sol);
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
%             nPerDofs = size(obj.bc.periodic_constrained, 1) + obj.perDirConst;
            perVector = zeros(nPerDofs, 1);
%             perVector = zeros(3,1);
            uD = obj.bc.dirichlet_values;
            fullRHS = [obj.RHS; perVector; uD]; 
        end

        function Ct = createConstraintMatrix(obj)
            %% OLD VERSION
            s.dirDOFs        = obj.bc.dirichlet;
            s.sizeK          = obj.sizeK;
            DirComputer      = DirichletComputer(s);
            [CtDir, sizeDir] = DirComputer.compute();
%                     [CtDir, sizeDir] = obj.computeDirichletCond();
            perDOFslave      = obj.bc.periodic_constrained;
            perDOFmaster     = obj.bc.periodic_free;
            obj.perDirConst = sizeDir/2;
            obj.nConstraints = sizeDir + obj.sizePer; %+ obj.perDirConst;
%             obj.nConstraints = obj.sizePer; 
            Ct               = zeros(obj.nConstraints, obj.sizeK);
            for i = 1:obj.sizePer
                masterNode = perDOFmaster(i);
                slaveNode = perDOFslave(i);
                Ct(i, [masterNode slaveNode]) = [1 -1];
            end
            % Periodic conditions Dirichlet
            %compute master-slave
%             masterPerDir = obj.bc.dirichlet(1:obj.perDirConst);
%             slavePerDir = obj.bc.dirichlet(obj.perDirConst+1:end);
%             for i = 1:obj.perDirConst
%                 masterNode = masterPerDir(i);
%                 slaveNode = slavePerDir(i);
%                 Ct(obj.sizePer+i, [masterNode slaveNode]) = [1 -1];
%             end
            %
%             Ct(obj.sizePer+obj.perDirConst+1:end, :) = CtDir;
            Ct(obj.sizePer+1:end, :) = CtDir;
            obj.C = Ct;

            %% NEW VERSION WITH MINIMUM CONDITIONS
%             s.dirDOFs        = obj.bc.dirichlet;
%             s.sizeK          = obj.sizeK;
%             DirComputer      = DirichletComputer(s);
%             [CtDir, sizeDir] = DirComputer.compute();
%             obj.nConstraints = sizeDir + 3;
%             Ct               = zeros(obj.nConstraints, obj.sizeK);
%             nEqperType       = obj.sizePer/4;
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
%             obj.C = Ct;
        end

%         function [CtDir, sizeDir] = computeDirichletCond(obj)
%             dirDOFs    = obj.bc.dirichlet;
%             sizeDir = size(dirDOFs, 1);
%             CtDir = zeros(sizeDir, obj.sizeK);
%             for i = 1:sizeDir
%                 CtDir(i,dirDOFs(i)) = 1; 
%             end
%         end

        function [L, stressHomog] = computeStressHomog(obj, sol)
            %% OLD VERSION POINT PER POINT
            nEqperType = obj.sizePer/4;
            sigmaX = 0;
            sigmaY = 0;
            tauXY  = 0;
            d1  = obj.sizeK+1;
            d2  = obj.sizeK + obj.sizePer;
            d3  = obj.sizeK + obj.sizePer + obj.perDirConst;
            L   = sol(d1:d2);
            LDir = sol(d2+1:end);
%             LPerDir = sol(d2+1:d3); 
            % first 8: xx
            for i = 1:nEqperType
                sigmaX = sigmaX + L(i);
%                  if i == 1
%                     sigXaux = L(i);
%                     sigmaX = sigmaX + sigXaux;
%                 end
            end
            % second 8: xy
            for i = nEqperType+1:2*nEqperType
                tauXY = tauXY + L(i);
%                 if i == nEqperType+1
%                     sigXYaux = L(i);
%                     tauXY = tauXY + sigXYaux;
%                 end
            end
            % third 8: xy
            for i = 2*nEqperType+1:3*nEqperType
                tauXY = tauXY + L(i);
%                 if i ==3*nEqperType
%                     sigXYaux = L(i);
%                     tauXY = tauXY + sigXYaux;
%                 end
            end
            % last 8: yy
            for i = 3*nEqperType+1:4*nEqperType
                sigmaY = sigmaY + L(i);
%                 if i == 3*nEqperType+1
%                     sigYaux = L(i);
%                     sigmaY = sigmaY + sigYaux;
%                 end
            end
            sigmaX = sigmaX + LDir(1) + LDir(5);
            sigmaY = sigmaY + LDir(2) + LDir(4);
            tauXY = (tauXY + LDir(6) + LDir(3) - LDir(7) - LDir(8))/2;
            stressHomog = [sigmaX; sigmaY; tauXY];

            %% NEW VERSION WITH MINIMUM CONDITIONS
%             d1  = obj.sizeK+1;
%             L = sol(d1:end);
%             nEqperType = obj.sizePer/4;
%             stressHomog = [L(1)*(nEqperType+1); L(3)*(nEqperType+1); L(2)*(nEqperType+1)];
%             stressHomog = [L(1)*9; L(3)*9; L(2)*9];
        end

    end
end