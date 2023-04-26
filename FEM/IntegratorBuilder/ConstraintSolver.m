classdef ConstraintSolver < handle 

    properties (Access = public)

    end

    properties (Access = private)
        solver
        bc
        K
        RHS
        sizeK
        vstrain
        sizePer
        nConstraints
        solType
        solMode
        scale
    end

    methods (Access = public)
        function obj = ConstraintSolver(cParams)
            obj.init(cParams);
        end

        function fullRHS = assembleGlobalRHS(obj)
            switch obj.solType 
                case 'MONOLITIC'
                    fullRHS = obj.createGeneralVector();
                case 'REDUCED'
                    R           = obj.computeReactions();
                    CompleteRHS = obj.RHS + R;
                    fullRHS     = obj.fullToReducedVector(CompleteRHS);
            end
        end

        function fullLHS = assembleGlobalLHS(obj)
            switch obj.solType 
                case 'MONOLITIC'
                    Ct      = obj.createConstraintMatrix();
                    fullLHS = obj.createGeneralMatrix(Ct);
                case 'REDUCED'
                    fullLHS = fullToReducedMatrix(obj);
            end
        end

        function [u, L] = solveSystem(obj, lhs, rhs)
            sol         = obj.solver.solve(lhs, rhs);
            switch obj.solType
                case 'REDUCED'
                    u = reducedToFullVector(obj, sol);
                    if strcmp(obj.scale, 'MACRO')
                        L = obj.getReactions(sol);
                    else
                        L = 0;
                    end
                case 'MONOLITIC'
                    u = sol(1:obj.sizeK, 1);
                    L = sol(obj.sizeK+1:end, 1);
            end
         
        end
    end

    methods (Access = private)
        function init(obj, cParams)
            obj.solver = cParams.solver;
            obj.bc     = cParams.bc; 
            obj.K      = cParams.LHS;
            obj.sizeK  = size(obj.K, 1);
            obj.RHS    = cParams.RHS;
            obj.scale  = cParams.scale;
            if isfield(cParams, 'vstrain')
                obj.vstrain = cParams.vstrain;
                perDOFslave = obj.bc.periodic_constrained;
                obj.sizePer = size(perDOFslave, 1);
            end 
            obj.solType = cParams.solType;
            obj.solMode = cParams.solMode;
        end

        function fullLHS = createGeneralMatrix(obj, Ct)
            C       = Ct';
            nC      = obj.nConstraints;
            Z       = zeros(nC);
            Km      = obj.K;
            fullLHS = [Km C; C' Z];
        end

        function Ct = createConstraintMatrix(obj)
            switch obj.solMode
                case 'DISP'
                    if ~isempty(obj.vstrain)
                        [CtDir, CtPerDir] = obj.bc.conditionSetter.setDirichletLhs();
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
                    else
                        s.dirDOFs        = obj.bc.dirichlet;
                        s.sizeK          = obj.sizeK;
                        DirComputer      = DirichletComputer(s);
                        [CtDir, sizeDir] = DirComputer.computeDirCond;
                        obj.nConstraints = sizeDir;
                        Ct = CtDir;
                    end

                case 'FLUC'
                    if isprop(obj, 'vstrain')
                        s.dirDOFs        = obj.bc.dirichlet;
                        s.sizeK          = obj.sizeK;
                        DirComputer      = DirichletComputer(s);
                        [CtDir, sizeDir] = DirComputer.computeDirCond();
                        perDOFslave      = obj.bc.periodic_constrained;
                        perDOFmaster     = obj.bc.periodic_free;
                        obj.nConstraints = sizeDir + obj.sizePer; 
                        Ct               = zeros(obj.nConstraints, obj.sizeK);
                        for i = 1:obj.sizePer
                            masterNode                    = perDOFmaster(i);
                            slaveNode                     = perDOFslave(i);
                            Ct(i, [masterNode slaveNode]) = [1 -1];
                        end
                        Ct(obj.sizePer+1:end, :) = CtDir;
                    else
                        s.dirDOFs        = obj.bc.dirichlet;
                        s.sizeK          = obj.sizeK;
                        DirComputer      = DirichletComputer(s);
                        [CtDir, sizeDir] = DirComputer.computeDirCond();
                        obj.nConstraints = sizeDir;
                        Ct = CtDir;
                    end

            end
            
        end

        function fullRHS = createGeneralVector(obj) %Disp

            switch obj.solMode
                case 'DISP'
                    if ~isempty(obj.vstrain)
                        nU           = size(obj.K, 1);
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
                        [RHSDir, RHSDirPer] = obj.bc.conditionSetter.setDirichletRhs();
                        fullRHS             = [der1Rhs; perVector; RHSDirPer; RHSDir];
                    else
                        uD        = obj.bc.dirichlet_values;
                        fullRHS   = [obj.RHS; uD];
                    end

                case 'FLUC'
                    if isprop(obj, 'vstrain')
                        nPerDofs  = size(obj.bc.periodic_constrained, 1);
                        perVector = zeros(nPerDofs, 1);
                        uD        = obj.bc.dirichlet_values;
                        fullRHS   = [obj.RHS; perVector; uD]; 
                    else
                        uD        = obj.bc.dirichlet_values;
                        fullRHS   = [obj.RHS; uD];
                    end
            end
        end

        function red = fullToReducedMatrix(obj)
            switch obj.scale
                case 'MACRO'
                    red = obj.reduceMatrixDirichlet(obj.K);
                case 'MICRO'
                    red = obj.reduceMatrixPeriodic(obj.K);
            end
        end

        function red = fullToReducedVector(obj, vec)
            switch obj.scale
                case 'MACRO'
                    red = obj.reduceVectorDirichlet(vec);
                case 'MICRO'
                    red = obj.reduceVectorPeriodic(vec);
            end
        end
        
        function full = reducedToFullVector(obj, vec)
            switch obj.scale
                case 'MACRO'
                    full = obj.expandVectorDirichlet(vec);
                case 'MICRO'
                    full = obj.expandVectorPeriodic(vec);
            end
        end

        function Ared = reduceMatrixDirichlet(obj,A)
            fr = obj.bc.free';
            Ared = A(fr,fr);
        end
        
        function b_red = reduceVectorDirichlet(obj,b)
            fr = obj.bc.free';
            b_red = b(fr);
        end

        function Ared = reduceMatrixPeriodic(obj,A)
            MS = obj.bc.masterSlave;
            vF = obj.bc.free;
            vP = obj.computePeriodicNodes(MS(:,1));
            vQ = obj.computePeriodicNodes(MS(:,2));
            vI = setdiff(vF,vP);
            
            A_II = A(vI,vI);
            A_IP = A(vI,vP) + A(vI,vQ);
            A_PI = A(vP,vI) + A(vQ,vI); 
            A_PP = A(vP,vP) + A(vP,vQ) + A(vQ,vP) + A(vQ,vQ); 

            Ared = [A_II, A_IP; A_PI, A_PP];
        end
        
        function b_red = reduceVectorPeriodic(obj,b)
            vF = obj.bc.free;
            vP = obj.bc.periodic_free;
            vQ = obj.bc.periodic_constrained;
            vI = setdiff(vF,vP);
            b_I = b(vI);
            b_P = b(vP) + b(vQ);
            b_red = [b_I; b_P];
        end

        function b = expandVectorDirichlet(obj,bfree)
            dir = obj.bc.dirichlet;
            uD  = obj.bc.dirichlet_values;
            fr  = obj.bc.free;
            nsteps = length(bfree(1,:));
            ndof = sum(obj.bc.ndofs);
            uD = repmat(uD,1,nsteps);
            
            b = zeros(ndof,nsteps);
            b(fr,:) = bfree;
            if ~isempty(dir)
                b(dir,:) = uD;
            end
        end

        function b = expandVectorPeriodic(obj,bfree)
            vF = obj.bc.free;
            vP = obj.bc.periodic_free;
            vC = obj.bc.periodic_constrained;
            vI = setdiff(vF,vP);
            b = zeros(obj.bc.ndofs,1);
            b(vI) = bfree(1:1:size(vI,2));
            b(vP) = bfree(size(vI,2)+1:1:size(bfree,1));
            b(vC) = b(vP);
        end

        function perDof = computePeriodicNodes(obj,perNodes)
            nunkn   = obj.bc.dim.ndimf;
            nlib    = size(perNodes,1);
            perDof  = zeros(nlib*nunkn,1);
            for iunkn = 1:nunkn
                indDof           = nlib*(iunkn - 1) + [1:nlib];
                perDof(indDof,1) = obj.nod2dof(obj.bc.dim.ndimf, perNodes,iunkn);
            end
        end

        function idof = nod2dof(obj, ndimf, inode, iunkn)
            idof(:,1) = ndimf*(inode - 1) + iunkn;
        end

        function R = computeReactions(obj)
            boundaryCond  = obj.bc;
            dirich        = boundaryCond.dirichlet;
            dirichV       = boundaryCond.dirichlet_values;
            if ~isempty(dirich)
                R = -obj.K(:,dirich)*dirichV;
            else
                R = zeros(sum(obj.bc.dim.ndofs(:)),1);
            end
        end

        function R = getReactions(obj, sol)
            R         = zeros(obj.sizeK, 1);
            dirich    = obj.bc.dirichlet;
            dirichV   = obj.bc.dirichlet_values;
            free      = obj.bc.free;
            dl        = sol;
            R(dirich) = obj.K(dirich, free)*dl + obj.K(dirich, dirich)*dirichV;
            R         = R(dirich);
        end

    end


end