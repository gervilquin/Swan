classdef FunSolver < handle

    % Analytical Solver (AX = b)

    properties (Access = private)
        BC
        mesh
    end

    properties (Access = private)
        free
    end

    methods (Access = public)

        function obj = FunSolver(cParams)
            obj.init(cParams);
        end

        function x = solve(obj,LHS,RHS)
            fr = obj.computeFreeDofs(LHS);
            lhs = obj.reduceMatrix(LHS,fr);
            rhs = obj.reduceVector(RHS,fr);
            x = lhs\rhs;
%             X = obj.expandVector(x,fr);
        end

    end

    methods (Access = private)

        function obj = init(obj,cParams)
            obj.BC   = cParams.BC;
            obj.mesh = cParams.mesh;
        end

        function a = reduceMatrix(obj,A,fr)
            a  = A(fr,fr);
        end
        
        function b = reduceVector(obj,B,fr)
            b = B(fr);
        end

        function b = expandVectorDirichlet(obj,bfree)
            dir = obj.dirichlet;
            uD  = obj.dirichlet_values;
            nsteps = length(bfree(1,:));
            ndof = sum(obj.ndofs);
            uD = repmat(uD,1,nsteps);
            
            b = zeros(ndof,nsteps);
            b(fr,:) = bfree;
            if ~isempty(dir)
                b(dir,:) = uD;
            end
        end

        function free = computeFreeDofs(obj, A)
            nNods = size(obj.mesh.coord,1);
            nDimf = size(A,1)/nNods;
            nDofs = nNods*nDimf;
            dofs = 1:1:nDofs;
            constr = 1:1:nDofs;
            coor = obj.mesh.coord';
            isDirichlet = obj.BC.dirichlet.domainFun(coor)';
            freeNod = find(~isDirichlet);
            for iDimf = 1:nDimf
                iDof = obj.nod2dof(nDimf, freeNod, iDimf);
                constr = setdiff(constr,iDof);
            end
            free = setdiff(dofs,constr);

        end

        function idof = nod2dof(obj, ndimf, inode, iunkn)
%             ndimf = obj.dim.ndimf;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

    end

end