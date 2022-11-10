classdef Projector_toP2 < Projector

    properties (Access = private)
        field
    end
    
    methods (Access = public)

        function obj = Projector_toP2(cParams)
            obj.init(cParams);
            obj.createField();
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.type    = obj.mesh.type;
            s.connec  = obj.field.connec;
            s.fValues = xProj;
            xFun = P2Function(s);
        end

    end

    methods (Access = private)

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1; 
            s.interpolationOrder = 'QUADRATIC';
            s.quadratureOrder    = 'QUADRATIC';
            obj.field = Field(s);
        end
        
        function LHS = computeLHS(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.field;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);
            obj.field.interpolation.computeShapeDeriv(xV);
            shapes = permute(obj.field.interpolation.shape,[1 3 2]);
            conne = obj.field.connec;

            nGaus = quad.ngaus;
            nFlds = fun.ndimf;
            nNode = size(conne,2);
            nDofs = obj.field.dim.nnodes;

            fGaus = fun.evaluate(xV);
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(iField,igaus,:));
                    for inode = 1:nNode
                        dofs = conne(:,inode);
                        Ni = shapes(inode,igaus);
                        int = Ni*fG.*dVg;
                        f(:,iField) = f(:,iField) + accumarray(dofs,int,[nDofs 1]);
                    end
                end
            end
            RHS = f;
        end

        function q = createRHSQuadrature(obj, fun)
%             ord = obj.determineQuadratureOrder(fun);
            ord = 'QUADRATIC'; % CUBIC/ QUADRATIC should be enough
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end
        
    end

end