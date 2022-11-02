classdef H1Projector_toPeriodicP1 < Projector

    properties (Access = private)
        fieldMass
        fieldStiffness
        unfMesh
        masterSlave
    end
    
    methods (Access = public)

        function obj = H1Projector_toPeriodicP1(cParams)
            obj.init(cParams);
            obj.masterSlave = cParams.masterSlave;
            if isprop(obj.mesh,'unfittedBoundaryMesh')
                obj.unfMesh = obj.mesh;
                obj.mesh    = obj.mesh.backgroundMesh;
            end
            obj.createFieldMass();
            obj.createFieldStiffness();
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = obj.solveSystem(LHS,RHS);
            s.type    = obj.mesh.type;
            s.connec  = obj.mesh.connec;
            s.fValues = xProj;
            xFun = P1Function(s);
        end

    end

    methods (Access = private)

        function createFieldMass(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1; 
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            obj.fieldMass = Field(s);
        end

        function createFieldStiffness(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1; 
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'CONSTANT';
            obj.fieldStiffness = Field(s);
        end
        
        function LHS = computeLHS(obj)
            LHSM    = obj.computeMassMatrix();
            LHSK    = obj.computeStiffnessMatrix();
            epsilon = obj.mesh.computeMeanCellSize();
            LHS     = LHSM + epsilon^2*LHSK;
        end

        function LHSM = computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.fieldMass;
            lhs = LHSintegrator.create(s);
            LHSM = lhs.compute();
        end

        function LHSK = computeStiffnessMatrix(obj)
            s.type  = 'StiffnessMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.fieldStiffness;
            lhs = LHSintegrator.create(s);
            LHSK = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
            if isempty(obj.unfMesh)
                quad = obj.createRHSQuadrature(fun);
                xV = quad.posgp;
                dV = obj.mesh.computeDvolume(quad);
                obj.mesh.interpolation.computeShapeDeriv(xV);
                shapes = permute(obj.mesh.interpolation.shape,[1 3 2]);
                conne = obj.mesh.connec;

                nGaus = quad.ngaus;
                nFlds = fun.ndimf;
                nNode = size(conne,2);
                nDofs = obj.mesh.nnodes;

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
            else
                npnod = obj.unfMesh.backgroundMesh.nnodes;
                f = ones(npnod,1);
                s.mesh = obj.unfMesh;
                s.type = 'Unfitted';
                integrator = RHSintegrator.create(s);
                RHS = integrator.integrateInDomain(f);
            end
        end

        function q = createRHSQuadrature(obj, fun)
            ord = obj.determineQuadratureOrder(fun);
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end

        function x = solveSystem(obj,LHS,RHS)
            x = zeros(size(RHS));
            dofs = 1:1:length(x);
            master = obj.masterSlave(:,1);
            slave  = obj.masterSlave(:,2);
            n      = length(master);
            for i = 1:n
                LHS(:,master(i)) = LHS(:,master(i)) + LHS(:,slave(i));
                LHS(master(i),:) = LHS(master(i),:) + LHS(slave(i),:);
                RHS(master(i))   = RHS(master(i))   + RHS(slave(i));
            end
            LHS(:,slave) = [];
            LHS(slave,:) = [];
            RHS(slave)   = [];
            xPeriodic = LHS\RHS;
            for i =1:length(slave)
                dofs(dofs==slave(i)) = [];
            end
            x(dofs) = xPeriodic;
            x(slave) = x(master);
        end

    end

end