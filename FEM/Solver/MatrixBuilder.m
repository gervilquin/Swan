classdef MatrixBuilder < handle
    % Class to build matrix for direct solution without reductions. 

    methods (Static)
        function genM = createGeneralMatrix(K, bc)
            % Method to create the general matrix
            sK = size(K, 1);
            Ct = MatrixBuilder.buildCMatrix(bc, sK);
            C = Ct';
            sC1 = size(C, 1);
            sC2 = size(C, 2);
            genM = zeros(sK+sC2);
            genM(1:sK, 1:sK) = K;
            genM(sK+1:end, 1:sC1) = C';
            genM(1:sC1, sK+1:end) = C;
        end

        function genFVector = createGeneralVector(RHS, bc)
            dirichletValues = bc.dirichlet_values;
            sRHS = size(RHS,1);
            sDirichlet = size(dirichletValues, 1);
            genFVector = zeros(sRHS+sDirichlet, 1);
            genFVector(1:sRHS) = RHS;
            genFVector(sRHS+1:end) = dirichletValues;
        end

        function Ct = buildCMatrix(bc, sK)
            % dirichlet: contains restricted dofs DOFr x 1
            % dirichlet_values: contains displacement values DOFr x 1
            % C is a squared matrix
            dirichletDOFs = bc.dirichlet;
            nDOFr = size(dirichletDOFs, 1);
            Ct = zeros(nDOFr, sK);
            for i = 1:nDOFr
                Ct(i,dirichletDOFs(i)) = 1;
            end        
        end
    end
end