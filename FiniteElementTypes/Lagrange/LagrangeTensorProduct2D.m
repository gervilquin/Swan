classdef LagrangeTensorProduct2D < handle
   
    properties (Access = public)
        polynomialOrder
        n_vertices
        vertices
        ndofs
        nodes
        lagrangePolynomials
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = LagrangeTensorProduct2D(polynomialOrder)
            obj.init(polynomialOrder);
        end
        
        function plotShapeFunctions(obj)
            set(groot,'defaulttextinterpreter','latex');
            
            s.coord = obj.vertices;
            s.connec = [1 2 3;2 3 4];
            m = Mesh(s);
            for i=1:3
                m = m.remesh(2);
            end

            for s = 1:obj.ndofs
                figure()
                m.plot();
                trisurf(m.connec,m.coord(:,1),m.coord(:,2),obj.shapeFunctions{s}(m.coord(:,1),m.coord(:,2)));
                xlim([0 1]); ylim([0 1]); zlim([-0.5 1]);
                xlabel('x'); ylabel('y'); zlabel('z');
                title("s:"+string(s-1))
                grid on
            end
            
        end
        
    end
    
    
    methods (Access = private)
       
        function init(obj,polynomialOrder)
            obj.polynomialOrder = polynomialOrder;
            obj.computeVertices();
            obj.computeNdof();
            obj.computeNodes();
            obj.computeShapeFunctions();
        end
        
        function computeVertices(obj)
            obj.n_vertices = 4;
            obj.vertices = [0,0;0,1;1,0;1,1];
        end
        
        function computeNdof(obj)
            k = obj.polynomialOrder;
            obj.ndofs = (k+1)^2;
        end
        
        function computeNodes(obj)
            k = obj.polynomialOrder;
            node = zeros(obj.ndofs,2);
            for i = 1:k+1
                for j = 1:k+1
                    s = obj.computeMonomialIndeces(i,j);
                    node(s,1)=(i-1)/k;
                    node(s,2)=(j-1)/k;
                end
            end
            obj.nodes = node;
        end
        
        function s = computeMonomialIndeces(obj,i,j)
            s = i + (j-1)*(obj.polynomialOrder+1);
        end
                
        function computeShapeFunctions(obj)
            syms x y
            basisMonomialFormSym = obj.computeBasisInMonomialForm();
            basisMonomialForm = matlabFunction(basisMonomialFormSym);
            shapeFuncs = cell(obj.ndofs,1);
            for s = 1:obj.ndofs
                a = obj.computeShapeFunctionCoefficients(basisMonomialForm,s);
                shapeFuncs{s} = matlabFunction(basisMonomialFormSym*a,'Vars',[x y]);
            end
            
            obj.shapeFunctions = shapeFuncs;
        end
        
        function basisMonomialFormSym = computeBasisInMonomialForm(obj)
            k = obj.polynomialOrder;
            syms x y
            for i = 1:(k+1)
                for j = 1:(k+1)
                    s = obj.computeMonomialIndeces(i,j);
                    basisMonomialFormSym(s) = x^(i-1)*y^(j-1);
                end
            end
        end
        
        function coefs = computeShapeFunctionCoefficients(obj,X,s)
            A = obj.applyLinearFormInMonomialForm(X);
            b = obj.computeLinearFormValues(s);
            coefs = A\b;
        end
        
        function B = computeLinearFormValues(obj,s)
            B = zeros(obj.ndofs,1);
            B(s) = 1;
        end
        
        function A = applyLinearFormInMonomialForm(obj,X)
            A = zeros(obj.ndofs);
            for s = 1:obj.ndofs
                        node(:) = obj.nodes(s,:);
                        A(s,:) = X(node(1),node(2));
            end
        end
        
    end
    
end