classdef AnalyticalFunction < L2Function
    
    properties (Access = public)
        ndimf
    end
    
    properties (Access = private)
        fHandle
        mesh
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = AnalyticalFunction(cParams)
            obj.init(cParams)
        end

        function fxV = evaluate(obj, xGLoc)
            xV = obj.mesh.computeXgauss(xGLoc);
            x = xV(1,:,:);
            y = xV(2,:,:);
            xVec = [x;y];
            fxV = obj.fHandle(xV);
        end

        function plot(obj)
            p.mesh   = obj.mesh;
            p.connec = obj.mesh.connec;
            proj = Projector_toP1(p);
            p1fun = proj.project(obj);
            p1fun.plot(obj.mesh)
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fHandle = cParams.fHandle;
            obj.ndimf   = cParams.ndimf;
            obj.mesh    = cParams.mesh;
        end
        
    end
    
end