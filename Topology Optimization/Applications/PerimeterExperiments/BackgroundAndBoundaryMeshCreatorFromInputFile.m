classdef BackgroundAndBoundaryMeshCreatorFromInputFile < handle
    
    properties (Access = public)
        backgroundMesh
        boundaryMesh
    end
    
    properties (Access = private)
        inputFile
        isRectangularBox
    end
    
    methods (Access = public)
        
        function obj = BackgroundAndBoundaryMeshCreatorFromInputFile(cParams)
            obj.init(cParams);
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();            
        end
    
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.inputFile        = cParams.inputFile; 
           obj.isRectangularBox = cParams.isBackgroundMeshRectangularBox;
        end
        
        function createBackgroundMesh(obj)
            a.fileName         = obj.inputFile;
            m                 = FemDataContainer(a);
            s.coord            = m.mesh.coord;
            s.connec           = m.mesh.connec;
            obj.backgroundMesh = Mesh(s);                  
        end          
        
        function createBoundaryMesh(obj)
            eval(obj.inputFile);
            if exist('External_border_nodes','var') && ~isempty(External_border_nodes)
                s.borderNodes    = External_border_nodes;
                s.borderElements = External_border_elements;
                s.backgroundMesh = obj.backgroundMesh;
                s.type = 'FromData';
                b = BoundaryMeshCreator.create(s);
                obj.boundaryMesh = b.create();
            elseif obj.isRectangularBox
                s.backgroundMesh = obj.backgroundMesh;
                s.dimension = 1:obj.backgroundMesh.ndim;
                s.type = 'FromReactangularBox';
                bC = BoundaryMeshCreator.create(s);
                obj.boundaryMesh = bC.create();
            else
                
            end  
        end
            
        
    end
        
  
    
end