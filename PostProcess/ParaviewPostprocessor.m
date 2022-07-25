classdef ParaviewPostprocessor < handle
    
    % Its output is a .vtu file. It can represent both simple fields (and
    % supposedly data at Gaussian points).
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        coord
        connec
        filename
        outputFile
        d_u
    end
    
    methods (Access = public)
        
        function obj = ParaviewPostprocessor(cParams)
            obj.init(cParams);
            obj.openFile();
            obj.createPiece();
%             obj.writeCoords();
%             obj.writeConnec();
%             obj.writeData();
            obj.saveFile();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh     = cParams.mesh;
            obj.coord    = cParams.mesh.coord;
            obj.connec   = cParams.mesh.connec;
            obj.filename = cParams.filename;
            obj.d_u = cParams.d_u;
        end
        
        function openFile(obj)
%             fullfile = strcat(obj.filename, '.vtu');
%             obj.outputFile = fopen(fullfile,'w');
        end
        
        function createPiece(obj)
            docNode = com.mathworks.xml.XMLUtils.createDocument('VTKFile');
            % Creating VTKFile
            fileN = docNode.getDocumentElement;
            fileN.setAttribute('type', 'UnstructuredGrid');
            fileN.setAttribute('version', '0.1');
            fileN.setAttribute('byte_order', 'LittleEndian');

            % Creating UnstructuredGrid
            gridN = docNode.createElement('UnstructuredGrid');
            fileN.appendChild(gridN);

            % Creating Piece
            pieceN = docNode.createElement('Piece');
            pieceN.setAttribute('NumberOfPoints', string(size(obj.coord,1))); %
            pieceN.setAttribute('NumberOfCells', string(size(obj.connec,1)));
            gridN.appendChild(pieceN);

            % Creating Points
            pointsN = docNode.createElement('Points');
            pieceN.appendChild(pointsN);

            % Creating Points DataArray
            coordStr = obj.writeCoords();
            pointsDAN = docNode.createElement('DataArray');
            pointsDAN.setAttribute('type', 'Float32');
            pointsDAN.setAttribute('Name', 'Points');
            pointsDAN.setAttribute('NumberOfComponents', '3'); % !!
            pointsDAN.setAttribute('format', 'ascii');
            pointsDAT = docNode.createTextNode(coordStr);
            pointsDAN.appendChild(pointsDAT);
            pointsN.appendChild(pointsDAN);

            % Creating Cells
            [connecStr, offsetStr, typesStr] = obj.writeConnec();
            cellsN = docNode.createElement('Cells');
            pieceN.appendChild(cellsN);

            % Creating Cells Connectivity DataArray
            connecDAN = docNode.createElement('DataArray');
            connecDAN.setAttribute('type', 'Int32');
            connecDAN.setAttribute('Name', 'connectivity');
            connecDAN.setAttribute('format', 'ascii');
            connecDAT = docNode.createTextNode(connecStr);
            connecDAN.appendChild(connecDAT);
            cellsN.appendChild(connecDAN);

            % Creating Cells Offsets DataArray
            offsetsDAN = docNode.createElement('DataArray');
            offsetsDAN.setAttribute('type', 'Int32');
            offsetsDAN.setAttribute('Name', 'offsets');
            offsetsDAN.setAttribute('format', 'ascii');
            offsetDAT = docNode.createTextNode(offsetStr);
            offsetsDAN.appendChild(offsetDAT);
            cellsN.appendChild(offsetsDAN);

            % Creating Cells Types DataArray
            typesDAN = docNode.createElement('DataArray');
            typesDAN.setAttribute('type', 'UInt8');
            typesDAN.setAttribute('Name', 'types');
            typesDAN.setAttribute('format', 'ascii');
            typesDAT = docNode.createTextNode(typesStr);
            typesDAN.appendChild(typesDAT);
            cellsN.appendChild(typesDAN);

            % Create PointData
            pointDataN = docNode.createElement('PointData');
            pieceN.appendChild(pointDataN);
            
            % Create Displacement DataArray
            dispStr = sprintf('%.4f %.4f %.4f\n', obj.d_u');
            displacementDAN = docNode.createElement('DataArray');
            displacementDAN.setAttribute('type', 'Float64');
            displacementDAN.setAttribute('Name', 'Displacement');
            displacementDAN.setAttribute('NumberOfComponents', '3'); % !!
            displacementDAN.setAttribute('format', 'ascii');
            displacementDAT = docNode.createTextNode(dispStr);
            displacementDAN.appendChild(displacementDAT);
            pointDataN.appendChild(displacementDAN);
            
            xmlwrite(docNode)
        end
        
        function node = createPointsNode(obj, docNode)
        end

        function coordStr = writeCoords(obj)
            nnodes = size(obj.coord,1);
            if (size(obj.coord,2) == 2)
                obj.coord = [obj.coord, zeros(nnodes,1)];
            end
            coordStr = sprintf('%.4f %.4f %.4f\n', obj.coord');
        end
        
        function [connecStr, offsetStr, typesStr] = writeConnec(obj)
            nelems = size(obj.connec,1);
            nnodEl = size(obj.connec,2);
            connecP = obj.connec-1;
            format = [repmat(' %d',1,nnodEl),'\n'];
            connecStr = sprintf(format, connecP');
            
            offset = nnodEl:nnodEl:nelems*nnodEl;
            offsetStr = sprintf('%d ', offset);
            
            vtkType = obj.getCellType();
            types = vtkType*ones(nelems,1);
            typesStr = sprintf('%d ', types);
        end
        
        function writeData(obj)
        end
        
        function saveFile(obj)
        end
        
        function t = getCellType(obj)
            switch obj.mesh.type
                case 'TRIANGLE'
                    t = 5;
                case 'QUAD'
                    t = 9;
                case 'TETRAHEDRA'
                    t = 10;
                case 'HEXAHEDRA'
                    t = 12;
            end
            
        end
        
    end
    
end

