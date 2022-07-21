classdef ParaviewPostprocessor < handle
    
    % Its output is a legacy .vtk file. It works for simple fields, but it
    % cannot represent data at Gaussian points, as it requires the newer
    % XML implementation.
    
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
            obj.writeCoords();
            obj.writeConnec();
            obj.writeData();
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
            pieceN.setAttribute('NumberOfPoints', '1'); % !!
            pieceN.setAttribute('NumberOfCells', '2'); % !!
            gridN.appendChild(pieceN);

            % Creating Points
            pointsN = docNode.createElement('Points');
            pieceN.appendChild(pointsN);

            % Creating Points DataArray
            pointsDAN = docNode.createElement('DataArray');
            pointsDAN.setAttribute('type', 'Float32');
            pointsDAN.setAttribute('Name', 'Points');
            pointsDAN.setAttribute('NumberOfComponents', '3'); % !!
            pointsDAN.setAttribute('format', 'ascii');
            pointsN.appendChild(pointsDAN);

            % Creating Cells
            cellsN = docNode.createElement('Cells');
            pieceN.appendChild(cellsN);

            % Creating Cells Connectivity DataArray
            connecDAN = docNode.createElement('DataArray');
            connecDAN.setAttribute('type', 'Int32');
            connecDAN.setAttribute('Name', 'connectivity');
            connecDAN.setAttribute('format', 'ascii');
            cellsN.appendChild(connecDAN);

            % Creating Cells Offsets DataArray
            offsetsDAN = docNode.createElement('DataArray');
            offsetsDAN.setAttribute('type', 'Int32');
            offsetsDAN.setAttribute('Name', 'offsets');
            offsetsDAN.setAttribute('format', 'ascii');
            cellsN.appendChild(offsetsDAN);

            % Creating Cells Types DataArray
            typesDAN = docNode.createElement('DataArray');
            typesDAN.setAttribute('type', 'UInt8');
            typesDAN.setAttribute('Name', 'types');
            typesDAN.setAttribute('format', 'ascii');
            cellsN.appendChild(typesDAN);


            xmlwrite(docNode)
        end
        
        function writeCoords(obj)
        end
        
        function writeConnec(obj)
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

