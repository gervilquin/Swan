classdef StrainComputer < handle
    
    properties (Access = private)
        dim
        mesh
        geometry
        quadrature
        displacement
        interpolation
        dofsInElem
    end
    
    methods (Access = public)
        
        function obj = StrainComputer(cParams)
            obj.init(cParams)
            obj.createGeometry();
        end

        function strain = compute(obj)
            strain = obj.computeStrain();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.dim                = cParams.dim;
            obj.mesh               = cParams.mesh;
            obj.quadrature         = cParams.quadrature;
            obj.displacement       = cParams.displacement;
            obj.interpolation      = cParams.interpolation;
            obj.dofsInElem         = cParams.dofsInElem;
        end
       
        function createGeometry(obj)
            q   = obj.quadrature;
            int = obj.interpolation;
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end
        
        function strain = computeStrain(obj)
            nstre = obj.dim.nstre;
            nelem = obj.dim.nelem;
            ngaus = obj.dim.ngaus;
            nnode = obj.dim.nnode;
            nunkn = obj.dim.ndimField;
            idx = obj.dofsInElem;
            d_u = obj.displacement;
            strain = zeros(nstre,nelem,ngaus);
            for igaus = 1:ngaus
                Bmat = obj.computeB(igaus);
                for istre=1:nstre
                    for ievab=1:nnode*nunkn
%                             ievab = nunkn*(inode-1)+idime;
%                             disp(ievab)
                            B = squeeze(Bmat(istre,ievab,:));
                            u = d_u(idx(ievab,:));
                            strain(istre,:,igaus)=strain(istre,:,igaus)+(B.*u)';
                        
                    end
                end
            end
            strain = permute(strain, [3 1 2]);
        end

        function Bmat = computeB(obj,igaus)
            s.dim          = obj.dim;
            s.geometry     = obj.geometry;
            s.globalConnec = [];
            s.dofsInElem = [];
            Bcomp = BMatrixComputer(s);
            Bmat = Bcomp.computeBmat(igaus);
        end
    end
    
end