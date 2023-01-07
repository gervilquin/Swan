classdef ObjectSaver < handle

    methods (Static)
        function s = saveObj(obj)
            s.fileName = obj.fileName;
            s.scale = obj.scale;
            s.dim = obj.dim;
            s.type = obj.type;
            s.nelem = obj.nelem;
            s.bc = obj.bc;
            s.mesh = obj.mesh;
            s.material = obj.material;
            s.ngaus = obj.ngaus;
            s.interpolationType = obj.interpolationType;
            s.loadedFile = obj.loadedFile;
        end

    end
end