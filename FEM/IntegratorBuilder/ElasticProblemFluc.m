classdef ElasticProblemFluc < ElasticProblem

    properties (Access = public)
        sizePer
    end

    methods
        
        function Ch = computeChomog(obj)
            nstre = obj.material.nstre;
            basis = diag(ones(nstre,1));
            Ch = zeros(nstre,nstre);
            switch obj.solType
                case {'REDUCED'}
                    nelem = size(obj.material.C,3);
                    npnod = obj.displacementField.dim.nnodes;
                    ndofs = npnod*obj.displacementField.dim.ndimf;
                    ngaus = obj.quadrature.ngaus;
                    tStrn  = zeros(nstre,ngaus,nstre,nelem);
                    tStrss = zeros(nstre,ngaus,nstre,nelem);
                    tDisp  = zeros(nstre,ndofs);
                    for istre=1:nstre
                        obj.vstrain = basis(istre,:);
                        obj.solve();
                        vars = obj.computeStressStrainAndCh();
                        Ch(:,istre)         = vars.stress_homog;
                        tStrn(istre,:,:,:)  = vars.strain;
                        tStrss(istre,:,:,:) = vars.stress;
                        tDisp(istre,:)      = vars.d_u;
                        obj.assignVarsToPrint(istre);
                    end
                    obj.variables.Chomog  = Ch;
                    obj.variables.tstrain = tStrn;
                    obj.variables.tstress = tStrss;
                    obj.variables.tdisp   = tDisp;

                case {'MONOLITIC'}
                    basis = diag(ones(nstre,1));
                    for istre=1:nstre
                        obj.vstrain = basis(istre,:);
                        obj.solve();
                  
                    perDOFslave = obj.boundaryConditions.periodic_constrained;
                    obj.sizePer = size(perDOFslave, 1);
                    nEqperType = obj.sizePer/4;
                    L   = obj.variables.LangMult;
                    sigmaX     = 0;
                    sigmaY     = 0;
                    tauXY      = 0;
                    d1         = obj.sizePer;
                    LPer          = L(1:d1);
                    LDir       = L(d1+1:end);
                    for i = 1:nEqperType
                        sigmaX = sigmaX + LPer(i);
                    end
                    for i = nEqperType+1:2*nEqperType
                        tauXY = tauXY + LPer(i);
                    end
                    for i = 2*nEqperType+1:3*nEqperType
                        tauXY = tauXY + LPer(i);
                    end
                    for i = 3*nEqperType+1:4*nEqperType
                        sigmaY = sigmaY + LPer(i);
                    end
                    sigmaX      = sigmaX + LDir(1) + LDir(5);
                    sigmaY      = sigmaY + LDir(2) + LDir(4);
                    tauXY       = (tauXY + LDir(6) + LDir(3) - LDir(7) - LDir(8))/2;
                    Ch(istre, :) = [sigmaX; sigmaY; tauXY];
%                         Ch(:, istre) = obj.variables.React;
                    end
                    obj.variables.Chomog = Ch;
            end
        end

        function [L,  LDir, stressHomog] = computeStressHomog(obj, sol)
            nEqperType = obj.sizePer/4;
            sigmaX     = 0;
            sigmaY     = 0;
            tauXY      = 0;
            d1         = obj.sizeK+1;
            d2         = obj.sizeK + obj.sizePer;
            L          = sol(d1:d2);
            LDir       = sol(d2+1:end);
            for i = 1:nEqperType
                sigmaX = sigmaX + L(i);
            end
            for i = nEqperType+1:2*nEqperType
                tauXY = tauXY + L(i);
            end
            for i = 2*nEqperType+1:3*nEqperType
                tauXY = tauXY + L(i);
            end
            for i = 3*nEqperType+1:4*nEqperType
                sigmaY = sigmaY + L(i);
            end
            sigmaX      = sigmaX + LDir(1) + LDir(5);
            sigmaY      = sigmaY + LDir(2) + LDir(4);
            tauXY       = (tauXY + LDir(6) + LDir(3) - LDir(7) - LDir(8))/2;
            stressHomog = [sigmaX; sigmaY; tauXY];
        end

    end
end