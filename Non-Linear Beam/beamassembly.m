function [KG,KGT,FG] = beamassembly(KE,KTE,FE,NNE,NDOF,DOF,KG,KGT,FG)


NDOFE = NNE*NDOF; % Number of degrees of freedom per element
    for j = 1:NDOFE
     FG(DOF(j), 1) = FG(DOF(j),1) + FE(j,1);   
        for k = 1:NDOFE
       KG(DOF(j), DOF(k)) = KG(DOF(j),DOF(k)) + KE(j,k);
       KGT(DOF(j), DOF(k)) = KGT(DOF(j),DOF(k)) + KTE(j,k);
        end
    end

    