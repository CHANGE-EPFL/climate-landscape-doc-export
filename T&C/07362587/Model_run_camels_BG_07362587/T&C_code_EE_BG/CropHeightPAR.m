function [hc_H,SAI_H,B_H,Nreserve_H,Preserve_H,Kreserve_H,...
    hc_L,SAI_L,B_L,Nreserve_L,Preserve_L,Kreserve_L] = ...
    ...
    CropHeightPAR...
    ...
    (Ccrown, hc_Htm1,SAI_Htm1,aSE_H,LAI_Htm1,LAIdead_H,B_Htm1,Nreserve_Htm1,Preserve_Htm1,Kreserve_Htm1,ManIH,Mpar_H,VegH_Param_Dyn,...
    hc_Ltm1,SAI_Ltm1,aSE_L,LAI_Ltm1,LAIdead_L,B_Ltm1,Nreserve_Ltm1,Preserve_Ltm1,Kreserve_Ltm1,ManIL,Mpar_L,VegL_Param_Dyn,OPT_SoilBiogeochemistry)

hc_H =  zeros(1, length(Ccrown));
hc_L =  zeros(1, length(Ccrown));
SAI_H =  zeros(1, length(Ccrown));
SAI_L =  zeros(1, length(Ccrown));

Nreserve_H =  zeros(1, length(Ccrown));
Preserve_H =  zeros(1, length(Ccrown));
Kreserve_H =  zeros(1, length(Ccrown));

Nreserve_L =  zeros(1, length(Ccrown));
Preserve_L =  zeros(1, length(Ccrown));
Kreserve_L =  zeros(1, length(Ccrown));
B_H = B_Htm1;
B_L = B_Ltm1;
for cl = 1:length(Ccrown)
    if aSE_H(cl) == 2
        [hc_H(1,cl)] = GrassHeight(LAI_Htm1(1, cl),LAIdead_H(1, cl));
    elseif aSE_H(cl) == 5
        [hc_H(1,cl),SAI_H(1,cl),B_H(1,:,:),~,Nreserve_H(1,cl),Preserve_H(1,cl),Kreserve_H(1,cl),~] =...
            CropHeightType(LAI_Htm1(1, cl),LAIdead_H(1, cl),cl,B_Htm1(1,:,:),...
            Ccrown,Nreserve_Htm1(1,cl),Preserve_Htm1(1,cl),Kreserve_Htm1(1,cl),...
            ManIH(1,:),Mpar_H,VegH_Param_Dyn,OPT_SoilBiogeochemistry);
    else
        hc_H(1,cl) = hc_Htm1(1,cl);
        SAI_H(1,cl) = SAI_Htm1(1,cl);
        B_H(1,cl,:) = B_Htm1(1,cl,:);
        Nreserve_H(1,cl) = Nreserve_Htm1(1,cl);
        Preserve_H(1,cl) = Preserve_Htm1(1,cl);
        Kreserve_H(1,cl) = Kreserve_Htm1(1,cl);
    end

    if aSE_L(cl) == 2
        [hc_L(1,cl)] = GrassHeight(LAI_Ltm1(1, cl),LAIdead_L(1, cl));
    elseif aSE_H(cl) == 5
        [hc_L(1,cl),SAI_L(1,cl),B_L(1,:,:),~,Nreserve_L(1,cl),Preserve_L(1,cl),Kreserve_L(1,cl),~] =...
            CropHeightType(LAI_Ltm1(1, cl),LAIdead_L(1, cl),cl,B_Ltm1(1,:,:),...
            Ccrown,Nreserve_Ltm1(1,cl),Preserve_Ltm1(1,cl),Kreserve_Ltm1(1,cl),...
            ManIL(1,:),Mpar_L,VegL_Param_Dyn,OPT_SoilBiogeochemistry);
    else
        hc_L(1,cl) = hc_Ltm1(1,cl);
        SAI_L(1,cl) = SAI_Ltm1(1,cl);
        B_L(1,cl,:) = B_Ltm1(1,cl,:);
        Nreserve_L(1,cl) = Nreserve_Ltm1(1,cl);
        Preserve_L(1,cl) = Preserve_Ltm1(1,cl);
        Kreserve_L(1,cl) = Kreserve_Ltm1(1,cl);
    end
end

