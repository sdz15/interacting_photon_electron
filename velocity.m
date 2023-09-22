function val = velocity(pMM,pMP,pPM,pPP)
    j_00 = j00(pMM,pMP,pPM,pPP);
    if (norm(j_00)>=1e-6)
        j_10 = j10(pMM, pMP, pPM, pPP);
        j_01 = j01(pMM, pMP, pPM, pPP);
        val = [j_10; j_01]./j_00;
        return
    end
    val = [0; 0];
end