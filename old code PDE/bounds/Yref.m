function Ycont = Yref(Z0cont,Z1cont,Z2cont,A0, A1,X0,X1,G_Hopf,const_G_Hopf,...
    Es0,Es1, const_Es0, const_Es1)

 
Ycont1 = YHopfCont(A0, (A0+A1)/2,X0,(X0+X1)/2,G_Hopf,const_G_Hopf,...
    Es0,(Es0+Es1)/2, const_Es0, (const_Es0+const_Es1)/2);
Ycont2 = YHopfCont((A0+A1)/2, A1,(X0+X1)/2,X1,G_Hopf,const_G_Hopf,....
    (Es0+Es1)/2,Es1, (const_Es0+const_Es1)/2, const_Es1);

if any(Ycont1>((Z0cont + Z1cont -1)^2/(4*Z2cont)))
    Ycont1 = Yref(Z0cont,Z1cont,Z2cont,A0, (A0+A1)/2,X0,(X0+X1)/2,...
        G_Hopf,const_G_Hopf,Es0,(Es0+Es1)/2, const_Es0, (const_Es0+const_Es1)/2);
end
if any(Ycont2>((Z0cont + Z1cont -1)^2/(4*Z2cont)))
    Ycont2 = Yref(Z0cont,Z1cont,Z2cont,(A0+A1)/2, A1,(X0+X1)/2,X1,...
        G_Hopf,const_G_Hopf,(Es0+Es1)/2,Es1, (const_Es0+const_Es1)/2, const_Es1);
end

Ycont = max(Ycont1,Ycont2);
fprintf('Ycont computed, %f\n',max(Ycont));
