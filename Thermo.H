const label NumberOfThermo=17;

/*scalarListList Param_Eos(NumberOfThermo);

forAll(Param_Eos, i){
//for(label i=0; i<NumberOfThermo; i++){
   Param_Eos[i].setSize(N);

}*/

double Param_Eos[NumberOfThermo][N];

/*const label Cv         =0;                        
const label Cp         =1;
const label Cp_MCv     =2;
const label Gamma      =3;
const label Pinf       =4;
const label Ener_Ref   =5;
// Values ​​calculated by code
const label Gamma_Mu   =6;
const label Gamma_Pinf =7;
const label Gamma_Pu   =8;
const label Rap_Gam_Pu =9;
const label Rap_Gam_Mu =10;
const label Covol      =11;
const label Emissivite =12;
const label Prandtl    =13;
const label Cond_Therm =14;
const label Visc_Dyn   =15;
const label Wmol       =16;*/

    /*CALL GetValue(config,nom_section,"Gamma",Param_Eos(Gamma,espece),found)
    CALL GetValue(config,nom_section,"Pinf",Param_Eos(Pinf,espece),found)
    CALL GetValue(config,nom_section,"Eref",Param_Eos(Ener_Ref,espece),found)
    CALL GetValue(config,nom_section,"Cv",Param_Eos(Cv,espece),found)      
    CALL GetValue(config,nom_section,"Masse Molaire",Param_Eos(Wmol,espece),found)    
    CALL GetValue(config,nom_section,"Covolume",Param_Eos(Covol,espece),found)*/

IOdictionary speciesDict
    (
        IOobject
        (
            "thermoData",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    );
    wordList speciesNamesCyp
    (
        speciesDict.lookup("species")
    );

   // forAll(speciesNamesCyp, s)
   for(int s = 0; s < N; s++)
    {
        word currentSpecie = speciesNamesCyp[s];
        Info<< "Doing stuff for specie: " << currentSpecie << endl;

        const dictionary& subDict = speciesDict.subDict(currentSpecie);

       dimensionedScalar Gamma_("Gamma", dimensionSet(0,0,0,0,0,0,0), scalar(0.0));
        dimensionedScalar Pinf_("Pinf", dimensionSet(0,0,0,0,0,0,0), scalar(0.0));
        dimensionedScalar Ener_Ref_("Ener_Ref", dimensionSet(0,0,0,0,0,0,0), scalar(0.0));
        dimensionedScalar Cv_("Cv", dimensionSet(0,0,0,0,0,0,0), scalar(0.0));
        dimensionedScalar Wmol_("Wmol", dimensionSet(0,0,0,0,0,0,0), scalar(0.0));
        dimensionedScalar Covol_("Covol", dimensionSet(0,0,0,0,0,0,0), scalar(0.0));

        Gamma_ = readScalar(subDict.lookup("Gamma"));
        Pinf_ = readScalar(subDict.lookup("Pinf"));
        Ener_Ref_ = readScalar(subDict.lookup("Ener_Ref"));
        Cv_ = readScalar(subDict.lookup("Cv"));
        Wmol_ = readScalar(subDict.lookup("Wmol"));
        Covol_ = readScalar(subDict.lookup("Covol"));

        Param_Eos[Gamma][s]=Gamma_.value();
        Param_Eos[Pinf][s]=Pinf_.value();
        Param_Eos[Ener_Ref][s]=Ener_Ref_.value();
        Param_Eos[Cv][s]=Cv_.value();
        Param_Eos[Wmol][s]=Wmol_.value();
        Param_Eos[Covol][s]=Covol_.value();
        
       // const dimensionedScalar Gamma("Uinf", dimensionSet(0,0,0,0,0,0,0), readScalar(speciesDict.lookup("Uinf")));
        
       /* Param_Eos[Gamma][s]=subDict.lookupEntry("Gamma");
        Param_Eos[Pinf][s]=subDict.lookup("Pinf");
        Param_Eos[Ener_Ref][s]=subDict.lookup("Ener_Ref");
        Param_Eos[Cv][s]=subDict.lookup("Cv");
        Param_Eos[Wmol][s]=subDict.lookup("Wmol");
        Param_Eos[Covol][s]=subDict.lookup("Covol");    */   
        
        

  /*Param_Eos(Cp,espece)=Param_Eos(Gamma,espece)*Param_Eos(Cv,espece)
    Param_Eos(Cp_MCv,espece)=Param_Eos(Cp,espece)-Param_Eos(Cv,espece)
    Param_Eos(Gamma_Mu,espece)=Param_Eos(Gamma,espece)-1._DP
    Param_Eos(Gamma_Pu,espece)=Param_Eos(Gamma,espece)+1._DP
    Param_Eos(Gamma_Pinf,espece)=Param_Eos(Gamma,espece)*Param_Eos(Pinf,espece)
    Param_Eos(Rap_Gam_Pu,espece)=0.5_DP*Param_Eos(Gamma_Pu,espece)/Param_Eos(Gamma,espece)
    Param_Eos(Rap_Gam_Mu,espece)=1._DP-Param_Eos(Rap_Gam_Pu,espece)*/

        Param_Eos[Cp][s]=Param_Eos[Gamma][s]*Param_Eos[Cv][s];
        Param_Eos[Cp_MCv][s]=Param_Eos[Cp][s]-Param_Eos[Cv][s];
        Param_Eos[Gamma_Mu][s]=Param_Eos[Gamma][s]-1.0;
        Param_Eos[Gamma_Pu][s]=Param_Eos[Gamma][s]+1.0;
        Param_Eos[Gamma_Pinf][s]=Param_Eos[Gamma][s]*Param_Eos[Pinf][s];
        Param_Eos[Rap_Gam_Pu][s]=0.5*Param_Eos[Gamma_Pu][s]/Param_Eos[Gamma][s];
        Param_Eos[Rap_Gam_Mu][s]=1.0-Param_Eos[Rap_Gam_Pu][s];

        Info<< "Gamma = " << Param_Eos[Gamma][s]<< endl;
        Info<< "Pinf= " << Param_Eos[Pinf][s] << endl;
        Info<< "Ener_Ref= " << Param_Eos[Ener_Ref][s] << endl;
        Info<< "Cv = " << Param_Eos[Cv][s] << endl;
        Info<< "Wmol= " << Param_Eos[Wmol][s] << endl;
        Info<< "Covol= " << Param_Eos[Covol][s] << endl;
        Info<< "Cp= " << Param_Eos[Cp][s] << endl;
        Info<< "Gamma_Mu= " << Param_Eos[Gamma_Mu][s] << endl;       
    }
