//calculate the mixture viscosity and heat conduct coefficients based on mass fraction

volScalarField nu = nu0*Foam::pow(T/T0,0.6759);   //mu=nu_ref*(temp/300).^0.69;

/*volScalarField MCpG=nuG;


       forAll(rho,cellI){   

           double YcellI[N];
        for(label j=0; j<N; j++){ 
         YcellI[j]=Y[j][cellI];
        }            

       for(label j=1; j<N; j++){

        if(T[cellI]>Tref){
        MCpG[cellI]+=YcellI[j]*varCp[j][cellI];
        }else{
        MCpG[cellI]+= Param_Eos[Cp][j]*YcellI[j];
        }

        }
 

                }


volScalarField kG = MCpG*nuG/Pr;     //lamda=Cp.*mu/Pr;


volScalarField nu=nuG;
volScalarField k=kG;*/

volScalarField k = nu;

        forAll(rho,cellI){   //caculate pressure from density and internal energy

           double YcellI[N];
        for(label j=0; j<N; j++){ 
         YcellI[j]=Y[j][cellI];
        }            

       /* double rhok[N];

        for(label i=0; i<N; i++){ 
        rhok[i]=(p[cellI]+Param_Eos[Pinf][i])/(Param_Eos[Cp_MCv][i]*T[cellI]+(p[cellI]+Param_Eos[Pinf][i])*Param_Eos[Covol][i]);
        }

        double alpha[N],alphaI[N];         
  
        for(int i= 0; i < N; i++){ 
        alpha[i]=rho[cellI]*YcellI[i]/rhok[i] ;
        }

       double sumA=0.0;

       for(int i= 0; i < N; i++){ 
       sumA+=alpha[i];
       }

       for(int i= 0; i < N; i++){ 
       alphaI[i]=alpha[i]/sumA;
       }

       double alphaLiquid=0.0;
       double alphaGas=0.0;

       alphaLiquid=alphaI[0];

       alphaGas=1.0-alphaLiquid;

       double nuL=1.965e-04;
       double kL=0.1518;

       nu[cellI]=nuL*alphaLiquid+nuG[cellI]*alphaGas;
        k[cellI]=kL*alphaLiquid+kG[cellI]*alphaGas;*/
        
        nu[cellI]=nu0.value()*YcellI[0]+nu1.value()*YcellI[1]+nu2.value()*YcellI[2]; //assume constant viscosity coefficient
        k[cellI]=k0.value()*YcellI[0]+k1.value()*YcellI[1]+k2.value()*YcellI[2];     //assume constant heat conduct coefficient
                }

//Info<<"mixture property"<<endl;

volTensorField tauMC("tauMC", nu*dev2(Foam::T(fvc::grad(U))));
volVectorField stress = fvc::laplacian(nu, U)+fvc::div(tauMC);
volScalarField viscoE = fvc::div(mesh.Sf()&fvc::interpolate(tauMC&U))
                        +fvc::div(mesh.Sf()&fvc::interpolate(nu*fvc::grad(U)&U))
                        +fvc::laplacian(k,T);
                        
Info<<"viscosity "<<endl;
