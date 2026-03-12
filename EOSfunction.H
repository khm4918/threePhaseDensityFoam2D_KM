



void EROPM        //Internal energy according to the density and pressure of the fluid
(                 //see IJNMF 2016  Eq8, 9
    scalar& interE,
    scalar rho,
    scalar pre,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
){

double v=0.0; //specific volume
double qt=0.0;
double a[N],b[N];

v=1.0/rho;

for(int i= 0; i < N; i++){
qt+=Y[i]*Param_Eos[Ener_Ref][i];
}

for(int i= 0; i < N; i++){
a[i]=Y[i]*Param_Eos[Cv][i]*(pre+Param_Eos[Gamma_Pinf][i])/(pre+Param_Eos[Pinf][i]);
b[i] =Y[i]*Param_Eos[Cp_MCv][i]/(pre+Param_Eos[Pinf][i]);
}

double sumA=0.0,sumB=0.0;
for(int i= 0; i < N; i++){
sumA+=a[i];
sumB+=b[i];
}

/* Info<< "EROPM " << endl;
 Info<< "Param_Eos[Cv][0] " << Param_Eos[Cv][0]<< endl;
 Info<< "Param_Eos[Pinf][0] " << Param_Eos[Pinf][0]<< endl;
 Info<< "Param_Eos[Ener_Ref][0] " << Param_Eos[Ener_Ref][0]<< endl;
  Info<< "b[0] " << b[0]<< endl;

 Info<< "EROPM " << endl;
 Info<< "Param_Eos[Cv][1] " << Param_Eos[Cv][1]<< endl;
 Info<< "Param_Eos[Pinf][1] " << Param_Eos[Pinf][1]<< endl;
 Info<< "Param_Eos[Ener_Ref][1] " << Param_Eos[Ener_Ref][1]<< endl;
 Info<< "b[1] " << b[1]<< endl;

 Info<< "EROPM " << endl;
 Info<< "Param_Eos[Cv][2] " << Param_Eos[Cv][2]<< endl;
 Info<< "Param_Eos[Pinf][2] " << Param_Eos[Pinf][2]<< endl;
 Info<< "Param_Eos[Ener_Ref][2] " << Param_Eos[Ener_Ref][2]<< endl;
Info<< "b[2] " << b[2]<< endl;

   Info<< "sumA " << sumA << endl;
   Info<< "sumB " << sumB << endl;*/


//SG
//interE=qt+v*sumA/sumB;

//NASG
double sumCovol=0.0;

for(int i= 0; i < N; i++){
sumCovol+=Y[i]*Param_Eos[Covol][i];
}

interE=qt+(v-sumCovol)*sumA/sumB;  //see 2.8 CAF
//test

//for(int i= 0; i < N; i++){

//}

}



void ETOPM        //Internal energy according to the density and pressure of the fluid
(                 //see IJNMF 2016  Eq8, 9
    scalar& interE,
    scalar temperature,
    scalar pre,        
    double Y[],
    const label  N,
    double Param_Eos[][comp]
){

double v=0.0; //specific volume
double qt=0.0;
double a[N],b[N];

//v=1.0/rho;

for(int i= 0; i < N; i++){
qt+=Y[i]*Param_Eos[Ener_Ref][i];
}

for(int i= 0; i < N; i++){
a[i]=Y[i]*Param_Eos[Cv][i]*(pre+Param_Eos[Gamma_Pinf][i])/(pre+Param_Eos[Pinf][i]);
//b[i] =Y[i]*Param_Eos[Cp_MCv][i]/(pre+Param_Eos[Pinf][i]);
}

double sumA=0.0,sumB=0.0;
for(int i= 0; i < N; i++){
sumA+=a[i];
//sumB+=b[i];
}




interE=qt+sumA*temperature;      

//test

//for(int i= 0; i < N; i++){

//}

}









void ROPTM        //Density as a function of fluid pressure and temperature
(                 //see IJNMF Eq 2 3
    scalar& density,
    scalar pre,
    scalar Temperature,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
){

  //ROPT_SG=(Pression+Param_Eos(Pinf,Fluide))/(Param_Eos(Cp_MCv,Fluide)*Temperature)
  // ROPT_SG_AD=(Pression+Param_Eos(Pinf,Fluide))/(Param_Eos(Cp_MCv,Fluide)*Temperature)

 //  ROPT_NASG=(Pression+Param_Eos(Pinf,Fluide))/(Param_Eos(Cp_MCv,Fluide)*& Temperature+(Pression//+Param_Eos(Pinf,Fluide))*Param_Eos(Covol,Fluide))
  double parD[N];

  for(int i= 0; i < N; i++){ 
  //SG
 // parD[i]=(pre+Param_Eos[Pinf][i])/(Param_Eos[Cp_MCv][i]*Temperature);  //see 2.7 CAF

//NASG
//ROPT_NASG=(Pression+Param_Eos(Pinf,Fluide))/(Param_Eos(Cp_MCv,Fluide)*Temperature+(Pression+Param_Eos(Pinf,Fluide))*Param_Eos(Covol,Fluide))

   parD[i]=(pre+Param_Eos[Pinf][i])/(Param_Eos[Cp_MCv][i]*Temperature+(pre+Param_Eos[Pinf][i])*Param_Eos[Covol][i]);

  }

  double somme=0.0;

    for(int i= 0; i < N; i++){ 
  somme=somme+Y[i]/parD[i];
  } 

  density=1.0/somme;
}



void PROEM      // pressure according to Internal energy and density // See CF Eq 2.10, 2.11
(                
    scalar& pressure,
    scalar density,
    scalar interE,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
){

double qt=0.0, cvt=0.0, cpt=0.0, a=0.0, b=0.0, c=0.0, delta=0.0;

  for(int i= 0; i < N; i++){ 
  qt+=Y[i]*Param_Eos[Ener_Ref][i];
  cvt+=Y[i]*Param_Eos[Cv][i];
  cpt+=Y[i]*Param_Eos[Cp][i];
  }

double sumCovol=0.0;

for(int i= 0; i < N; i++){
sumCovol+=Y[i]*Param_Eos[Covol][i];
}

double v=0.0;

v=1.0/density;

double X[N];  

  for(int i= 0; i < N; i++){ 

  //SG
  //X[i] = Y[i]*Param_Eos[Cp_MCv][i]/cvt*(density*(interE-qt)-Param_Eos[Pinf][i]);


  //NASG
  X[i] = Y[i]*Param_Eos[Cp_MCv][i]/cvt*((1.0/(v-sumCovol))*(interE-qt)-Param_Eos[Pinf][i]);

  }  


double sumX=0.0;

double sumXG=0.0;

  for(int i= 0; i < N; i++){ 
  sumX+=X[i];
  }

  for(int i= 1; i < N; i++){ 
  sumXG+=X[i];
  }

  a=1.0;

  b=-(sumX-Param_Eos[Pinf][0]);

  c=-Param_Eos[Pinf][0]*sumXG;

  delta=b*b-4.0*a*c;

  pressure = (-b+Foam::sqrt(delta))/(2.0*a);

}





void TROPM        //temperature as a function of density and pressure
(                 //see CF Eq 2.8
    scalar& temperature,
    scalar density,
    scalar pressure,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
){
double v;
double a[N];
double sumA;

v=1.0/density;


  for(int i= 0; i < N; i++){ 

  a[i] = Y[i]*Param_Eos[Cp_MCv][i]/(pressure+Param_Eos[Pinf][i]);

  } 

  sumA=0.0;

  for(int i= 0; i < N; i++){ 

  sumA+=a[i];

  }

double sumCovol=0.0;

for(int i= 0; i < N; i++){
sumCovol+=Y[i]*Param_Eos[Covol][i];
}


//SG
//  temperature=v/sumA;


//NASG

temperature=(v-sumCovol)/sumA;

}

void APROM       
(                 
    double alphaI[],
    scalar pressure,
    scalar density,
    scalar Temperature,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
){

double rhok[N];
double cs[N];
double alpha[N];

//double alphaI[N];

//ROPT_SG=(Pression+Param_Eos(Pinf,Fluide))/(Param_Eos(Cp_MCv,Fluide)*Temperature)
for(int i= 0; i < N; i++){ 

//  ROPT_NASG=(Pression+Param_Eos(Pinf,Fluide))/(Param_Eos(Cp_MCv,Fluide)*&
//             Temperature+(Pression+Param_Eos(Pinf,Fluide))*Param_Eos(Covol,Fluide))

//SG
//rhok[i]=(pressure+Param_Eos[Pinf][i])/(Param_Eos[Cp_MCv][i]*Temperature);
//NASG
rhok[i]=(pressure+Param_Eos[Pinf][i])/(Param_Eos[Cp_MCv][i]*Temperature+(pressure+Param_Eos[Pinf][i])*Param_Eos[Covol][i]);

}

for(int i= 0; i < N; i++){ 
alpha[i]=density*Y[i]/rhok[i] ;
}

double sumA=0.0;

for(int i= 0; i < N; i++){ 
sumA+=alpha[i];
}

for(int i= 0; i < N; i++){ 
alphaI[i]=alpha[i]/sumA;
}

//Info << " Param_Eos[Gamma][0]" << Param_Eos[Gamma][0] << endl;
//Info << " Param_Eos[Gamma][1]" << Param_Eos[Gamma][1] << endl;
//Info << " Y[1]" << Y[1] << endl;
}




void CPROM        //calculate sound speed, see IJNMF Eq 3
(                 
    scalar& soundSpeed,
    scalar pressure,
    scalar density,
    scalar Temperature,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
){

double rhok[N];
double cs[N];
double alpha[N];

double alphaI[N];

//ROPT_SG=(Pression+Param_Eos(Pinf,Fluide))/(Param_Eos(Cp_MCv,Fluide)*Temperature)
for(int i= 0; i < N; i++){ 

//  ROPT_NASG=(Pression+Param_Eos(Pinf,Fluide))/(Param_Eos(Cp_MCv,Fluide)*&
//             Temperature+(Pression+Param_Eos(Pinf,Fluide))*Param_Eos(Covol,Fluide))

//SG
//rhok[i]=(pressure+Param_Eos[Pinf][i])/(Param_Eos[Cp_MCv][i]*Temperature);
//NASG
rhok[i]=(pressure+Param_Eos[Pinf][i])/(Param_Eos[Cp_MCv][i]*Temperature+(pressure+Param_Eos[Pinf][i])*Param_Eos[Covol][i]);

}

//CPRO_SG=SQRT(Param_Eos(Gamma,Fluide)*(Pression+Param_Eos(Pinf,Fluide))/Densite)
for(int i= 0; i < N; i++){ 

//SG
cs[i]=Foam::sqrt(Param_Eos[Gamma][i]*(pressure+Param_Eos[Pinf][i])/rhok[i]);

//NASG   //see Eq 8 general formulation
//  CPRO_NASG=SQRT(Param_Eos(Gamma,Fluide)*(Pression+Param_Eos(Pinf,Fluide))/(Densite*(1.0_DP-Densite*Param_Eos(Covol,Fluide))))
cs[i]=Foam::sqrt(Param_Eos[Gamma][i]*(pressure+Param_Eos[Pinf][i])/(rhok[i]*(1.0-rhok[i]*Param_Eos[Covol][i])));
}

for(int i= 0; i < N; i++){ 
alpha[i]=density*Y[i]/rhok[i] ;
}

double sumA=0.0;

for(int i= 0; i < N; i++){ 
sumA+=alpha[i];
}

for(int i= 0; i < N; i++){ 
alphaI[i]=alpha[i]/sumA;
}

double sumCS=0.0;

for(int i= 0; i < N; i++){ 
sumCS+=alphaI[i]/(rhok[i]*cs[i]*cs[i]);
}

soundSpeed=1.0/(Foam::sqrt(density*sumCS));


//Info << " Param_Eos[Gamma][0]" << Param_Eos[Gamma][0] << endl;
//Info << " Param_Eos[Gamma][1]" << Param_Eos[Gamma][1] << endl;
//Info << " Y[1]" << Y[1] << endl;
}




void TSAT
(                 
    scalar& TS,
    const scalar pressure,        
    const label  N,
    double Param_Eos[][comp]
){
double a,b,c;

//a=4.6543;
//b=1435.264;
//c=-64.848;

//oxygen
a=3.9523;
b=340.024;
c=-4.144;

TS=b/(a-Foam::log10(pressure*1e-5))-c;


}



void PSAT
(                 
    scalar& PS,
    const scalar temperature,        
    const label  N,
    double Param_Eos[][comp]
){

double a,b,c;

//a=4.6543;
//b=1435.264;
//c=-64.848;

//oxygen
a=3.9523;
b=340.024;
c=-4.144;

double abc;

abc=a-b/(temperature+c);

PS=Foam::pow(10.0,abc)*1e5;



}


