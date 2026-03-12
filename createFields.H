    Info<< "Reading thermophysical properties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

  /*  dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );*/

    dimensionedScalar gamma
    (
        transportProperties.lookup("gamma")
    );

    dimensionedScalar T0
    (
        transportProperties.lookup("T0")
    );


    dimensionedScalar nu0
    (
        transportProperties.lookup("nu0")
    );
    
    dimensionedScalar nu1
    (
        transportProperties.lookup("nu1")
    );

    dimensionedScalar nu2
    (
        transportProperties.lookup("nu2")
    );
    
    dimensionedScalar k0
    (
        transportProperties.lookup("nu0")
    );
    
    dimensionedScalar k1
    (
        transportProperties.lookup("nu1")
    );

    dimensionedScalar k2
    (
        transportProperties.lookup("nu2")
    );    
    
 /*   dimensionedScalar R
    (
        transportProperties.lookup("R")
    );

    dimensionedScalar Cv
    (
        transportProperties.lookup("Cv")
    );*/

    Info<< "Reading field pressure\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T\n" << endl;
  volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field rhoU\n" << endl;
    volVectorField M
    (
        IOobject
        (
            "M",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        rho*U,
        U.boundaryField().types()
    );

    Info<< "Calculating field total energy\n" << endl;
    volScalarField E
    (
        IOobject
        (
            "E",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
       p/(gamma.value()-1.0)+(rho*U&U)/2.0,
       p.boundaryField().types()
    );

    volScalarField interE
    (
        IOobject
        (
            "interE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
       (E-(rho*U&U)/2.0)/rho,
       p.boundaryField().types()
    );

    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        gamma*p/rho,
        zeroGradientFvPatchField<scalar>::typeName
    );

    volScalarField Mach
    (
        IOobject
        (
            "Mach",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mag(U)/sqrt(C),
        zeroGradientFvPatchField<scalar>::typeName
    );

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless, -1.0)
);







    Info<< "Reading/calculating face flux field phi\n" << endl;
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        linearInterpolate(rho*U) & mesh.Sf()
    );

    surfaceScalarField rhoFlux
    (
        IOobject
        (
            "rhoFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        // fluxes in mass conservation equation: \varrho \vec{u}
        // only initialisation!
        (linearInterpolate(rho*U) & mesh.Sf())
    );

    surfaceVectorField rhoUFlux
    (
        IOobject
        (
            "rhoUFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        // fluxes in momentum equation: \varrho \vec{u} \vec{u}
        // only initialisation!
        rhoFlux*linearInterpolate(U)
    );

    surfaceScalarField rhoEFlux
    (
        IOobject
        (
            "rhoEFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        // fluxes in total energy equation: \varrho E \vec{u}
        // only initialisation!
        (linearInterpolate(E*U) & mesh.Sf())
    );

    ///////// initial components
const    label N = 3;    
//const    label comp = 3;     
 
PtrList<volScalarField> Y(N);

for(int count = 0; count < N; count++)
{
    word nameTi ("Y" + name(count));
    Y.set
    (
        count,
        new volScalarField 
    (
            IOobject
            (
                nameTi,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

PtrList<volScalarField> rhoY(N);

for(int count = 0; count < N; count++)
{
    word nameTi ("rhoY" + name(count));
    rhoY.set
    (
        count,
        new volScalarField 
    (
            IOobject
            (
                nameTi,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            rho*Y[count]
        )
    );
}

PtrList<volScalarField> rhoYI(N);

for(int count = 0; count < N; count++)
{
    word nameTi ("rhoYI" + name(count));
    rhoYI.set
    (
        count,
        new volScalarField 
    (
            IOobject
            (
                nameTi,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rhoY[count]
        )
    );
}


PtrList<volScalarField> alpha(N);

for(int count = 0; count < N; count++)
{
    word nameTi ("alpha" + name(count));
    alpha.set
    (
        count,
        new volScalarField 
    (
            IOobject
            (
                nameTi,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Y[count]
        )
    );
}

PtrList<surfaceScalarField> rhoYFlux(N);

for(int count = 0; count < N; count++)
{
    word nameTi ("rhoYFlux" + name(count));
    rhoYFlux.set
    (
        count,
        new surfaceScalarField
    (
            IOobject
            (
                nameTi,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (linearInterpolate(rhoY[count]*U) & mesh.Sf())
        )
    );
}

PtrList<volScalarField> flv_rhoY(N);

for(int count = 0; count < N; count++)
{
    word nameTi ("flv_rhoY" + name(count));
    flv_rhoY.set
    (
        count,
        new volScalarField 
    (
            IOobject
            (
                nameTi,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::surfaceIntegrate(rhoYFlux[count])
        )
    );
}

  #include "Thermo.H"
Info<< "Finish  reading thermo "<< endl;  
 // #include "EOSfunction.H"

    ////////////end initial multicomponent field

scalar r0      = 2.3e-03;
scalar a0      = 0.1e-03;
int n          = 12;

scalar Ps      = 550e6;			//shock peak
scalar t_star  = 0.9e-6;		//shock decay time
scalar u_shock = 1900;

scalar cs      = 1500;

scalar r_arr   = 15.5e-03;
scalar r_shock = r_arr - u_shock*t_star;

scalar x_off = 0;
scalar y_off = 0;

const scalar SMALL    = 1e-12;
const scalar U_CAP    = 1.0e5;
const scalar pAmbient = 1.01e5;
const scalar TAmbient = 300.0;						
const scalar rhoAmbient = 1000;						//scalar rhoAmbient = 1030; according to J-boy
const scalar rho_air = 1.293;

scalar gamma_w    = Param_Eos[Gamma][0];
scalar pInf_w     = Param_Eos[Pinf][0];
scalar Ener_Ref_w = Param_Eos[Ener_Ref][0];
scalar Cv_w       = Param_Eos[Cv][0];
scalar Wmol_w     = Param_Eos[Wmol][0];
scalar Covol_w 	  = Param_Eos[Covol][0];

// ------------------------------------------------------------------------- //
// Correct stiffened-gas post-shock density

scalar c0 = Foam::sqrt(gamma_w * ((pAmbient + pInf_w)/(rhoAmbient*(1-(rhoAmbient*Covol_w)))));
scalar Ms = u_shock/c0;

Info<< "  gamma   = " << gamma_w << nl
    << "  pInf_w    = " << pInf_w  << nl
    << "  speed of sound    = " << c0  << nl
    << "  Water Covol    = " << Covol_w  << nl	
    << "  Cv      = " << Cv_w << nl << endl;

// ------------------------------------------------------------------------- //
// Initialize fields with debug protections

label count_rprime  = 0;
label count_rmag    = 0;
label count_pfried  = 0;
label count_rhoadbt = 0;


   scalar time=runTime.startTime().value();
   scalar dTime=runTime.deltaTValue();
    if(time<dTime) 
    {
          label np=40;


 Info<< "Finish initialization field" << endl;

           forAll(U,cellI)
        {
            scalar sumY0=0.0;
            scalar sumY1=0.0;
            scalar sumY2=0.0;
            scalar sumP=0.0;
            //vector sumV=vector::zero;
            scalar sumNum=0.0;
            if(mixMeshType[cellI]==0)
            {
                label pt[4];
                pt[0]=mixVertLabels[cellI][0];
                pt[1]=mixVertLabels[cellI][1];
                pt[2]=mixVertLabels[cellI][2];
                pt[3]=mixVertLabels[cellI][3];

                scalar xM[4],yM[4];
                for(label m=0;m<4;m++)
                {
                    xM[m]=mesh.points()[pt[m]].x(); 
                    yM[m]=mesh.points()[pt[m]].y();
                }

                scalar N[4];
                for(int i=0; i<=np; i++)
                {
                    for(int j=0; j<=np; j++)
                    {
                        scalar wgI = -1.0+double(i)*2.0/double(np);
                        scalar wgJ = -1.0+double(j)*2.0/double(np);
                        N[0]=(1.0-wgI)*(1.0-wgJ)/4.0;
                        N[1]=(1.0+wgI)*(1.0-wgJ)/4.0;
                        N[2]=(1.0+wgI)*(1.0+wgJ)/4.0;
                        N[3]=(1.0-wgI)*(1.0+wgJ)/4.0;    
                        scalar xx = N[0]*xM[0]+N[1]*xM[1]+N[2]*xM[2]+N[3]*xM[3];
                        scalar yy = N[0]*yM[0]+N[1]*yM[1]+N[2]*yM[2]+N[3]*yM[3];
						
						// Interface perturbation
						scalar azim  = Foam::atan2(yy, xx);
						scalar pertR = r0 - a0*Foam::cos(n*azim);
						scalar r2    = xx*xx + yy*yy;
						
                        if(r2 <= pertR*pertR)		//change to match interface
                        {
                            sumY0+=1.0e-8;
                            sumY1+=1.0e-8;
                            sumY2+=0.99999998;
                            sumP+=1.0e5;
							rho[cellI] = rho_air;
                        }
                        else
                        {
                            sumY0+=0.99999998;
                            sumY1+=1.0e-8;
                            sumY2+=1.0e-8;
                            sumP+=1.0e5;
							rho[cellI] = rhoAmbient;
                        }
                        sumNum+=1.0;
                    }
                }
            }
            else
            {
                label pt[3];
                pt[0]=mixVertLabels[cellI][0];
                pt[1]=mixVertLabels[cellI][1];
                pt[2]=mixVertLabels[cellI][2];

                scalar xM[3],yM[3];
                for(label m=0;m<3;m++)
                {
                   xM[m]=mesh.points()[pt[m]].x(); 
                   yM[m]=mesh.points()[pt[m]].y();
                } 

                scalar N[3];
               /* for(int i=0;i<=np;i++)
                {
                    for(int j=0;j<=np;j++)
                    {
                        vector w=vector(0,0,0)+vector(i,j,0)*1.0/double(np);
                        N[0]=w[0];
                        N[1]=w[1];
                        N[2]=1.0-w[0]-w[1];
                        if(N[2]<0.0) continue;                        
                        scalar xx = N[0]*xM[1]+N[1]*xM[2]+N[2]*xM[0];
                        scalar yy = N[0]*yM[1]+N[1]*yM[2]+N[2]*yM[0];
                        if(xx<=x0+(yy/sqrt(3.0)))
                        {
                            sumP+=116.5;
                            sumR+=8.0;
                        }
                        else
                        {
                            sumP+=1.0;
                            sumR+=1.4;
                        }
                        sumNum+=1.0;
                    }
                }*/
            }

             Y[0][cellI] = sumY0/sumNum;
             Y[1][cellI] = sumY1/sumNum;
             Y[2][cellI] = sumY2/sumNum;
            p[cellI]=sumP/sumNum;
            T[cellI]=TAmbient; 
            U[cellI]=vector(0.0,0.0,0.0); 
            //rho[cellI] = sumR/sumNum;
        }
        
       // p.correctBoundaryConditions();
       // rho.correctBoundaryConditions();                  
 
 scalar gabba = 0;                   
 scalar rhoPost = 0;                  
    
    forAll(p,cellI){
    scalar xx=mesh.C()[cellI].x();
    scalar yy=mesh.C()[cellI].y();
    
	// Position relative to shock center
    scalar x_shift = xx - x_off;
    scalar y_shift = yy - y_off;
    scalar r_mag   = Foam::sqrt(x_shift*x_shift + y_shift*y_shift);

    // --- Converging shock: r_prime positive inside shock front ---
    scalar r_prime = r_mag - r_shock;      // positive for r < r_shock
    scalar tEff    = r_prime / Foam::max(u_shock, SMALL);

    // Friedlander pressure
    scalar p_friedlander = pAmbient;
    if (tEff >= 0.0)
    {
        scalar arg       = tEff / Foam::max(t_star, SMALL);
        scalar friedTerm = (1.0 - arg) * Foam::exp(-arg);
        p_friedlander    = Ps * friedTerm;
    }

    // Debug and clamp
    if (r_prime < 1e-8 && count_rprime < 10)
    {
        Info << "Warning: r_prime very small at cell " << cellI << ": " << r_prime << endl;
        count_rprime++;
    }
    if (r_mag < 1e-8 && count_rmag < 10)
    {
        Info << "Warning: r_mag very small at cell " << cellI << ": " << r_mag << endl;
        count_rmag++;
    }
    if (p_friedlander < pAmbient && count_pfried < 10)
    {
        Info << "Warning: p_friedlander below ambient at cell " << cellI << ": " << p_friedlander << endl;
        count_pfried++;
        p_friedlander = pAmbient; // clamp
    }

    if (p_friedlander > pAmbient + SMALL)
    {
        // Initial post-shock density
		if (gabba == 0)
		{
			rhoPost += rhoAmbient / (1 - (p_friedlander - pAmbient)/(rhoAmbient * Foam::pow(u_shock, 2)));
			Info << " Setting max rho Post: " << rhoPost << endl;
			gabba += 1;
		}
		
        // Other post-shock density
        scalar rho_adbt = rhoAmbient / (1 - (p_friedlander - pAmbient)/(rhoAmbient * Foam::pow(u_shock, 2)));
		//Info << "  rho_shock   = " << rho_adbt << endl;
		//scalar rho_adbt = rhoAmbient / (1 - (Ps - pAmbient)/(rhoAmbient * Foam::pow(u_shock, 2)));

        if (rho_adbt < SMALL && count_rhoadbt < 10)
        {
            Info << "Warning: rho_adbt very small at cell " << cellI << ": " << rho_adbt << endl;
            count_rhoadbt++;
            rho_adbt = SMALL; // clamp
        }

        // Safe radii
        scalar safe_rprime = Foam::max(Foam::mag(r_prime), 1e-6);
        scalar safe_rmag   = Foam::max(Foam::mag(r_mag), 1e-6);

        // Velocity magnitude
		scalar U_mag = Foam::sqrt( (Ps - pAmbient) * (1.0/rhoAmbient - 1.0/rhoPost) );

        // Smooth near focus
        //if (r_mag < 0.5*r_shock)
        //    U_mag *= (r_mag / (0.5*r_shock));

        // Cap
        if (U_mag > U_CAP)
            U_mag = U_CAP;

        // Internal energy and temperature
		//Info << "Estimated Rho: " << rho_adbt << endl;
        scalar T_adbt = (1-(rho_adbt*Covol_w)) * (p_friedlander + pInf_w) / ((Cv_w * rho_adbt) * (gamma_w-1));
		scalar T_ambo = (1-(rhoAmbient*Covol_w)) * (pAmbient + pInf_w) / ((Cv_w * rhoAmbient) * (gamma_w-1));
		//Info << " T_ambo        = " << T_ambo << endl;
		//Info << "  T_shock    = " << T_adbt << endl;

        //scalar T_adbt = (1-(rho_adbt*Covol_w)) * (Ps + pInf_w) / ((Cv_w * rho_adbt) * (gamma_w-1));

        // Assign fields
        p[cellI] = p_friedlander;
		rho[cellI] = rho_adbt;
		T[cellI] = TAmbient + T_adbt - T_ambo;		

        scalar ux = -x_shift / safe_rmag;
        scalar uy = -y_shift / safe_rmag;
		U[cellI] = vector(U_mag*ux, U_mag*uy, 0.0);
    
    }  
    
    }                 
////////////////////////////////////////////Y


////////////////////////////////////////////Y







       forAll(rho,cellI){   //caculate rho from T and P for r < 10mm
	
        double YcellI[N];
        for(label j=0; j<N; j++){ 
         YcellI[j]=Y[j][cellI];
        }

   ROPTM       
    (                 
    rho[cellI],
    p[cellI],
    T[cellI],         
    YcellI,
    N,
    Param_Eos
   );

        }

  /*forAll(rho,cellI){   	  //calculate T from rho and p for r > 10mm
           
		   double YcellI[N];
        for(label j=0; j<N; j++){ 
         YcellI[j]=Y[j][cellI];
        }

       TROPM      
    (                 
    T[cellI],
    rho[cellI],
    p[cellI],         
    YcellI,
    N,
    Param_Eos
    );

        }*/




 Info<< "\nYa\n" << endl;
        forAll(Y,i){

             rhoY[i]=Y[i]*rho;

                    }
 Info<< "\nYb\n" << endl;
/////////////////////////////////////////////Y

        M = rho*U;
        M.correctBoundaryConditions();
        M.replace(vector::Z,0.0);

        forAll(rho,cellI){

        double YcellI[N];
        for(label j=0; j<N; j++){ 
         YcellI[j]=Y[j][cellI];
        }

    EROPM        
    (                 
    interE[cellI],
    rho[cellI],
    p[cellI],         
    YcellI,
    N,
    Param_Eos
   );



        }
        
        
                    forAll(rho,cellI){   //update alpha

           double YcellI[N];
        for(label j=0; j<N; j++){ 
         YcellI[j]=Y[j][cellI];
        }
        
        double AcellI[N];
        
        APROM      
     (                 
    AcellI,
    p[cellI],
    rho[cellI],
    T[cellI],         
    YcellI,
    N,
    Param_Eos
    );

    for(label j=0; j<N; j++){ 
        alpha[j][cellI]=AcellI[j];
      }

    }


      //  E = p/(gamma.value()-1.0)+(rho*U&U)/2.0;

       E=rho*(interE+0.5*(U&U));
       E.correctBoundaryConditions();



        U.write();
        p.write();
        rho.write();
        E.write();
        interE.write();
        T.write();
        Y[0].write();
		Y[1].write();
        Y[2].write();
    }
//exit(0);
