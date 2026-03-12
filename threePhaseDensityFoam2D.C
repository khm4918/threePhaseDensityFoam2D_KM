/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    sonicFoam

Description
    Transient solver for trans-sonic/supersonic, laminar or turbulent flow
    of a compressible gas.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "upwind.H"
#include "syncTools.H"
#include "vectorList.H"
#include "hexMatcher.H"
#include "prismMatcher.H"
#include "simpleControl.H"
#include "scalarMatrices.H"
#include "slicedSurfaceFields.H"
#include "volPointInterpolation.H"
#include "cyclicAMIPolyPatch.H"
#include "cyclicAMILduInterface.H"
#include <vector>

#include "directionInterpolate.H"

//#include "simpleControl.H"

//#include "psiChemistryCombustionModel.H"
#include "multivariateScheme.H"
//0.0 for others 1.0 for double mach
#define eps_M   2.5
//kpre  for precondition
#define kpre    1.0
//free stream mach number,   to specify when mach is small than 1.0
#define Mapre   3.0

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//globle
 //   #include "createTime.H"
//    #include "createMesh.H"
//const    label N = 3; 
//#include "Thermo.H"

const    label comp = 3;

void EROPM        //Internal energy according to the density and pressure of the fluid
(                 //see IJNMF 2016  Eq8, 9
    scalar& interE,
    scalar rho,
    scalar pre,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
); 

void ROPTM        //Density as a function of fluid pressure and temperature
(                 
    scalar& density,
    scalar pre,
    scalar Temperature,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
);

void PROEM      //Internal energy according to the density and pressure of the fluid
(              
    scalar& pressure,
    scalar density,
    scalar interE,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
);

void TROPM        //temperature as a function of density and pressure
(                 //see CF Eq 2.8
    scalar& temperature,
    scalar density,
    scalar pressure,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
);

void CPROM        //calculate sound speed, see IJNMF Eq 3
(                 
    scalar& soundSpeed,
    scalar pressure,
    scalar density,
    scalar Temperature,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
);

void ETOPM        //Internal energy according to the density and pressure of the fluid
(                 //see IJNMF 2016  Eq8, 9
    scalar& interE,
    scalar temperature,
    scalar pre,        
    double Y[],
    const label  N,
    double Param_Eos[][comp]
);

void APROM        //calculate VOF
(                 
    double alphaI[],
    scalar pressure,
    scalar density,
    scalar Temperature,         
    double Y[],
    const label  N,
    double Param_Eos[][comp]
);


void TSAT
(                 
    scalar& TS,
    const scalar pressure,        
    const label  N,
    double Param_Eos[][comp]
);

void PSAT
(                 
    scalar& PS,
    const scalar temperature,        
    const label  N,
    double Param_Eos[][comp]
);


void evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    double rhoYFluxFace[],
    
    const scalar pLeft,
    const scalar pRight,
    const vector ULeft,
    const vector URight,

    const scalar rhoLeft,
    const scalar rhoRight,

    const scalar TLeft,
    const scalar TRight,
    
    double YLeft[],
    double YRight[],

    const scalar pown,
    const scalar pnei,
    const vector Uown,
    const vector Unei,

    const scalar rhoown,
    const scalar rhonei,

    const scalar Town,
    const scalar Tnei,
    
    double Yown[],
    double Ynei[],

    const scalar kappaLeft,
    const scalar kappaRight,

    const vector Sf,
    const scalar magSf,
    const vector dotX,
    const scalar Konstant,
    const label  N,
    double Param_Eos[][comp]
);



void evaluateFluxOutlet
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    double rhoYFluxFace[],
    
    const scalar pLeft,
   scalar pRight,
    const vector ULeft,
   vector URight,

    const scalar rhoLeftG,
    scalar rhoRightG,

    const scalar TLeft,
   scalar TRight,
    
    double YLeft[],
    double YRight[],

    const scalar pown,
    const scalar pnei,
    const vector Uown,
    const vector Unei,

    const scalar rhoown,
    const scalar rhonei,

    const scalar Town,
    const scalar Tnei,
    
    double Yown[],
    double Ynei[],

    const scalar kappaLeft,
    const scalar kappaRight,

    const vector Sf,
    const scalar magSf,
    const vector dotX,
    const scalar Konstant,
    const label  N,
    double Param_Eos[][comp]
);

#include "thermoParaNumber.H"

#include "EOSfunction.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "storeInformation.H"

    #include "createFields.H"
    #include "initContinuityErrs.H"
    
    #include "createTimeControls.H"

   //  simpleControl simple(mesh);

    typedef List<vectorList> vectorListList;

     volVectorField gradR(fvc::grad(rho));

      volVectorField gradT(fvc::grad(T));
   // volVectorField gradRca(fvc::grad(rhoca));
    
    volVectorField gradE(fvc::grad(E));
    volVectorField gradP(fvc::grad(p));
    volTensorField gradM(fvc::grad(M));
    volTensorField gradU(fvc::grad(U));

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



   // label component=Y.size();
   // Info<< "test\n" <<component<<endl;


 PtrList<volVectorField> gradY(N);

for(int count = 0; count < N; count++)
{
    word nameTi ("gradY" + name(count));
    gradY.set
    (
        count,
        new volVectorField 
    (
            IOobject
            (
                nameTi,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (fvc::grad(Y[count]))
        )
    );
}



 typedef List<scalarField> scalarFieldList;

    scalarFieldList qlYLeft(N);
    scalarFieldList qlYRght(N);

   
for(int j=0; j<N; j++){

   qlYLeft[j].setSize(mesh.nFaces());
   qlYRght[j].setSize(mesh.nFaces());

}



   typedef List<vectorField> vectorFieldList;

   vectorFieldList rdY(N);

   //vectorField rdR(mesh.nPoints(),vector::zero);
for(int j=0; j<N; j++){

   rdY[j].setSize(mesh.nPoints());

}

    
    volScalarField rhod(rho);
volScalarField Td(T);
    volVectorField Ud(U);
    volScalarField Ed(E);
    volScalarField pd(p);


//volScalarField Yi(Y[0]);


PtrList<volScalarField> Yd(N);

for(int count = 0; count < N; count++)
{
    word nameTi ("Yd" + name(count));
    Yd.set
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
            Y[count]
        )
    );
}



    Info<< "\nStarting time loop\n" << endl;
    
    

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readTimeControls.H"       
        #include "CourantNo.H"
        #include "setDeltaT.H"

        volScalarField rho0(rho);      
        volVectorField M0(M);
        volScalarField E0(E);  

    
 for(label j=0; j<N; j++)
    {
        
    rhoYI[j]=rhoY[j];
	
    }

   
        for(int r=0; r<2; r++)    //rk loop
        {
            U = M/rho;
            U.correctBoundaryConditions();

            Ud = U;
            pd = p;
            rhod = rho;

             
              forAll(rho,cell){

             if(rho[cell]<=0.0){Info<< "\nnegtive\n"<<cell << endl;}

              }


          for(int j=0; j<N; j++){ 
          Yd[j]=Y[j];
            }

         
         volScalarField Y0(Y[0]);
         volScalarField Y1(Y[1]);
         volScalarField Y2(Y[2]);
  
        /* #include "Reconstruction.H"

         #include "ReconstructionT.H"
  
         #include "ReconstructionForY.H"*/
         
        // #include "ReconstructionForY.H"
         
	surfaceScalarField qlRLeft(interpolate(rho, pos,rho.name())); //rhoLeft
	surfaceScalarField qlRRght(interpolate(rho, neg,rho.name())); //rhoRght

	surfaceVectorField qlULeft(interpolate(U, pos, U.name()));
	surfaceVectorField qlURght(interpolate(U, neg, U.name()));

	surfaceScalarField qlPLeft(interpolate(p, pos, p.name()));
	surfaceScalarField qlPRght(interpolate(p, neg, p.name()));         
         
	surfaceScalarField qlTLeft(interpolate(T, pos, T.name()));
	surfaceScalarField qlTRght(interpolate(T, neg, T.name()));     
	
	
        surfaceScalarField qlYLeft0(interpolate(Y0, pos,Y0.name()));
        surfaceScalarField qlYRght0(interpolate(Y0, neg,Y0.name()));
	       
        surfaceScalarField qlYLeft1(interpolate(Y1, pos,Y1.name()));
        surfaceScalarField qlYRght1(interpolate(Y1, neg,Y1.name()));    
        
        surfaceScalarField qlYLeft2(interpolate(Y2, pos,Y2.name()));
        surfaceScalarField qlYRght2(interpolate(Y2, neg,Y2.name()));     

         #include "Riemann.H"


            volScalarField flv_R = fvc::surfaceIntegrate(rhoFlux);            
            volVectorField flv_M = fvc::surfaceIntegrate(rhoUFlux);
            volScalarField flv_E = fvc::surfaceIntegrate(rhoEFlux);


           // Info<< "\ne\n" << endl;
            for(label j=0; j<N; j++)
            {
               flv_rhoY[j]=fvc::surfaceIntegrate(rhoYFlux[j]);

             }
           //  Info<< "\nf\n" << endl;
            
           /*  forAll(p,i){
             Info << " flv_rhoY[0]" << flv_rhoY[0][i] << endl;
              }*/


        #include "viscosity.H"    
            
            
        // #include "surfaceTensionSplit.H"    
            

        dimensionedScalar dt=runTime.deltaT();      
        
                                 
           forAll(rho,cell){                             
            if(r==0){ rho[cell]=rho[cell]-dt.value()*flv_R[cell];}
            if(r==1){rho[cell]=0.5*(rho0[cell]+rho[cell]-dt.value()*flv_R[cell]);}
          
            }                                 
            

          forAll(rho,cell){ 
            if(r==0){ M[cell]=M[cell]-dt.value()*flv_M[cell]+dt.value()*stress[cell];}
            if(r==1){ M[cell]=0.5*(M0[cell]+M[cell]-dt.value()*flv_M[cell]+dt.value()*stress[cell]);}
           
            }

         forAll(rho,cell){ 
          if(r==0){ E[cell]=E[cell]-dt.value()*flv_E[cell]+dt.value()*viscoE[cell];}
          if(r==1){ E[cell]=0.5*(E0[cell]+E[cell]-dt.value()*flv_E[cell]+dt.value()*viscoE[cell]);}
       
         }

            for(label j=0; j<N; j++){ 
 	forAll(rho,cell){         

        if(r==0){ rhoY[j][cell]=rhoY[j][cell]-dt.value()*flv_rhoY[j][cell];}

	if(r==1){ rhoY[j][cell]=0.5*(rhoYI[j][cell]+rhoY[j][cell]-dt.value()*flv_rhoY[j][cell]); }


            }
                           }



            rho.correctBoundaryConditions();
   
            for(label j=0; j<N; j++){ 
         rhoY[j].correctBoundaryConditions();
            }
            


            M.correctBoundaryConditions();
            E.correctBoundaryConditions();
          // M.replace(vector::Y,0.0); // 1D problem
            M.replace(vector::Z,0.0);  // 2D problem
           
            
            
            U=M/rho;
            
            U.correctBoundaryConditions();


             forAll(Y,j){

             Y[j]=rhoY[j]/rho;                      //add to check negtive value
 

            }


            for(label j=0; j<N; j++){ 
         Y[j].correctBoundaryConditions();
            }

            
           interE=(E-(rho*U&U)/2.0)/rho;

           interE.correctBoundaryConditions();
            
         //    p = (gamma.value() - 1.0)*(E - 0.5*magSqr(M)/rho);

        
   forAll(rho,cellI){   //caculate pressure from density and internal energy

           double YcellI[N];
        for(label j=0; j<N; j++){ 
         YcellI[j]=Y[j][cellI];
        }

       PROEM      
    (                 
    p[cellI],
    rho[cellI],
    interE[cellI],         
    YcellI,
    N,
    Param_Eos
    );


        }


            p.correctBoundaryConditions();
            
            
           // T=p/R/rho;
            
       //     T.boundaryField()=p.boundaryField()/rho.boundaryField();

          forAll(rho,cellI){   //caculate pressure from density and internal energy

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


        }

  
  //#include "enforceInlet.H"
  
        rhoPhi = rhoFlux;
       T.correctBoundaryConditions();
       E.correctBoundaryConditions();

          forAll(rho,cellI){   //caculate pressure from density and internal energy

           double YcellI[N];
        for(label j=0; j<N; j++){ 
         YcellI[j]=Y[j][cellI];
        }

        CPROM      
(                 
    C[cellI],
    p[cellI],
    rho[cellI],
    T[cellI],         
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




 //exit(0);


        } //RK finish

      //  Info << " Min(totalH) = " << min(totalH).value() << " Max(totalH) = " << max(totalH).value() << endl;
             dimensionedScalar CSMALL("CSMALL",C.dimensions(),SMALL);
      //  volScalarField CM = gamma*p/rho;
       // C = max(CSMALL,CM);
        C.correctBoundaryConditions(); 
        Info << " Min(C) = " << min(C).value() << " Max(C) = " << max(C).value() << endl;
        Mach = mag(U)/C;
        Info << " Min(Mach) = " << min(Mach).value() << " Max(Mach) = " << max(Mach).value() << endl;
        volScalarField div = fvc::div(rhoPhi);
        Info << " Min(div) = " << min(div).value() << " Max(div) = " << max(div).value() << endl;

        Info << " Min(rho) = " << min(rho).value() << " Max(rho) = " << max(rho).value() << endl;
        
       // Info << " Min(ca) = " << min(ca).value() << " Max(ca) = " << max(ca).value() << endl;
        Info << " Min(T) = " << min(T).value() << " Max(T) = " << max(T).value() << endl;
        
        Info << " Min(U) = " << min(U).value() << " Max(U) = " << max(U).value() << endl;
        Info << " Min(E) = " << min(E).value() << " Max(E) = " << max(E).value() << endl;
        Info << " Min(p) = " << min(p).value() << " Max(p) = " << max(p).value() << endl;

        for(label j=0; j<N; j++){ 
      Info << " Min(Y) = " << min(Y[j]).value() << " Max(Y) = " << max(Y[j]).value() << endl;
       Info << " Min(rhoY) = " << min(rhoY[j]).value() << " Max(rhoY) = " << max(rhoY[j]).value() << endl;

// exit(0);




        }


  //   #include "phaseChange.H"



       
        runTime.write();
        



        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

//        runTime.writeNow();
  //   exit(0);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //




//---SLAU2
void evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    double rhoYFluxFace[],
    
    const scalar pLeft,
    const scalar pRight,
    const vector ULeft,
    const vector URight,

    const scalar rhoLeftG,
    const scalar rhoRightG,

    const scalar TLeft,
    const scalar TRight,
    
    double YLeft[],
    double YRight[],

    const scalar pown,
    const scalar pnei,
    const vector Uown,
    const vector Unei,

    const scalar rhoown,
    const scalar rhonei,

    const scalar Town,
    const scalar Tnei,
    
    double Yown[],
    double Ynei[],

    const scalar kappaLeft,
    const scalar kappaRight,

    const vector Sf,
    const scalar magSf,
    const vector dotX,
    const scalar Konstant,
    const label  N,
    double Param_Eos[][comp]
)
{
    // bounding variables
    const scalar rhoMin = SMALL;

    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

 
  
  const scalar Q=0.0;
    

double YTilde[3];

  for(label j=0; j<N; j++){ 
   YTilde[j]=0.5*(YLeft[j]+YRight[j]);

   }

scalar rhoLeft,rhoRight;

   ROPTM       
    (                 
    rhoLeft,
    pLeft,
    TLeft,         
    YLeft,
    N,
    Param_Eos
   );

   ROPTM       
    (                 
    rhoRight,
    pRight,
    TRight,         
    YRight,
    N,
    Param_Eos
   );

if(fabs(rhoLeftG-rhoRightG)<fabs(rhoLeft-rhoRight)){

rhoLeft=rhoLeftG;

rhoRight=rhoRightG;


}

/*scalar TLeft,TRight;

TROPM      
(         
    TLeft,
    rhoLeft,
    pLeft,         
    YLeft,
    N,
    Param_Eos
);

TROPM      
(         
    TRight,
    rhoRight,
    pRight,         
    YRight,
    N,
    Param_Eos
);*/


scalar rhoTilde=0.5*(rhoLeft+rhoRight);
      
scalar interELeft, interERight;

 ETOPM      
( 
    interELeft,
    TLeft,
    pLeft,         
    YLeft,
    N,
    Param_Eos
);

 ETOPM      
( 
    interERight,
    TRight,
    pRight,         
    YRight,
    N,
    Param_Eos
);


    const scalar HLeft = (rhoLeft*interELeft + (rhoLeft*ULeft&ULeft)/2.0+pLeft)/rhoLeft;
    const scalar HRight = (rhoRight*interERight+ (rhoRight*URight&URight)/2.0+pRight)/rhoRight; 

// const scalar rhoELeft  = rhoLeft*interELeft + (rhoLeft*ULeft&ULeft)/2.0;           
// const scalar rhoERight = rhoRight*interERight+ (rhoRight*URight&URight)/2.0; 


scalar aLeft,aRight;

CPROM      
(                 
    aLeft,
    pLeft,
    rhoLeft,
    TLeft,         
    YLeft,
    N,
    Param_Eos
);

CPROM      
(                 
    aRight,
    pRight,
    rhoRight,
    TRight,         
    YRight,
    N,
    Param_Eos
);

scalar aTilde = 0.5*(aLeft+aRight);


/*
CPROM      
(                 
    aTilde,
    0.5*(pRight+pLeft),
    0.5*(rhoRight+rhoLeft),
    0.5*(TRight+TLeft),         
    YTilde,
    N,
    Param_Eos
);*/

    const scalar UnLeft=ULeft.x()*normalVector.x()+ULeft.y()*normalVector.y()+ULeft.z()*normalVector.z();
    const scalar UnRight=URight.x()*normalVector.x()+URight.y()*normalVector.y()+URight.z()*normalVector.z();
    
    const scalar abUnLeft=fabs(UnLeft);
    const scalar abUnRight=fabs(UnRight);

/*
scalar aStarL, aStarR;

aStarL=aLeft*aLeft/max(aLeft,abUnLeft);

aStarR=aRight*aRight/max(aRight,abUnRight);

aTilde=min(aStarL,aStarR);*/




 double rhoYLeft[N];
 double rhoYRight[N];

    for(label j=0; j<N; j++){ 

    rhoYLeft[j]=rhoLeft*YLeft[j];
    rhoYRight[j]=rhoRight*YRight[j];

   /* rhoYLeft[j]=rhoown*Yown[j];
    rhoYRight[j]=rhonei*Ynei[j];*/


    }


    const scalar MaLeft=(ULeft.x()*normalVector.x()+ULeft.y()*normalVector.y()+ULeft.z()*normalVector.z())/aTilde;
    const scalar MaRight=(URight.x()*normalVector.x()+URight.y()*normalVector.y()+URight.z()*normalVector.z())/aTilde;
    
    const scalar MaParTran=Foam::sqrt(0.5*(ULeft.x()*ULeft.x()+ULeft.y()*ULeft.y()+ULeft.z()*ULeft.z()+URight.x()*URight.x()+URight.y()*URight.y()+URight.z()*URight.z()))/aTilde;
    const scalar MaPar=min(1.0, MaParTran);
    
    const scalar xbar=(1.0-MaPar)*(1.0-MaPar);
    
    //calculate mass flux
    const scalar gLeft=-max(min(MaLeft,0.0),-1.0);
    const scalar gRight=min(max(MaRight,0.0),1.0);
    const scalar g=gLeft*gRight;
    
    const scalar meanNormalU=(rhoLeft*abUnLeft+rhoRight*abUnRight)/(rhoLeft+rhoRight);
    
    const scalar abMeanULeft=(1.0-g)*meanNormalU+g*abUnLeft;
    const scalar abMeanURight=(1.0-g)*meanNormalU+g*abUnRight;
    
   // const scalar massFlux=0.5*(rhoLeft*(UnLeft+abMeanULeft)+rhoRight*(UnRight-abMeanURight)-xbar*(pRight-pLeft)/aTilde);


    //another way to calculate massflux


   const scalar massFlux=0.5*(rhoLeft*UnLeft+rhoRight*UnRight-meanNormalU*(rhoRight-rhoLeft))*(1.0-g)-0.5*xbar/aTilde*(pRight-pLeft);




    //calculate pressure flux
    const scalar Ma1PlusLeft   = 0.5*(1.0+sign(MaLeft));
    const scalar Ma1MinusRight = 0.5*(1.0-sign(MaRight));
   // Info<<"MaLeft"<<MaLeft<<"signMaLeft"<<sign(MaLeft)<<endl;
    
    
    const scalar Ma2PlusLeft = 0.25*(MaLeft+1.0)*(MaLeft+1.0)*(2.0-MaLeft);
    const scalar Ma2MinusRight=0.25*(MaRight-1.0)*(MaRight-1.0)*(2.0+MaRight);
    
    scalar fpLeft;
    scalar fpRight;
    
    if(fabs(MaLeft)>=1.0){
       fpLeft=Ma1PlusLeft;
    }
    else{
       fpLeft=Ma2PlusLeft;
    }
    
       if(fabs(MaRight)>=1.0){
       fpRight=Ma1MinusRight;
    }
    else{
       fpRight= Ma2MinusRight;
    } 


  // const scalar pFlux=0.5*(pLeft+pRight)+0.5*(fpLeft-fpRight)*(pLeft-pRight)+0.5*(1.0-xbar)*(fpLeft+fpRight-1.0)*(pLeft+pRight);
  
  //slau2:

//    const scalar pFlux=0.5*(pLeft+pRight)+0.5*(fpLeft-fpRight)*(pLeft-pRight)+0.5*MaParTran*(fpLeft+fpRight-1.0)*(pLeft+pRight);
  
//slau2  real gas
 //  const scalar pFlux=0.5*(pLeft+pRight)+0.5*(fpLeft-fpRight)*(pLeft-pRight)+0.5*MaParTran*(fpLeft+fpRight-1.0)*(aTilde*rhoTilde);

  const scalar pFlux=0.5*(pLeft+pRight)+0.5*(fpLeft-fpRight)*(pLeft-pRight)+1.0*MaParTran*aTilde*(fpLeft+fpRight-1.0)*(aTilde*rhoTilde);


//slau1 real gas
//const scalar pFlux=0.5*(pLeft+pRight)+0.5*(fpLeft-fpRight)*(pLeft-pRight)+0.5*(1.0-xbar)*(fpLeft+fpRight-1.0)*(aTilde*rhoTilde);

    rhoFlux= (0.5*(massFlux+fabs(massFlux))+0.5*(massFlux-fabs(massFlux)))*magSf;
    rhoUFlux=(0.5*(massFlux+fabs(massFlux))*ULeft+0.5*(massFlux-fabs(massFlux))*URight+pFlux*normalVector)*magSf;
    rhoEFlux=(0.5*(massFlux+fabs(massFlux))*HLeft+0.5*(massFlux-fabs(massFlux))*HRight)*magSf;

   for(label j=0; j<N; j++){ 

    rhoYFluxFace[j]=(0.5*(massFlux+fabs(massFlux))*YLeft[j]+0.5*(massFlux-fabs(massFlux))*YRight[j])*magSf;

    }

/*
    rhoFlux  = (convectionSpeed*rhoState)*magSf;
    
  //  rhocaFlux  = (convectionSpeed*rhocaState)*magSf;

    rhoUFlux = (convectionSpeed*rhoUState+pState*normalVector)*magSf;
    rhoEFlux = (convectionSpeed*(rhoEState+pState)+pState*(dotX & normalVector))*magSf;  


   for(label j=0; j<N; j++){ 

    rhoYFluxFace[j]=(convectionSpeed*rhoYState[j])*magSf;

    }*/

                 
}

void evaluateFluxOutlet
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    double rhoYFluxFace[],
    
    const scalar pLeft,
   scalar pRight,
    const vector ULeft,
   vector URight,

    const scalar rhoLeftG,
    scalar rhoRightG,

    const scalar TLeft,
   scalar TRight,
    
    double YLeft[],
    double YRight[],

    const scalar pown,
    const scalar pnei,
    const vector Uown,
    const vector Unei,

    const scalar rhoown,
    const scalar rhonei,

    const scalar Town,
    const scalar Tnei,
    
    double Yown[],
    double Ynei[],

    const scalar kappaLeft,
    const scalar kappaRight,

    const vector Sf,
    const scalar magSf,
    const vector dotX,
    const scalar Konstant,
    const label  N,
    double Param_Eos[][comp]
)
{
    // bounding variables
    const scalar rhoMin = SMALL;

    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

  double rhoYLeft[N];
 double rhoYRight[N];
  
  const scalar Q=0.0;
    
 // scalar scal=mag(ULeft);
scalar scal= ULeft & normalVector;


  scalar pOut=3.0e6;
  scalar TOut=150.0;
  scalar y0=1e-8; scalar y1=1e-8; scalar y2=0.99999998;

  double YOut[N];

   YOut[0]=y0; YOut[1]=y1; YOut[2]=y2; 

scalar rhoLeft,rhoRight;


   ROPTM       
    (                 
    rhoLeft,
    pLeft,
    TLeft,         
    YLeft,
    N,
    Param_Eos
   );

scalar aLeft;

CPROM      
(                 
    aLeft,
    pLeft,
    rhoLeft,
    TLeft,         
    YLeft,
    N,
    Param_Eos
);

scalar interELeft,interERight;

 ETOPM      
( 
    interELeft,
    TLeft,
    pLeft,         
    YLeft,
    N,
    Param_Eos
);

scalar rhoELeft  = rhoLeft*interELeft + (rhoLeft*ULeft&ULeft)/2.0; 





      /*  Wr(VP_Rho) = Wl(VP_Rho)+(pstar-Wl(VP_P))/((CPROM(Wl(VP_P),Wl(VP_Rho),Wl(VP_Y0+1:VP_Y0+Nbre_Yk)))**2)
     	Wr(VP_U0+1:VP_U0+N_dim)=Wl(VP_U0+1:VP_U0+N_dim)+(Wl(VP_P)-pstar)  &
            /(Wl(VP_Rho)*CPROM(Wl(VP_P),Wl(VP_Rho),Wl(VP_Y0+1:VP_Y0+Nbre_Yk)))*norm(:)
        scal2=DOT_PRODUCT(Wr(VP_U0+1:VP_U0+N_dim),norm)*/

rhoRight=rhoLeft+(pOut-pLeft)/(aLeft*aLeft);

URight=ULeft+(pLeft-pOut)/(rhoLeft*aLeft)*(normalVector);

//scalar scal2=mag(URight);
scalar scal2=URight & normalVector;

    scalar convectionSpeed = 0.0;
    scalar rhoState = 0.0;
    vector rhoUState = vector::zero;
    scalar rhoEState = 0.0;

     scalar pState=0.0;

double rhoYState[N];

 vector URelLeft  = ULeft  - dotX;
  vector URelRight = URight - dotX;

    // compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
  scalar qLeft  = (URelLeft  & normalVector);
 scalar qRight = (URelRight & normalVector);

 vector rhoULeft  = rhoLeft *ULeft;
   vector rhoURight = rhoRight*URight;


if(scal>aLeft){

/*URight=ULeft;
rhoRight=rhoLeft;

 ETOPM      
( 
    interERight,
    TOut,
    pOut,         
    YOut,
    N,
    Param_Eos
);*/
           /*Etr=Wr(VP_E)+0.5_DP*SUM(Wl(VP_U0+1:VP_U0+N_dim)**2)
           fxchapo(VC_roY0+1:VC_roY0+Nbre_Yk) = Wl(VP_Rho)*Wl(VP_Y0+1: VP_Y0+Nbre_Yk)*scal
           fxchapo(VC_roscal0+1) = Wl(VP_scal0+1)*scal
           fxchapo(VC_roscal0+2) = Wl(VP_Rho)*Wl(VP_scal0+2)*scal
 	   fxchapo(VC_Masse) = Wl(VP_Rho)*scal
   	   fxchapo(VC_QDM0+1: VC_QDM0+N_dim)  = (Wl(VP_Rho)*scal*Wl(VP_U0+1:VP_U0+N_dim)+Wl(VP_P)*norm(:))
           fxchapo(VC_E)     = (Wl(VP_Rho)*Etr+Wl(VP_P))*scal*/
    /*rhoFlux  = (convectionSpeed*rhoState)*magSf;
    
    rhoUFlux = (convectionSpeed*rhoUState+pState*normalVector)*magSf;
    rhoEFlux = (convectionSpeed*(rhoEState+pState)+pState*(dotX & normalVector))*magSf;  


   for(label j=0; j<N; j++){ 

    rhoYFluxFace[j]=(convectionSpeed*rhoYState[j])*magSf;

    }*/

convectionSpeed=scal;

rhoState=rhoLeft;

rhoUState=rhoULeft;

pState=pLeft;

rhoEState=rhoELeft;

    for(label j=0; j<N; j++){ 

    rhoYLeft[j]=rhoLeft*YLeft[j];
 //   rhoYRight[j]=rhoRight*YRight[j];

   /* rhoYLeft[j]=rhoown*Yown[j];
    rhoYRight[j]=rhonei*Ynei[j];*/

    }

    for(label j=0; j<N; j++){ 

    rhoYState[j]=rhoYLeft[j];

    }


}
else if(scal2<0.0){

  /*          Wr(VP_P) = pstar
           Wr(VP_Rho) = ROPTM(CLPOut%P,CLPOut%T,CLPOut%Yk)
           Wr(VP_Y0+1:VP_Y0+Nbre_Yk) = CLPOut%Yk
           Wr(VP_scal0+1:VP_scal0+Nbre_Yk) = CLPOut%scal
           Wr(VP_E) = EROPM(Wr(VP_Rho),Wr(VP_P),Wr(VP_Y0+1:VP_Y0+Nbre_Yk))
           Wr(VP_U0+1:VP_U0+N_dim)= Wl(VP_U0+1:VP_U0+N_dim)-2.0_DP*scal*norm(:)
!           Wr(:) = Wl(:)
!          scal=DOT_PRODUCT(Wl(VP_U0+1:VP_U0+N_dim),norm(:))
!           Wr(VP_U0+1:VP_U0+N_dim)=Wl(VP_U0+1:VP_U0+N_dim)-2.0_DP*scal*norm(:)

           Etr=Wr(VP_E)+0.5_DP*SUM(Wr(VP_U0+1:VP_U0+N_dim)**2)

  	fxchapo(VC_roY0+1:VC_roY0+Nbre_Yk) = Wr(VP_Rho)*Wr(VP_Y0+1:VP_Y0+Nbre_Yk)*scal2
  	fxchapo(VC_roscal0+1) = Wr(VP_scal0+1)*scal2
  	fxchapo(VC_roscal0+2) = Wr(VP_Rho)*Wr(VP_scal0+2)*scal2
   	fxchapo(VC_Masse) = Wr(VP_Rho)*scal2
   	fxchapo(VC_QDM0+1: VC_QDM0+N_dim)  = (Wr(VP_Rho)*scal2*Wr(VP_U0+1:VP_U0+N_dim)+Wr(VP_P)*norm(:))
   	fxchapo(VC_E)     = (Wr(VP_Rho)*Etr+Wr(VP_P))*scal2*/
  pRight=pOut;

   ROPTM       
    (                 
    rhoRight,
    pRight,
    TOut,         
    YOut,
    N,
    Param_Eos
   );

  YRight=YOut;

   ETOPM      
( 
    interERight,
    TOut,
    pRight,         
    YRight,
    N,
    Param_Eos
);

   URight=ULeft-2.0*scal*normalVector;  //????????????????

scalar rhoERight = rhoRight*interERight+ (rhoRight*URight&URight)/2.0; 


 URelRight = URight - dotX;

 qRight = (URelRight & normalVector);

rhoURight = rhoRight*URight;


//convectionSpeed=qRight;
 convectionSpeed=scal2;

rhoState=rhoRight;

rhoUState=rhoURight;

pState=pRight;

rhoEState=rhoERight;

    for(label j=0; j<N; j++){ 

    rhoYRight[j]=rhoRight*YRight[j];
 //   rhoYRight[j]=rhoRight*YRight[j];

   /* rhoYLeft[j]=rhoown*Yown[j];
    rhoYRight[j]=rhonei*Ynei[j];*/

    }

    for(label j=0; j<N; j++){ 

    rhoYState[j]=rhoYRight[j];

    }

}
else{

      /*  Wr(VP_P) = pstar
        Wr(VP_Y0+1:VP_Y0+Nbre_Yk) = Wl(VP_Y0+1:VP_Y0+Nbre_Yk)
        Wr(VP_scal0+1:VP_scal0+Nbre_scal) = Wl(VP_scal0+1:VP_scal0+Nbre_scal)
        Wr(VP_E) = EROPM(Wr(VP_Rho),Wr(VP_P),Wr(VP_Y0+1:VP_Y0+Nbre_Yk))

 
           Etr=Wr(VP_E)+0.5_DP*SUM(Wr(VP_U0+1:VP_U0+N_dim)**2)

  	fxchapo(VC_roY0+1:VC_roY0+Nbre_Yk) = Wr(VP_Rho)*Wr(VP_Y0+1:VP_Y0+Nbre_Yk)*scal2
        fxchapo(VC_roscal0+1) = Wr(VP_scal0+1)*scal2
        fxchapo(VC_roscal0+2) = Wr(VP_Rho)*Wr(VP_scal0+2)*scal2
   	fxchapo(VC_Masse) = Wr(VP_Rho)*scal2
   	fxchapo(VC_QDM0+1: VC_QDM0+N_dim)  = (Wr(VP_Rho)*scal2*Wr(VP_U0+1:VP_U0+N_dim)+Wr(VP_P)*norm(:))
   	fxchapo(VC_E)     = (Wr(VP_Rho)*Etr+Wr(VP_P))*scal2*/

//rhoRight, and URight are calculated at the beginning


pRight=pOut;
YRight=YLeft;

EROPM        //Internal energy according to the density and pressure of the fluid
(                 //see IJNMF 2016  Eq8, 9
     interERight,
    rhoRight,
    pRight,         
    YRight,
    N,
  Param_Eos
);

scalar rhoERight = rhoRight*interERight+ (rhoRight*URight&URight)/2.0; 


 /// URight=ULeft;  ///?????????????????????????????

 URelRight = URight - dotX;

 qRight = (URelRight & normalVector);

rhoURight = rhoRight*URight;


convectionSpeed=scal2;

rhoState=rhoRight;

rhoUState=rhoURight;

pState=pRight;

rhoEState=rhoERight;

    for(label j=0; j<N; j++){ 

    rhoYRight[j]=rhoRight*YRight[j];
 //   rhoYRight[j]=rhoRight*YRight[j];

   /* rhoYLeft[j]=rhoown*Yown[j];
    rhoYRight[j]=rhonei*Ynei[j];*/

    }

    for(label j=0; j<N; j++){ 

    rhoYState[j]=rhoYRight[j];

    }

}


  



    rhoFlux  = (convectionSpeed*rhoState)*magSf;
    
    rhoUFlux = (convectionSpeed*rhoUState+pState*normalVector)*magSf;
    rhoEFlux = (convectionSpeed*(rhoEState+pState)+pState*(dotX & normalVector))*magSf;  


   for(label j=0; j<N; j++){ 

    rhoYFluxFace[j]=(convectionSpeed*rhoYState[j])*magSf;

    }





                 
}

