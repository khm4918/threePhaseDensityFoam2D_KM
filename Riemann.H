          



            forAll(rhoFlux,faceI)
            {    
                 label own = owner[faceI];
                label nei = neigh[faceI];

              double qLYLeftFace[N];
              double qlYRghtFace[N];
              
              
              /*qlYLeft[0][faceI]=qlYLeft0[faceI];
              qlYRght[0][faceI]=qlYRght0[faceI];
              qlYLeft[1][faceI]=qlYLeft1[faceI];
              qlYRght[1][faceI]=qlYRght1[faceI];              
              qlYLeft[2][faceI]=qlYLeft2[faceI];
              qlYRght[2][faceI]=qlYRght2[faceI];  */            
              
              for(label j=0; j<N; j++){ 
     
             //qLYLeftFace[j]=qlYLeft[j][faceI];
             //qlYRghtFace[j]=qlYRght[j][faceI];

            // qLYLeftFace[j]=min(max(qlYLeft[j][faceI],1e-08),1.0);
            // qlYRghtFace[j]=min(max(qlYRght[j][faceI],1e-08),1.0);
              }


             qLYLeftFace[0]=qlYLeft0[faceI];
             qlYRghtFace[0]=qlYRght0[faceI];
             qLYLeftFace[1]=qlYLeft1[faceI];
             qlYRghtFace[1]=qlYRght1[faceI];
             qLYLeftFace[2]=qlYLeft2[faceI];
             qlYRghtFace[2]=qlYRght2[faceI];

              double qLYown[N];
              double qlYnei[N];
              
              for(label j=0; j<N; j++){ 
     
              qLYown[j]=Y[j][own];
              qlYnei[j]=Y[j][nei];

              }
              
             double rhoYFluxFace[N]; 
             
           //  Info<<"rhoLeft "<< qlRLeft[faceI] << "rhoown "<< rho[own]<< endl;
           
          /* Info<<"pLeft "<< qlPLeft[faceI] << "pown "<< p[own]<< endl;
           Info<<"pRght "<< qlPRght[faceI] << "pnei "<< p[nei]<< endl;
           
           Info<<"TLeft "<< qlTLeft[faceI] << "Town "<< T[own]<< endl;
           Info<<"TRght "<< qlTRght[faceI] << "Tnei "<< T[nei]<< endl;

           Info<<"ULeft "<< qlULeft[faceI] << "Uown "<< U[own]<< endl;
           Info<<"URght "<< qlURght[faceI] << "Unei "<< U[nei]<< endl;
           
           Info<<"Y0Left "<< qlYLeft[0][faceI] << "Y0own "<< qLYown[0]<< endl;
           Info<<"Y0Rght "<< qlYRght[0][faceI] << "Y0nei "<< qlYnei[0]<< endl;
           
           */
           
          // Info<<"Y1Left "<< qlYLeft[1][faceI] << "Y1own "<< qLYown[1]<< endl;
          // Info<<"Y2Left "<< qlYLeft[2][faceI] << "Y2own "<< qLYown[2]<< endl;
                
              evaluateFlux
                (
                    rhoFlux[faceI],
                    rhoUFlux[faceI],
                    rhoEFlux[faceI],
                    rhoYFluxFace,                     /////////////////////////////////////////////////             
                    max(qlPLeft[faceI],SMALL),
                   max(qlPRght[faceI],SMALL),
                    qlULeft[faceI],
                    qlURght[faceI],
                    max(qlRLeft[faceI],SMALL),
                    max(qlRRght[faceI],SMALL),      
                    max(qlTLeft[faceI],SMALL),
                    max(qlTRght[faceI],SMALL),            
                    qLYLeftFace,                      /////////////////////////////////////////////////
                    qlYRghtFace,                      /////////////////////////////////////////////////
                    p[own],
                    p[nei],
                    U[own],
                    U[nei],
                    rho[own],
                    rho[nei],     
                    T[own],
                    T[nei],           
                    qLYown,                      /////////////////////////////////////////////////
                    qlYnei, 
                    gamma.value(),
                    gamma.value(),
                    mesh.Sf()[faceI],      
                    mesh.magSf()[faceI],   
                    vector::zero,
                    0.05,
                    N,
                    Param_Eos
                );  

              for(label j=0; j<N; j++){ 
     
              rhoYFlux[j][faceI]=rhoYFluxFace[j];

              }

                            

            }
 
 
           p.correctBoundaryConditions();
           T.correctBoundaryConditions();
           U.correctBoundaryConditions();
           rho.correctBoundaryConditions();

             for(label j=0; j<N; j++){ 
             rhoY[j].correctBoundaryConditions();
            }
         //  rhoca.correctBoundaryConditions();
           
           
 

            forAll(p.boundaryField(), patchI)
            {
                 fvsPatchScalarField& pRhoFlux = rhoFlux.boundaryFieldRef()[patchI];
                                
                 fvsPatchVectorField& pRhoUFlux = rhoUFlux.boundaryFieldRef()[patchI];
                 fvsPatchScalarField& pRhoEFlux = rhoEFlux.boundaryFieldRef()[patchI];
                
                double rhoYFluxFace[N]; 
                
                
                const fvPatchScalarField& pp = p.boundaryField()[patchI];
                const fvPatchScalarField& pE = E.boundaryField()[patchI];
                const fvPatchVectorField& pM = M.boundaryField()[patchI];
                const fvPatchVectorField& pU = U.boundaryField()[patchI];
                const fvPatchScalarField& prho = rho.boundaryField()[patchI];
                
            //    const fvPatchScalarField& prhoca = rhoca.boundaryField()[patchI];
        const fvsPatchScalarField&   qlRLeftpf = qlRLeft.boundaryField()[patchI];
	const fvsPatchVectorField&   qlULeftpf = qlULeft.boundaryField()[patchI];
	const fvsPatchScalarField&   qlPLeftpf = qlPLeft.boundaryField()[patchI];
	const fvsPatchScalarField&   qlTLeftpf = qlTLeft.boundaryField()[patchI];
	
        const fvsPatchScalarField&   qlRRghtpf = qlRRght.boundaryField()[patchI];
	const fvsPatchVectorField&   qlURghtpf = qlURght.boundaryField()[patchI];
	const fvsPatchScalarField&   qlPRghtpf = qlPRght.boundaryField()[patchI];
	const fvsPatchScalarField&   qlTRghtpf = qlTRght.boundaryField()[patchI];	
	
        const fvsPatchScalarField&   qlYLeft0pf = qlYLeft0.boundaryField()[patchI];   
        const fvsPatchScalarField&   qlYLeft1pf = qlYLeft1.boundaryField()[patchI];
        const fvsPatchScalarField&   qlYLeft2pf = qlYLeft2.boundaryField()[patchI];
        
        const fvsPatchScalarField&   qlYRght0pf = qlYRght0.boundaryField()[patchI];   
        const fvsPatchScalarField&   qlYRght1pf = qlYRght1.boundaryField()[patchI];
        const fvsPatchScalarField&   qlYRght2pf = qlYRght2.boundaryField()[patchI];        
                
              //  const fvPatchScalarField& ptotalH = totalH.boundaryField()[patchI];
                const fvsPatchVectorField& pSf = mesh.Sf().boundaryField()[patchI];
                const fvsPatchScalarField& pMagSf = mesh.magSf().boundaryField()[patchI];

               if (pp.coupled())
                {

                   forAll(pp,pfaceI)
                    {
                        const coupledPolyPatch& couplePatch(refCast<const coupledPolyPatch>(mesh.boundaryMesh()[patchI]));
                        label faceI=mesh.boundary()[patchI].start()+pfaceI;
                  

                      // label faceI=mesh.boundary()[patchI].start()+pfaceI;
                    
                                  label own = owner[faceI];
                                  label nei = owner[faceI];   /////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              double qLYLeftFace[N];
              double qlYRghtFace[N];
 
 
             /* qlYLeft[0][faceI]=qlYLeft0[faceI];
              qlYRght[0][faceI]=qlYRght0[faceI];
              qlYLeft[1][faceI]=qlYLeft1[faceI];
              qlYRght[1][faceI]=qlYRght1[faceI];              
              qlYLeft[2][faceI]=qlYLeft2[faceI];
              qlYRght[2][faceI]=qlYRght2[faceI]; */
                           
              for(label j=0; j<N; j++){ 
     
             // qLYLeftFace[j]=qlYLeft[j][faceI];
             // qlYRghtFace[j]=qlYRght[j][faceI];

              }
      
           /*  qLYLeftFace[0]=qlYLeft0[faceI];
             qlYRghtFace[0]=qlYRght0[faceI];
             qLYLeftFace[1]=qlYLeft1[faceI];
             qlYRghtFace[1]=qlYRght1[faceI];
             qLYLeftFace[2]=qlYLeft2[faceI];
             qlYRghtFace[2]=qlYRght2[faceI];*/
  
             qLYLeftFace[0]=qlYLeft0pf[pfaceI];
             qlYRghtFace[0]=qlYRght0pf[pfaceI];
             qLYLeftFace[1]=qlYLeft1pf[pfaceI];
             qlYRghtFace[1]=qlYRght1pf[pfaceI];
             qLYLeftFace[2]=qlYLeft2pf[pfaceI];
             qlYRghtFace[2]=qlYRght2pf[pfaceI];
              

              double qLYown[N];
              double qlYnei[N];
              
              for(label j=0; j<N; j++){ 
     
              qLYown[j]=Y[j][own];
              qlYnei[j]=Y[j][nei];

              }

                           evaluateFlux
                            (
                                pRhoFlux[pfaceI],
                                pRhoUFlux[pfaceI],
                                pRhoEFlux[pfaceI],
                           rhoYFluxFace,
                                
                   max(qlPLeftpf[pfaceI],SMALL),
                   max(qlPRghtpf[pfaceI],SMALL),
                    qlULeftpf[pfaceI],
                    qlURghtpf[pfaceI],
                    max(qlRLeftpf[pfaceI],SMALL),
                    max(qlRRghtpf[pfaceI],SMALL),
                    max(qlTLeftpf[pfaceI],SMALL),
                    max(qlTRghtpf[pfaceI],SMALL),
                    
                    qLYLeftFace,                      /////////////////////////////////////////////////
                    qlYRghtFace,  
                    p[own],
                    p[nei],
                    U[own],
                    U[nei],
                    rho[own],
                    rho[nei],     
                    T[own],
                    T[nei],           
                    qLYown,                      /////////////////////////////////////////////////
                    qlYnei,                         

                                gamma.value(),
                                gamma.value(),
                                pSf[pfaceI],
                                pMagSf[pfaceI],
                                vector::zero,
                                0.5,
                                N,
                                Param_Eos
                            );

              for(label j=0; j<N; j++){     
              rhoYFlux[j].boundaryFieldRef()[patchI][pfaceI]=rhoYFluxFace[j];
              }

                      }
                  
                          
                  
                }
                else
                {

                forAll(pp,pfaceI)
                    {
                                           // const coupledPolyPatch& couplePatch(refCast<const coupledPolyPatch>(mesh.boundaryMesh()[patchI]));
                       label faceI=mesh.boundary()[patchI].start()+pfaceI;
                    
                                  label own = owner[faceI];
                                  
                                  
                                  vector Ub;

              double qLYLeftFace[N];
              double qlYRghtFace[N];
              
              for(label j=0; j<N; j++){ 
     
              qLYLeftFace[j]=Y[j][own];
              qlYRghtFace[j]=Y[j][own];

              }

                              vector uLeft=U[own];
                              vector uRight=U[own];

               /* if(U.boundaryField().types()[patchI]=="symmetryPlane"){

                       vector nHat=pSf[pfaceI]/pMagSf[pfaceI];
                    //   qlMLeft[faceI]=pM.patchInternalField()()[pfaceI];
                    //   qlMRght[faceI]=transform(I-2.0*sqr(nHat),qlMLeft[faceI]);
                    //   qlMLeft[faceI]=M[own];
                    //   qlMRght[faceI]=transform(I-2.0*sqr(nHat),qlMLeft[faceI]);
                       uLeft=U[own];
                     // uRight=transform(I-2.0*sqr(nHat),uLeft);
                       uRight=vector(-uLeft.x(), uLeft.y(),0.0);
                                              }*/
               //  if(U.boundaryField().types()[patchI]=="symmetryPlane"){

                  if(mesh.boundary()[patchI].name()=="outletP"){

// Info<<"find"<<endl;                              
                              uLeft=U[own];
                            uRight=U[own];

                      evaluateFluxOutlet
                        (
                            pRhoFlux[pfaceI],
                            pRhoUFlux[pfaceI],
                            pRhoEFlux[pfaceI],
                            rhoYFluxFace,
                            
  
                    p[own],
                    p[own],
                    uLeft,
                    uRight,
                    rho[own],
                    rho[own],
                    T[own],
                    T[own],
                 qLYLeftFace,
                 qlYRghtFace,  

                   p[own],
                    p[own],
                    uLeft,
                    uRight,
                    rho[own],
                    rho[own],
                    T[own],
                    T[own],
                 qLYLeftFace,
                 qlYRghtFace,  
                         
                            gamma.value(),
                            gamma.value(),
                            pSf[pfaceI],
                            pMagSf[pfaceI],
                            vector::zero,
                            0.5,
                            N,
                            Param_Eos
                        );

                   }
                else if(U.boundaryField().types()[patchI]=="symmetryPlane" || U.boundaryField().types()[patchI]=="slip"){

                      vector nHat=pSf[pfaceI]/pMagSf[pfaceI];

                      uLeft=U[own];
                      uRight=transform(I-2.0*sqr(nHat),uLeft);

                     evaluateFlux
                        (
                            pRhoFlux[pfaceI],
                            pRhoUFlux[pfaceI],
                            pRhoEFlux[pfaceI],
                            rhoYFluxFace,
                            
  
                    p[own],
                    p[own],
                    uLeft,
                    uRight,
                    rho[own],
                    rho[own],
                    T[own],
                    T[own],
                 qLYLeftFace,
                 qlYRghtFace,  

                   p[own],
                    p[own],
                    uLeft,
                    uRight,
                    rho[own],
                    rho[own],
                    T[own],
                    T[own],
                 qLYLeftFace,
                 qlYRghtFace,  
                         
                            gamma.value(),
                            gamma.value(),
                            pSf[pfaceI],
                            pMagSf[pfaceI],
                            vector::zero,
                            0.5,
                            N,
                            Param_Eos
                        );


                  }
                  else{

                    //pp[pfaceI],
                  //  pp[pfaceI],
    
                  //  pU[pfaceI],
                   // pU[pfaceI],
               uLeft=U[own];
                            uRight=U[own];

                    evaluateFlux
                        (
                            pRhoFlux[pfaceI],
                            pRhoUFlux[pfaceI],
                            pRhoEFlux[pfaceI],
                            rhoYFluxFace,
                            
  
                    p[own],
                    p[own],
                    uLeft,
                    uRight,
                    rho[own],
                    rho[own],
                    T[own],
                    T[own],
                 qLYLeftFace,
                 qlYRghtFace,  

                   p[own],
                    p[own],
                    uLeft,
                    uRight,
                    rho[own],
                    rho[own],
                    T[own],
                    T[own],
                 qLYLeftFace,
                 qlYRghtFace,  
                         
                            gamma.value(),
                            gamma.value(),
                            pSf[pfaceI],
                            pMagSf[pfaceI],
                            vector::zero,
                            0.5,
                            N,
                            Param_Eos
                        );
         
                   }
                        
                       //   revise it later
                                     for(label j=0; j<N; j++){     
              rhoYFlux[j].boundaryFieldRef()[patchI][pfaceI]=rhoYFluxFace[j];
              }
     
                        
                    }
      
              
            }



          }
          
          

           Info<<"Rieman"<<endl;
