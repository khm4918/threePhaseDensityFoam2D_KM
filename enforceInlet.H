       //inlet condition //////////////////////////////////////////////////////
forAll(p,i){
    scalar x=mesh.C()[i].x();
    scalar y=mesh.C()[i].y();  
 
        if(x<=-0.019 && y>=0.0025 ){
        T[i]=150.0;
        U[i].x()=200.0;
        U[i].y()=0.0;
	p[i]=3.0e6;
        }

        if(x<=-0.019 && y<=-0.0025 ){
        T[i]=150.0;
        U[i].x()=200.0;
        U[i].y()=0.0;
	p[i]=3.0e6;
        }

        if(x<=-0.019 && y>=-0.00175 && y<=0.00175 ){
        T[i]=100.0;
        U[i].x()=30.0;
        U[i].y()=0.0;
	p[i]=3.0e6;
        }

                     }

////////////////////////////////////////////Y

forAll(Y,i){
                   // volScalarField& Yi = Y[i];
   forAll(p,cell){
    scalar x=mesh.C()[cell].x();
    scalar y=mesh.C()[cell].y();  
        


       if(x<=-0.019 && y>=0.0025 ){
        if(i==0){Y[i][cell]=1e-8;}
        if(i==1){Y[i][cell]=1e-8;}
        if(i==2){Y[i][cell]=0.99999998;}
     //   if(i==3){Y[i][cell]=1e-8;}
        }

        if(x<=-0.019 && y<=-0.0025 ){
        if(i==0){Y[i][cell]=1e-8;}
        if(i==1){Y[i][cell]=1e-8;}
        if(i==2){Y[i][cell]=0.99999998;}
     //   if(i==3){Y[i][cell]=1e-8;}
        }

        if(x<=-0.019 && y>=-0.00175 && y<=0.00175 ){
        if(i==0){Y[i][cell]=0.99999998;}
        if(i==1){Y[i][cell]=1e-8;}
        if(i==2){Y[i][cell]=1e-8;}
     //   if(i==3){Y[i][cell]=1e-8;}
        }



                     }
                    

                    }


                forAll(rho,cellI){   //caculate rho from T and P

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

                forAll(Y,i){

             rhoY[i]=Y[i]*rho;
                    }

        // T.correctBoundaryConditions();
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
            rhoPhi = rhoFlux;
  
    M = rho*U;

         //   totalH=E/rho+p/rho;
            E.correctBoundaryConditions();
  //inlet condition //////////////////////////////////////////////////////
