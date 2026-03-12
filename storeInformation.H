  label noe=mesh.nCells();
    labelListList mixFaceLabels(noe);
    labelListList mixVertLabels(noe);
    labelList mixMeshType(noe);

    hexMatcher hex;
    prismMatcher prism;
    cellShape shape;
    const pointField& faceCentres = mesh.faceCentres();
    forAll(mixFaceLabels,cellI)
    {
       bool isHex = hex.matches(mesh,cellI, shape);
       bool isPrism = prism.matches(mesh,cellI, shape);
       if(isHex)
       {
          mixMeshType[cellI]=0;
          label faceNum=hex.faceLabels().size();
          label vertNum=hex.vertLabels().size();
          mixFaceLabels[cellI].setSize(faceNum);
          mixVertLabels[cellI].setSize(vertNum);
          forAll(hex.faceLabels(),cfI) mixFaceLabels[cellI][cfI]=hex.faceLabels()[cfI];
          forAll(hex.vertLabels(),cpI) mixVertLabels[cellI][cpI]=hex.vertLabels()[cpI];

          label pt[8],fe[6];
          forAll(mixVertLabels[cellI],pI) pt[pI]=mixVertLabels[cellI][pI];
          forAll(mixFaceLabels[cellI],fI) fe[fI]=mixFaceLabels[cellI][fI];

          label Face[2]={4,5};
          label i=0;
          for(label m=0; m<6; m++)
          {
             if(mag(faceCentres[fe[m]].z()-mesh.C()[cellI].z())>SMALL)
             {
                if(i>1) 
                {
                   Info<<"wrong face arrangement or not in X-Y plane!!"<<endl;
                   exit(0);
                }
                Face[i++]=m;
             }
          }
    
          if(Face[0]==0 && Face[1]==1)
          {
             mixVertLabels[cellI][0]=pt[0];
             mixVertLabels[cellI][1]=pt[3];
             mixVertLabels[cellI][2]=pt[7];
             mixVertLabels[cellI][3]=pt[4];
             mixVertLabels[cellI][4]=pt[1];
             mixVertLabels[cellI][5]=pt[2];
             mixVertLabels[cellI][6]=pt[6];
             mixVertLabels[cellI][7]=pt[5];
             mixFaceLabels[cellI][0]=fe[2];
             mixFaceLabels[cellI][1]=fe[3];
             mixFaceLabels[cellI][2]=fe[4];
             mixFaceLabels[cellI][3]=fe[5];
             mixFaceLabels[cellI][4]=fe[0];
             mixFaceLabels[cellI][5]=fe[1];
          }
          else if(Face[0]==2 && Face[1]==3)
          {
             mixVertLabels[cellI][0]=pt[0];
             mixVertLabels[cellI][1]=pt[1];
             mixVertLabels[cellI][2]=pt[5];
             mixVertLabels[cellI][3]=pt[4];
             mixVertLabels[cellI][4]=pt[3];
             mixVertLabels[cellI][5]=pt[2];
             mixVertLabels[cellI][6]=pt[6];
             mixVertLabels[cellI][7]=pt[7];
             mixFaceLabels[cellI][0]=fe[0];
             mixFaceLabels[cellI][1]=fe[1];
             mixFaceLabels[cellI][2]=fe[4];
             mixFaceLabels[cellI][3]=fe[5];
             mixFaceLabels[cellI][4]=fe[2];
             mixFaceLabels[cellI][5]=fe[3];
          }
       }
       else if(isPrism)
       {
          mixMeshType[cellI]=1;
          label faceNum=prism.faceLabels().size();
          label vertNum=prism.vertLabels().size();
          mixFaceLabels[cellI].setSize(faceNum);
          mixVertLabels[cellI].setSize(vertNum);
          forAll(prism.faceLabels(),cfI) mixFaceLabels[cellI][cfI]=prism.faceLabels()[cfI];
          forAll(prism.vertLabels(),cpI) mixVertLabels[cellI][cpI]=prism.vertLabels()[cpI];

          label pt[6],fe[5];
          forAll(mixVertLabels[cellI],pI) pt[pI]=mixVertLabels[cellI][pI];
          forAll(mixFaceLabels[cellI],fI) fe[fI]=mixFaceLabels[cellI][fI];

          mixVertLabels[cellI][0]=pt[2];
          mixVertLabels[cellI][1]=pt[0];
          mixVertLabels[cellI][2]=pt[1];
          mixVertLabels[cellI][3]=pt[5];
          mixVertLabels[cellI][4]=pt[3];
          mixVertLabels[cellI][5]=pt[4];
          mixFaceLabels[cellI][0]=fe[2];
          mixFaceLabels[cellI][1]=fe[3];
          mixFaceLabels[cellI][2]=fe[4];
          mixFaceLabels[cellI][3]=fe[0];
          mixFaceLabels[cellI][4]=fe[1];
       }
       else
       {
          Info<<"this scheme only support Tri and Qua mesh element for 2d!!"<<endl;
          exit(0);
       }
    }
    

    const labelUList& owner = mesh.owner();
    const labelUList& neigh = mesh.neighbour();

    //Boundary Point_Search
    boolList isPatchPoint_(mesh.nPoints());
    isPatchPoint_ = false;

    boolList isPatchFace_(mesh.nFaces());
    isPatchFace_ = false;

    boolList isPatchCell_(mesh.nCells());
    isPatchCell_ = false;

    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& meshPatch = mesh.boundaryMesh()[patchI];

        if(!isA<emptyPolyPatch>(meshPatch) && !mesh.magSf().boundaryField()[patchI].coupled()) 
        {        
            const labelList& pfaceCells = mesh.boundary()[patchI].faceCells();         
            forAll(meshPatch, pfaceI)
            {
                isPatchCell_[pfaceCells[pfaceI]] = true;
                label faceI = mesh.boundary()[patchI].start()+pfaceI;
                isPatchFace_[faceI] = true;
                const face& facePoint=mesh.faces()[faceI];
                forAll(facePoint, facePointI) isPatchPoint_[facePoint[facePointI]] = true;
            }
        }
    }
    syncTools::syncPointList(mesh, isPatchPoint_, orEqOp<bool>(),false);

  /*  typedef List<vectorList> vectorListList;
    vectorListList Ma(noe);
    vectorListList MA(noe);   
    forAll(mesh.cells(),cellI)
    {
        label num=mixVertLabels[cellI].size();
        Ma[cellI].setSize(num);
        MA[cellI].setSize(num);
    }
    forAll(mesh.cells(),cellI)
    { 
        label pt[8];
        forAll(mixVertLabels[cellI],pI) pt[pI]=mixVertLabels[cellI][pI];

        if(mixMeshType[cellI]==0)
        {
            scalar xM[4],yM[4];
            for(label m=0;m<4;m++)
            {
                xM[m]=mesh.points()[pt[m]].x(); 
                yM[m]=mesh.points()[pt[m]].y();
            }

            scalar wa[8],wb[8];
            wa[0]=wa[4]=-1.0;
            wa[1]=wa[5]= 1.0;
            wa[2]=wa[6]= 1.0;
            wa[3]=wa[7]=-1.0;
            wb[0]=wb[4]=-1.0;
            wb[1]=wb[5]=-1.0;
            wb[2]=wb[6]= 1.0;
            wb[3]=wb[7]= 1.0;

            for(int n=0; n<8; n++)
            {
                Ma[cellI][n][0]=(-xM[0]+xM[1]+xM[2]-xM[3]+(xM[0]-xM[1]+xM[2]-xM[3])*wb[n])/4.0;
                Ma[cellI][n][1]=(-yM[0]+yM[1]+yM[2]-yM[3]+(yM[0]-yM[1]+yM[2]-yM[3])*wb[n])/4.0;
                MA[cellI][n][0]=(-xM[0]-xM[1]+xM[2]+xM[3]+(xM[0]-xM[1]+xM[2]-xM[3])*wa[n])/4.0;
                MA[cellI][n][1]=(-yM[0]-yM[1]+yM[2]+yM[3]+(yM[0]-yM[1]+yM[2]-yM[3])*wa[n])/4.0;
            }
        }
        else
        {
            scalar xM[3],yM[3];
            for(label m=0;m<3;m++)
            {
                xM[m]=mesh.points()[pt[m]].x(); 
                yM[m]=mesh.points()[pt[m]].y();
            }

            for(int n=0; n<6; n++)
            {
                Ma[cellI][n][0]=xM[1]-xM[0]; //sx
                Ma[cellI][n][1]=yM[1]-yM[0]; //sy
                MA[cellI][n][0]=xM[2]-xM[0]; //tx
                MA[cellI][n][1]=yM[2]-yM[0]; //ty
            } 
        }
    }

    tensorField l6(mesh.nPoints(),tensor::zero);
    forAll(mesh.cells(),cellI)
    {
        labelList& pt=mixVertLabels[cellI];
 
        for(label m=0; m<pt.size(); m++)
        {
            vector dd=mesh.C()[cellI]-mesh.points()[pt[m]];
            l6[pt[m]][0]+=1.0; //l11[p[m]]
            l6[pt[m]][1]+=dd[0]; //l12[p[m]] 
            l6[pt[m]][2]+=dd[1]; //l13[p[m]]
            l6[pt[m]][3]+=dd[0]*dd[0]; //l22[p[m]]
            l6[pt[m]][4]+=dd[0]*dd[1]; //l23[p[m]]
            l6[pt[m]][5]+=dd[1]*dd[1]; //l33[p[m]]
        }
    }

    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& meshPatch = mesh.boundaryMesh()[patchI];

        if(!isA<emptyPolyPatch>(meshPatch) && !mesh.magSf().boundaryField()[patchI].coupled()) 
        {     
            forAll(meshPatch, pfaceI)
            {
                label faceI = mesh.boundary()[patchI].start()+pfaceI;
                const face& facePoint=mesh.faces()[faceI];
                forAll(facePoint, facePointI)
                {
                    vector dd=faceCentres[faceI]-mesh.points()[facePoint[facePointI]];
                    l6[facePoint[facePointI]][0]+=1.0; //l11[p[m]]
                    l6[facePoint[facePointI]][1]+=dd[0]; //l12[p[m]]
                    l6[facePoint[facePointI]][2]+=dd[1]; //l13[p[m]]
                    l6[facePoint[facePointI]][3]+=dd[0]*dd[0]; //l22[p[m]]
                    l6[facePoint[facePointI]][4]+=dd[0]*dd[1]; //l23[p[m]]
                    l6[facePoint[facePointI]][5]+=dd[1]*dd[1]; //l33[p[m]]              
                }
            }
        }
    }
    syncTools::syncPointList(mesh, l6, plusEqOp<tensor>(),tensor::zero);

    List<scalarRectangularMatrix> leastSquareSVD(mesh.nPoints());
    forAll(mesh.points(),pointI)
    {       
        scalar l11=l6[pointI][0];
        scalar l12=l6[pointI][1];
        scalar l13=l6[pointI][2];
        scalar l22=l6[pointI][3];
        scalar l23=l6[pointI][4];
        scalar l33=l6[pointI][5];

        scalarRectangularMatrix matrix(3,3);
        scalar dd=2.0*l13*l12*l23-l13*l13*l22-l12*l12*l33+l11*l22*l33-l11*l23*l23;
        matrix[0][0]=(l22*l33-l23*l23)/dd;   //a11;
        matrix[1][1]=(l11*l33-l13*l13)/dd;   //a22
        matrix[2][2]=(l11*l22-l12*l12)/dd;   //a33
        matrix[0][1]=matrix[1][0]=(l13*l23-l12*l33)/dd;   //a12
        matrix[0][2]=matrix[2][0]=(l12*l23-l13*l22)/dd;   //a13
        matrix[1][2]=matrix[2][1]=(l13*l12-l11*l23)/dd;   //a23
  
        leastSquareSVD[pointI]=matrix;
    }


    vectorField cellNeigh(mesh.nFaces(),vector::zero);
    vectorField cellOwner(mesh.nFaces(),vector::zero);
    forAll(mesh.magSf(),faceI)
    {
        label own = owner[faceI];
        label nei = neigh[faceI];
        cellNeigh[faceI]=mesh.C()[nei];
        cellOwner[faceI]=mesh.C()[own];
    }

    forAll(mesh.boundary(), patchI)
    {
        const fvsPatchScalarField& magSfPf = mesh.magSf().boundaryField()[patchI];
        const pointField& faceCentres = mesh.faceCentres();
        const labelList& pfaceCells = mesh.boundary()[patchI].faceCells();
        forAll(magSfPf, pfaceI)
        {
            label pfCellI = pfaceCells[pfaceI];
            label faceI = mesh.boundary()[patchI].start()+pfaceI;
            if (magSfPf.coupled())
            {
                cellNeigh[faceI]=mesh.C()[pfCellI]-faceCentres[faceI];
                cellOwner[faceI]=mesh.C()[pfCellI]-faceCentres[faceI];                      
            }
            else
            {
                cellNeigh[faceI]=faceCentres[faceI];
                cellOwner[faceI]=mesh.C()[pfCellI]; 
            }
        }
    }
    syncTools::swapFaceList(mesh,cellNeigh);

    List<scalarRectangularMatrix> cellGradientInv(mesh.nCells());
    forAll(mesh.cells(),cellI)
    {
        labelList& pt=mixVertLabels[cellI];
        labelList& fe=mixFaceLabels[cellI];
 
        if(mixMeshType[cellI]==0)
        {
            vector cellPoints[8];
            for(label n=0; n<4; n++)
            {
                cellPoints[n]=(mesh.points()[pt[n]]-mesh.C()[cellI]);
                if(mesh.isInternalFace(fe[n]))
                { 
                    if(cellI==owner[fe[n]]) cellPoints[n+4]=(cellNeigh[fe[n]]-cellOwner[fe[n]]);
                    else cellPoints[n+4]=(cellOwner[fe[n]]-cellNeigh[fe[n]]);
                }
                else
                {
                    label iPatch=mesh.boundaryMesh().whichPatch(fe[n]);
                    const surfaceScalarField& magSf = mesh.magSf();
                    if(magSf.boundaryField()[iPatch].coupled()) cellPoints[n+4]=(cellNeigh[fe[n]]-cellOwner[fe[n]]);
                    else cellPoints[n+4]=(faceCentres[fe[n]]-mesh.C()[cellI]);
                }
            }

            scalarRectangularMatrix leastSquare(8,5);
            for(label n=0; n<8; n++)
            {        
               vector dd=cellPoints[n];
               leastSquare[n][0]=dd[0];
               leastSquare[n][1]=dd[1];
               leastSquare[n][2]=dd[0]*dd[0];
               leastSquare[n][3]=dd[0]*dd[1];
               leastSquare[n][4]=dd[1]*dd[1];
            }

            cellGradientInv[cellI]=SVDinv(leastSquare);
        }
        else
        {
            vector cellPoints[6];
            for(label n=0; n<3; n++)
            {
                cellPoints[n]=(mesh.points()[pt[n]]-mesh.C()[cellI]);
                if(mesh.isInternalFace(fe[n]))
                { 
                    if(cellI==owner[fe[n]]) cellPoints[n+3]=(cellNeigh[fe[n]]-cellOwner[fe[n]]);
                    else cellPoints[n+3]=(cellOwner[fe[n]]-cellNeigh[fe[n]]);
                }
                else
                {
                    label iPatch=mesh.boundaryMesh().whichPatch(fe[n]);
                    const surfaceScalarField& magSf = mesh.magSf();
                    if(magSf.boundaryField()[iPatch].coupled()) cellPoints[n+3]=(cellNeigh[fe[n]]-cellOwner[fe[n]]);
                    else cellPoints[n+3]=(faceCentres[fe[n]]-mesh.C()[cellI]);
                }
            }

            scalar matrix[2][6];
            scalar l11,l22,l12;
            l11=l22=l12=0.0;
            for(label n=0; n<6; n++)
            {        
                vector dd=cellPoints[n];
                matrix[0][n]=dd[0];
                matrix[1][n]=dd[1]; 

                l11+=dd[0]*dd[0]; 
                l12+=dd[0]*dd[1]; 
                l22+=dd[1]*dd[1]; 
            }         
     
            scalar symMatrix[2][2];
            scalar det=l11*l22-l12*l12;
            symMatrix[0][0]=l22/det;  
            symMatrix[1][1]=l11/det;
            symMatrix[0][1]=symMatrix[1][0]=-l12/det; 

            scalarRectangularMatrix matrixInv(2,6);
            for(int n=0; n<2; n++)
            {
                for(int m=0; m<6; m++)
                {
                    matrixInv[n][m]=0.0;
                    for(int k=0; k<2; k++)
                    {
                        matrixInv[n][m]+=symMatrix[n][k]*matrix[k][m];
                    }
                }           
            }       
     
            cellGradientInv[cellI]=matrixInv;
        }
    }


      List<scalarRectangularMatrix> fvmGradientInv(mesh.nCells());
    forAll(mesh.cells(),cellI)
    {
        labelList& fe=mixFaceLabels[cellI];

        label numFe=4;
        if(mixMeshType[cellI]==1) numFe=3;
 
        vector cellPoints[numFe];
        for(label n=0; n<numFe; n++)
        {
            if(mesh.isInternalFace(fe[n]))
            { 
                if(cellI==owner[fe[n]]) cellPoints[n]=(cellNeigh[fe[n]]-cellOwner[fe[n]]);
                else cellPoints[n]=(cellOwner[fe[n]]-cellNeigh[fe[n]]);
            }
            else
            {
                label iPatch=mesh.boundaryMesh().whichPatch(fe[n]);
                const surfaceScalarField& magSf = mesh.magSf();
                if(magSf.boundaryField()[iPatch].coupled()) cellPoints[n]=(cellNeigh[fe[n]]-cellOwner[fe[n]]);
                else cellPoints[n]=(cellNeigh[fe[n]]-cellOwner[fe[n]]);
            }
        }

        scalar matrix[2][numFe];
        scalar l11,l22,l12;
        l11=l22=l12=0.0;
        for(label n=0; n<numFe; n++)
        {        
            vector dd=cellPoints[n];
            matrix[0][n]=dd[0];
            matrix[1][n]=dd[1];

            l11+=dd[0]*dd[0]; 
            l12+=dd[0]*dd[1]; 
            l22+=dd[1]*dd[1]; 
        }         
 
        scalar symMatrix[2][2];
        scalar det=l11*l22-l12*l12;
        symMatrix[0][0]=l22/det;  
        symMatrix[1][1]=l11/det;
        symMatrix[0][1]=symMatrix[1][0]=-l12/det; 

        scalarRectangularMatrix matrixInv(2,numFe);
        for(int n=0; n<2; n++)
        {
            for(int m=0; m<numFe; m++)
            {
                matrixInv[n][m]=0.0;
                for(int k=0; k<2; k++)
                {
                    matrixInv[n][m]+=symMatrix[n][k]*matrix[k][m];
                }
            }           
        }       
 
        fvmGradientInv[cellI]=matrixInv;
    }*/


