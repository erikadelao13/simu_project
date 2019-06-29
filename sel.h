void showMatrix(Matrix K){
    for(int i=0;i<K.at(0).size();i++){
        cout << "[\t";
        for(int j=0;j<K.size();j++){
            cout << K.at(i).at(j) << "\t";
        }
        cout << "]\n";
    }
}

void showKs(vector<Matrix> Ks){
    for(int i=0;i<Ks.size();i++){
        cout << "K del elemento "<< i+1 << ":\n";
        showMatrix(Ks.at(i));
        cout << "*************************************\n";
    }
}

void showVector(Vector b){
    cout << "[\t";
    for(int i=0;i<b.size();i++){
        cout << b.at(i) << "\t";
    }
    cout << "]\n";
}

void showbs(vector<Vector> bs){
    for(int i=0;i<bs.size();i++){
        cout << "b del elemento "<< i+1 << ":\n";
        showVector(bs.at(i));
        cout << "*************************************\n";
    }
}

float calculateLocalD(int i,mesh m){
    float D,a,b,c,d;

    element e = m.getElement(i);

    node n1 = m.getNode(e.getNode1()-1); //retorna una lista de nodos valor del nodo menos 1
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

    a=n2.getX()-n1.getX();b=n2.getY()-n1.getY(); //|x2-x1    y2-y1|
    c=n3.getX()-n1.getX();d=n3.getY()-n1.getY(); //|x3 -x1   y3-y1|
    D = a*d - b*c; //determinante

    return D;
}

float calculateMagnitude(float v1, float v2){
    return sqrt(pow(v1,2)+pow(v2,2));
}

float calculateLocalArea(int i,mesh m){
    float A,s,a,b,c;
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1());
    node n2 = m.getNode(e.getNode2());
    node n3 = m.getNode(e.getNode3());

    a = calculateMagnitude(n2.getX()-n1.getX(),n2.getY()-n1.getY()); //raizcuadradade((x2-x1)`2 + (y2-y1)`2)
    b = calculateMagnitude(n3.getX()-n2.getX(),n3.getY()-n2.getY());//raizcuadradade((x3-x2)`2 + (y3-y2)`2)
    c = calculateMagnitude(n3.getX()-n1.getX(),n3.getY()-n1.getY());//raizcuadradade((x3-x1)`2 + (y3-y1)`2)
    s = (a+b+c)/2; //area del triangulo

    A = sqrt(s*(s-a)*(s-b)*(s-c));// equivalente a la integral de de area con -sin
    return A;
}

/*void calculateLocalA(int i,Matrix &A,mesh m){
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1());
    node n2 = m.getNode(e.getNode2());
    node n3 = m.getNode(e.getNode3());
    A.at(0).at(0) = n3.getY()-n1.getY(); A.at(0).at(1) = n1.getY()-n2.getY(); //|y3-y1 y1-y2|
    A.at(1).at(0) = n1.getX()-n3.getX(); A.at(1).at(1) = n2.getX()-n1.getX(); //|x1-x3 x2-x1|
}*/

/*void calculateB(Matrix &B){
    B.at(0).at(0) = -1; B.at(0).at(1) = 1; B.at(0).at(2) = 0; //|-1 1 0|
    B.at(1).at(0) = -1; B.at(1).at(1) = 0; B.at(1).at(2) = 1;// |-1 0 1|
}*/ 
void calculateProductNNt(Matrix &NNt){
    NNt.at(0).at(0) = 1/2; 
    NNt.at(0).at(1) = -1/3; 
    NNt.at(0).at(2) = -1/3;
    NNt.at(0).at(3) = -1/3; 
    NNt.at(1).at(0) = -1/3; 
    NNt.at(1).at(1) = 1/3; 
    NNt.at(1).at(2) = 1/4;
    NNt.at(1).at(3) = 1/4;
    NNt.at(2).at(0) = -1/3; 
    NNt.at(2).at(1) = 1/4; 
    NNt.at(2).at(2) = 1/3;
    NNt.at(2).at(3) = 1/4;
    NNt.at(3).at(0) = -1/3; 
    NNt.at(3).at(1) = 1/4; 
    NNt.at(3).at(2) = 1/4;
    NNt.at(3).at(3) = 1/3;           
}



Matrix createLocalK(int element,mesh &m){
    //float D,Ae,k = m.getParameter(THERMAL_CONDUCTIVITY); //k = 0.5
    Matrix K,A,NNt,Bt,At;
    float k,Th,Tc,L = m.getParameter(THERMAL_CONDUCTIVITY); 


    D = calculateLocalD(element,m); //elemento actual y el objeto mesh det
    Ae = calculateLocalArea(element,m); //Area del elemento

    zeroes(A,2); // da formato 0.0 a los valores de la matriz de area  matriz 2x2
    zeroes(B,2,3);  // da formato de 0.0 a B, matriz 2x3
    //calculateLocalA(element,A,m); //elemento A
    //calculateProductNNt(NNt); //multiplicacion de NNt
    //transpose(A,At); //A transpuesta. |y3-y1 x1-x3|
									//|y1-y2 x2-x1|
    //transpose(B,Bt); //B transpuesta |y3-y1 y1-y2|
    							   //|x1-x3 x2-x1|

    productRealMatrix(calculateLocalJ(element,m)*2*k*(Th-Tc/L),calculateProductNNt(NNt),K);
    //k minuscula es 0.5, K mayuscula es una matriz vacia por el momento, ahi se almacenara el resultado de la multiplicacion

    return K;
}

float calculateLocalJ(int i,mesh m){
    float J,a,b,c,d,k,f,g,h,i;
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1());
    node n2 = m.getNode(e.getNode2());
    node n3 = m.getNode(e.getNode3());
    node n4 = m.getNode(e.getNode4());

    a=n2.getX()-n1.getX();
    b=n3.getX()-n1.getX();
    c=n4.getX()-n1.getX();
    d=n2.getY()-n1.getY();
    k=n3.getY()-n1.getY();
    f=n4.getY()-n1.getY();
    g=n2.getZ()-n1.getZ();
    h=n3.getZ()-n1.getZ();
    i=n4.getZ()-n1.getZ();
    J = a*((k*i) - (h*f)) - b*((d*i) - (g*f)) + c*((d*h) - (g*k));
    return J;
}

/* Vector createLocalb(int element,mesh &m){
    Vector b;

    float Q = m.getParameter(HEAT_SOURCE),J;
    J = calculateLocalJ(element,m);

    b.push_back(0); b.push_back(Q*J*0.5); b.push_back(Q*J*0.5);
    //esto equivale a: QJ[ 0 ]
    //					 |0.5|
    //					 [0.5]

    return b;
} */

void crearSistemasLocales(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs){
    for(int i=0;i<m.getSize(ELEMENTS);i++){ //52 elementos
        localKs.push_back(createLocalK(i,m));//K local
        //localbs.push_back(createLocalb(i,m)); //b Local
    }
}

void assemblyK(element e,Matrix localK,Matrix &K){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = e.getNode3() - 1;
    int index4 = e.getNode4() - 1;

    K.at(index1).at(index1) += localK.at(0).at(0);
    K.at(index1).at(index2) += localK.at(0).at(1);
    K.at(index1).at(index3) += localK.at(0).at(2);
    K.at(index1).at(index4) += localK.at(0).at(3);
    K.at(index2).at(index1) += localK.at(1).at(0);
    K.at(index2).at(index2) += localK.at(1).at(1);
    K.at(index2).at(index3) += localK.at(1).at(2);
    K.at(index2).at(index4) += localK.at(1).at(3);
    K.at(index3).at(index1) += localK.at(2).at(0);
    K.at(index3).at(index2) += localK.at(2).at(1);
    K.at(index3).at(index3) += localK.at(2).at(2);
    K.at(index3).at(index4) += localK.at(2).at(3);
    K.at(index4).at(index1) += localK.at(3).at(0);
    K.at(index4).at(index2) += localK.at(3).at(1);
    K.at(index4).at(index3) += localK.at(3).at(2);
    K.at(index4).at(index4) += localK.at(3).at(3);
}

void assemblyb(element e,Vector localb,Vector &b){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = e.getNode3() - 1;
    int index4 = e.getNode4() - 1;

    b.at(index1) += 0;
    b.at(index2) += 0;
    b.at(index3) += 0;
    b.at(index4) += 0;
}

void ensamblaje(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs,Matrix &K,Vector &b){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        element e = m.getElement(i);
        assemblyK(e,localKs.at(i),K);
        assemblyb(e,localbs.at(i),b);
    }
}

// void applyNeumann(mesh &m,Vector &b){
//     for(int i=0;i<m.getSize(NEUMANN);i++){
//         condition c = m.getCondition(i,NEUMANN);
//         b.at(c.getNode1()-1) += c.getValue();
//     }
// }

void applyDirichlet(mesh &m,Matrix &K,Vector &b){
    for(int i=0;i<m.getSize(DIRICHLET);i++){
        condition c = m.getCondition(i,DIRICHLET);
        int index = c.getNode1()-1;

        K.erase(K.begin()+index);
        b.erase(b.begin()+index);

        for(int row=0;row<K.size();row++){
            float cell = K.at(row).at(index);
            K.at(row).erase(K.at(row).begin()+index);
            b.at(row) += -1*c.getValue()*cell;
        }
    }
}

void calculate(Matrix &K, Vector &b, Vector &T){
    cout << "Iniciando calculo de respuesta...\n";
    Matrix Kinv;
    cout << "Calculo de inversa...\n";
    inverseMatrix(K,Kinv);
    cout << "Calculo de respuesta...\n";
    productMatrixVector(Kinv,b,T);
}
