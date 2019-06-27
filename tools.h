#include <fstream>
#include "string.h"

void obtenerDatos(istream &file,int nlines,int n,int mode,item* item_list){
    string line;
    file >> line;
    if(nlines==DOUBLELINE) file >> line;

    for(int i=0;i<n;i++){
        switch(mode){
        case INT_FLOAT:
            int e0; float r0;
            file >> e0 >> r0;
            item_list[i].setValues(0,0,0,e0,0,0,r0); //solo valores en el nodo 1
            break;
        case INT_FLOAT_FLOAT:
            int e; float r,rr,rrr;
            file >> e >> r >> rr >> rrr;
            item_list[i].setValues(e,r,rr,rrr,0,0,0,0); //id, x, y, z. rrr es z 
            break;
        case INT_INT_INT_INT:
            int e1,e2,e3,e4,e5;
            file >> e1 >> e2 >> e3 >> e4 >> e5 ;
            item_list[i].setValues(e1,0,0,e2,e3,e4,,e5,0); //id, nodo1, nodo2, nodo 3, e5 es el nodo 4
            break;
        }
    }
}

void correctConditions(int n,condition *list,int *indices){
	//n es 6, list es lo que esta de la linea 100 a la 105 en el archivo, e indices es la posicion en el arreglo de cada uno
	//por ejemplo la condicion con valor 29 tiene indice[0]
    for(int i=0;i<n;i++)
        indices[i] = list[i].getNode1(); //asignamos el valor del nodo 1 a su indice correspondiente

    for(int i=0;i<n-1;i++){ 
        int pivot = list[i].getNode1();
        for(int j=i;j<n;j++)
            //Si la condici�n actual corresponde a un nodo posterior al nodo eliminado por
            //aplicar la condici�n anterior, se debe actualizar su posici�n.
            if(list[j].getNode1()>pivot)
                list[j].setNode1(list[j].getNode1()-1);
    }
}

void addExtension(char *newfilename,char *filename,char *extension){
    int ori_length = strlen(filename);
    int ext_length = strlen(extension);
    int i;
    for(i=0;i<ori_length;i++)
        newfilename[i] = filename[i];
    for(i=0;i<ext_length;i++)
        newfilename[ori_length+i] = extension[i];
    newfilename[ori_length+i] = '\0';
}

void leerMallayCondiciones(mesh &m,char *filename){ //nombre del proyecto en github, el objeto mesh, 
    char inputfilename[150];
    ifstream file;
    float k,Q; //
    int nnodes,neltos,ndirich,nneu; //numero de elementos de cada uno, dependen del archivo

    /*do{
        cout << "Ingrese el nombre del archivo que contiene los datos de la malla: ";
        cin >> filename;
        addExtension(inputfilename,filename,".dat");
        file.open(inputfilename);
    }while(!file);*/

    addExtension(inputfilename,filename,".dat"); //al nombre del proyecto en gid se le a;ade el punto dat
    file.open(inputfilename);

    file >> k >> Q; //primero las condiciones de K y Q, son la primer linea del archivo prueb2.dat
    //cout << "k y Q: "<<k<<" y "<<Q<<"\n";
    file >> nnodes >> neltos >> ndirich >> nneu; //segunda linea de condiciones, son los 4 datos que contienen el tama;o de 
    //las coordenadas(son 37) los elementos(son 52) 6 condiciones de dirichlet y 6 de neumann
    //cout << "sizes: "<<nnodes<<" y "<<neltos<<" y "<<ndirich<<" y "<<nneu<<"\n";

    m.setParameters(k,Q); // m es la instancia de una clase, que que se llama mesh, aqui envio [0.5] y [15] y en esa clase se guardan como arreglos
    m.setSizes(nnodes,neltos,ndirich,nneu); // en la clase mesh hay un arreglo llamado sizes con tama;o 4, 
	//perfecto para guardar los tama;os de las condiciones
    m.createData(); // me crea 5 arreglos de clase: element, nod, 2 tipo condition y uno tipo int para los indices de dirichlet

    obtenerDatos(file,SINGLELINE,nnodes,INT_FLOAT_FLOAT,m.getNodes()); //leo el archivo desde la linea 5 a la 41
    obtenerDatos(file,DOUBLELINE,neltos,INT_INT_INT_INT,m.getElements()); //leo archivo desde la 45 a la 96
    obtenerDatos(file,DOUBLELINE,ndirich,INT_FLOAT,m.getDirichlet()); // leo archivo de la linea 100 a la 105
    obtenerDatos(file,DOUBLELINE,nneu,INT_FLOAT,m.getNeumann()); // leo archivo de la linea 109 a la 114

    file.close();

    correctConditions(ndirich,m.getDirichlet(),m.getDirichletIndices());
}

bool findIndex(int v, int s, int *arr){
    for(int i=0;i<s;i++)
        if(arr[i]==v) return true;
    return false;
}

void writeResults(mesh m,Vector T,char *filename){ //toda la respuesta queda guardada en T 
    char outputfilename[150];
    int *dirich_indices = m.getDirichletIndices();
    condition *dirich = m.getDirichlet();
    ofstream file;

    addExtension(outputfilename,filename,".post.res");
    file.open(outputfilename);

    file << "GiD Post Results File 1.0\n";
    file << "Result \"Temperature\" \"Load Case 1\" 1 Scalar OnNodes\nComponentNames \"T\"\nValues\n";

    int Tpos = 0;
    int Dpos = 0;
    int n = m.getSize(NODES);
    for(int i=0;i<n;i++){
        if(findIndex(i+1,n,dirich_indices)){ 
            file << i+1 << " " << dirich[Dpos].getValue() << "\n";
            Dpos++;
        }else{
            file << i+1 << " " << T.at(Tpos) << "\n";
            Tpos++;
        }
    }

    file << "End values\n";

    file.close();
}
