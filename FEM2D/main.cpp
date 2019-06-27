#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "sel.h"
//por ejemplo si yo quiero tener dos archivos en gid, argc deberia de valer 3, el nombre del ejecutable, y los dos nombres de los otros dos archivos
int main(int argc, char *argv[]) //argc es la cantidad de argumentos que se puserion en la linea de comandos. argv es un conjunto de cadenas
{
    char filename[150]; //argv0 siempre es el nombre del ejecutable
    strcpy(filename,argv[1]);
    //en argc van el nombre del ejecutable y el nombre de proyecto en gid va en argv
    vector<Matrix> localKs;
    vector<Vector> localbs;
    Matrix K;
    Vector b;
    Vector T;
    //leo los archivos para guardar todo en el objeto mesh 
    cout << "IMPLEMENTACI"<<char(224)<<"N DEL M"<<char(144)<<"TODO DE LOS ELEMENTOS FINITOS\n"
         << "\t- TRANSFERENCIA DE CALOR\n" << "\t- 2 DIMENSIONES\n"
         << "\t- FUNCIONES DE FORMA LINEALES\n" << "\t- PESOS DE GALERKIN\n"
         << "*********************************************************************************\n\n";

    mesh m; //clase
    leerMallayCondiciones(m,filename); //es decir leer el archiv con los nodos, cantidad de elementos etc
    cout << "Datos obtenidos correctamente\n********************\n";

    crearSistemasLocales(m,localKs,localbs); // mandamos el objeto mesh y dos vectores, uno para k y otro para B recordando que al final lo que queremos es:
    //kT=b
    showKs(localKs); showbs(localbs); // imprimimos lo que hay en k y b
    cout << "******************************\n";

    zeroes(K,m.getSize(NODES)); //formato a K
    zeroes(b,m.getSize(NODES)); //formato a B
    ensamblaje(m,localKs,localbs,K,b);
    showMatrix(K); showVector(b);
    cout << "******************************\n";
    //cout << K.size() << " - "<<K.at(0).size()<<"\n";
    //cout << b.size() <<"\n";

    applyNeumann(m,b);
    showMatrix(K); 
    showVector(b);
    cout << "******************************\n";
    //cout << K.size() << " - "<<K.at(0).size()<<"\n";
    //cout << b.size() <<"\n";

    applyDirichlet(m,K,b);
    showMatrix(K); showVector(b);
    cout << "******************************\n";
    //cout << K.size() << " - "<<K.at(0).size()<<"\n";
    //cout << b.size() <<"\n";

    zeroes(T,b.size());
    calculate(K,b,T);

    cout << "La respuesta es: \n";
    showVector(T);

    writeResults(m,T,filename);

    return 0;
}
