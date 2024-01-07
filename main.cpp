#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "string"
#include <limits>
using namespace std;

struct GlobalData {
    int simulationTime;
    int simulationStepTime;
    int conductivity;
    int alfa;
    int tot;
    int initialTemp;
    int density;
    int specificHeat;
    int nodesNumber;
    int elementsNumber;
};

struct DnDxDnDy {
    double** dnDx;
    double** dnDy;
    double* x;
    double* y;
    double* detJ;
    int* boundaryEdges;
    double* detJ1D;
};

struct shapeFactors {
    double** eta_array; //dn/deta
    double** ksi_array;// pochodne dn/dksi
    double** shape_func;  //N1-N4
    double*** shape_func_boundary;
    int nodes = 4;
    int points; //pc
    int edge_points;
    vector <double> w_array;
    vector <double> w_array1D;
};

struct shapeFactors calculateDEtaDKsi (int nodes, int n) {
    int points;
    double** eta_array; //dn/deta
    double** ksi_array; //dn/dksi
    vector <double> w_array;
    vector <double> w_array1D;

    //weight
    double** shape_func;
    double*** shape_func_boundary = new double** [4]; //dla 4 krawedzi

    vector <double> eta; // wspolrzedne pc 2D
    vector <double> ksi;
    vector <double> eta_boundary[4]; // wspolrzedne pc 1D
    vector <double> ksi_boundary[4];

    double a; //wspolczynniki
    double b;
    double w1; //wagi
    double w2;

    int edgePoints; //wspolrzedne punktow brzegowych

    switch (n) {
        case 1:
            points = 4;
            a = 1/sqrt(3);

            edgePoints = 2;

            eta.insert (eta.end(), {-a, -a, a, a});
            ksi.insert  (ksi.end(), {-a, a, -a, a});

            eta_boundary[0].insert(eta_boundary[0].end(), {-1, -1 });  //wspolrzedne dla dolnej krawedzi
            ksi_boundary[0].insert(ksi_boundary[0].end(), {-a, a});
            eta_boundary[1].insert(eta_boundary[1].end(), {-a, a });
            ksi_boundary[1].insert(ksi_boundary[1].end(), {1, 1});  //wspolrzedne dla prawej krawedzi
            eta_boundary[2].insert(eta_boundary[2].end(), {1, 1});
            ksi_boundary[2].insert(ksi_boundary[2].end(), {a, -a}); //gorna
            eta_boundary[3].insert(eta_boundary[3].end(), {-a, a});
            ksi_boundary[3].insert(ksi_boundary[3].end(), {-1, -1}); //lewa

            w_array.insert(w_array.end(), {1, 1, 1, 1}); //tablica wag

            w_array1D.insert(w_array1D.end(), {1, 1});
            break;
        case 2:
            points = 9;
            a = sqrt(3./5);
            edgePoints = 3;
            eta.insert (eta.end(), {-a, -a, -a, 0, 0, 0, a, a, a});
            ksi.insert  (ksi.end(),{-a, 0, a, -a, 0, a, -a, 0, a});

            eta_boundary[0].insert(eta_boundary[0].end(), {-1, -1, -1});  //wspolrzedne dla dolnej krawedzi
            ksi_boundary[0].insert(ksi_boundary[0].end(), {-a, 0, a});
            eta_boundary[1].insert(eta_boundary[1].end(), {-a, 0, a });
            ksi_boundary[1].insert(ksi_boundary[1].end(), {1, 1, 1});  //wspolrzedne dla prawej krawedzi
            eta_boundary[2].insert(eta_boundary[2].end(), {1, 1, 1});
            ksi_boundary[2].insert(ksi_boundary[2].end(), {a, 0, -a}); //gorna
            eta_boundary[3].insert(eta_boundary[3].end(), {-a, 0, a});
            ksi_boundary[3].insert(ksi_boundary[3].end(), {-1, -1, -1}); //lewa

            w1 = 5./9;
            w2 = 8./9;

            w_array.insert(w_array.end(),{w1 * w1, w2 * w1, w1 * w1, w1 * w2, w2 * w2, w1 * w2, w1 * w1, w2 * w1, w1 * w1}); //wagi zielone punkty

            w_array1D.insert(w_array1D.end(), {w1, w2, w1}); //wagi na krawedzi

            break;
        case 3:
            points = 16;
            a = sqrt(3./7 + 2./7 * sqrt(6./5));
            b = sqrt(3./7 - 2./7 * sqrt(6./5));
            edgePoints = 4;
            eta.insert (eta.end(),{-a, -a, -a, -a, -b, -b, -b,- b, b, b ,b ,b, a ,a, a, a });
            ksi.insert  (ksi.end(), {-a, -b, b, a, -a, -b, b, a, -a, -b, b, a, -a, -b, b, a});

            eta_boundary[0].insert(eta_boundary[0].end(), {-1, -1, -1, -1});  //wspolrzedne dla dolnej krawedzi hbc
            ksi_boundary[0].insert(ksi_boundary[0].end(), {-a, -b, b, a});
            eta_boundary[1].insert(eta_boundary[1].end(), {-a, -b, b, a});
            ksi_boundary[1].insert(ksi_boundary[1].end(), {1, 1, 1, 1});  //wspolrzedne dla prawej krawedzi
            eta_boundary[2].insert(eta_boundary[2].end(), {1, 1, 1, 1});
            ksi_boundary[2].insert(ksi_boundary[2].end(), {a, b, -b, -a}); //gorna
            eta_boundary[3].insert(eta_boundary[3].end(), {-a, -b, b, a});
            ksi_boundary[3].insert(ksi_boundary[3].end(), {-1, -1, -1, -1}); //lewa krawedz
            w1 = 0.347855;
            w2 = 0.652145;

            w_array.insert(w_array.end(), {w1 * w1, w1 * w2, w1 * w2, w1 * w1, w2 * w1, w2 * w2, w2 * w2, w2 * w1, w2 * w1, w2 * w2, w2 * w2, w2 * w1, w1 * w1, w1 * w2, w1 * w2, w1 * w1});

            w_array1D.insert(w_array1D.end(), {w1, w2, w2, w1});

            break;
    }

    eta_array = new double* [nodes];
    ksi_array = new double* [nodes];
    shape_func = new double* [nodes];
    for (int i = 0; i < nodes; ++i) {
        eta_array[i] = new double[points];// dn/deta
        ksi_array[i] = new double[points]; //dn / dksi
        shape_func[i] = new double[points]; //N
    }

    for (int j = 0; j < points; ++j) {
        shape_func [0][j] = 0.25 * (1 - eta[j])*( 1 - ksi[j] ); //N1
        eta_array[0][j] = -0.25 *  (1 - ksi[j]); // dN1/deta
        ksi_array[0][j] = -0.25 * (1 - eta[j]); // dN1/dksi
        shape_func [1][j] = 0.25 * (1 - eta[j])*( 1 + ksi[j] ); // N2
        eta_array[1][j] = -0.25 *  (1 + ksi[j]); //dN2/deta
        ksi_array[1][j] = 0.25 * (1 - eta[j]); //dN2/dksi
        shape_func [2][j] = 0.25 * (1 + eta[j])*( 1 + ksi[j] ); // N3
        eta_array[2][j] = 0.25 *  (1 + ksi[j]); //dN3/deta
        ksi_array[2][j] = 0.25 * (1 + eta[j]); //dN3/dksi
        shape_func[3][j] = 0.25 * (1 + eta[j])*( 1 - ksi[j] ); // N4
        eta_array[3][j] = 0.25 *  (1 - ksi[j]); //dN4/deta
        ksi_array[3][j] = -0.25 * (1 + eta[j]);//dN4/dksi
    }

    for (int i = 0; i < 4; ++i) {
        shape_func_boundary[i] = new double* [edgePoints];
        for (int j = 0; j < edgePoints; ++j) {
            shape_func_boundary[i][j] = new double [4];
            shape_func_boundary [i][j][0] = 0.25 * (1 - eta_boundary[i][j])*( 1 - ksi_boundary[i][j] );//N1 BC
            shape_func_boundary [i][j][1] = 0.25 * (1 - eta_boundary[i][j])*( 1 + ksi_boundary[i][j] );//N2
            shape_func_boundary [i][j][2] = 0.25 * (1 + eta_boundary[i][j])*( 1 + ksi_boundary[i][j] );//N3
            shape_func_boundary [i][j][3] = 0.25 * (1 + eta_boundary[i][j])*( 1 - ksi_boundary[i][j] );//N4
        }
    }

    for (int i = 0; i < points; ++i) {
        cout << "Punkt nr " << i+1 << endl;
        cout << "ksi = "<< ksi[i] << " eta = " << eta[i] << endl;
    }
    cout << endl;

    cout << "dn / deta: " << endl;
    for (int i = 0; i < points; ++i) {
        for (int j = 0; j < nodes; ++j) {
            cout <<  eta_array[j][i] << " ";
        }
        cout << endl;
    }

    cout << "dn / dksi: " << endl;
    for (int i = 0; i < points; ++i) {
        for (int j = 0; j < nodes; ++j) {
            cout <<  ksi_array[j][i] << " ";
        }
        cout << endl;
    }
    cout << endl;

    struct shapeFactors results; //zamkniecie tego w strukture elementu uniwersalnego
    results.eta_array = eta_array;
    results.ksi_array = ksi_array;
    results.shape_func = shape_func;
    results.nodes = nodes;
    results.points = points;
    results.w_array = w_array;

    results.edge_points = edgePoints;
    results.shape_func_boundary = shape_func_boundary;
    results.w_array1D = w_array1D;

    return results;
}

struct DnDxDnDy* calculateDnDxDnDy(double** points, int** elements, struct shapeFactors factors, struct GlobalData gd, int* BC){ //liczenie pochodnych jakobianow, sprawdzanie czy krawedz jest brzegowa, liczenie dn/dx dn/dy dla danego elementu
    struct DnDxDnDy* results = new struct DnDxDnDy [gd.elementsNumber]; //lista dndx/dndy tyle ile elementow

    for (int i = 0; i < gd.elementsNumber; ++i) {
        results[i].dnDx = new double* [factors.points];
        results[i].dnDy = new double* [factors.points];
        results[i].x = new double [factors.points]; //x dla punktu calkowania
        results[i].y = new double [factors.points];
        results[i].detJ = new double[factors.points];
        results[i].detJ1D = new double [4];
        results[i].boundaryEdges = new int [4]; //prawdziwe tylko dla czworokątnych elementów
        double* x = new double [factors.nodes];
        double* y = new double [factors.nodes];
        for (int j = 0; j < factors.nodes; ++j) { //przepisanie wspolrzednych wezlow w elemencie
            x[j] = points[elements[i][j]][0];
            y[j] = points[elements[i][j]][1];// j-ty punkt w i-tym elemencie
        }

        //wyznaczanie bc na danej krawedzi
        if( BC[elements[i][0]] == 1 && BC[elements[i][1]] == 1 ){
            results[i].boundaryEdges[0] = 1;
        }

        if( BC[elements[i][1]] == 1 && BC[elements[i][2]] == 1 ){
            results[i].boundaryEdges[1] = 1;
        }
        if( BC[elements[i][2]] == 1 && BC[elements[i][3]] == 1 ){
            results[i].boundaryEdges[2] = 1;
        }
        if( BC[elements[i][3]] == 1 && BC[elements[i][0]] == 1 ){
            results[i].boundaryEdges[3] = 1;                          //wyznacza ktore krawedzie sa BC
        }

        double J [2][2] ={{0,0},{0,0}};
        for (int j = 0; j < factors.points; ++j) { // j-ty punkt calkowania

            results[i].dnDx[j] = new double [factors.nodes];
            results[i].dnDy[j] = new double [factors.nodes];
            results[i].x[j] = 0;
            results[i].y[j] = 0;

            J[0][0]=0;
            J[0][1]=0;
            J[1][0]=0;
            J[1][1]=0;

            for (int k = 0; k < factors.nodes; ++k) {
                //znajdowanie wspolrzednych punktow calkowania znajac wspolrzedne punktow elementu
                results[i].x[j] += factors.shape_func[k][j] * x[k]; //x pc j
                results[i].y[j] += factors.shape_func[k][j] * y[k]; //y pc j
                //obliczanie elementów macierzy jakobiego
                J[0][0] += factors.ksi_array[k][j] * x[k]; //dx/dksi  = dN/dksi * x
                J[0][1] += factors.ksi_array[k][j] * y[k];
                J[1][0] += factors.eta_array[k][j] * x[k];
                J[1][1] += factors.eta_array[k][j] * y[k]; // dy/deta = dn/deta * y
            }
            double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0]; //detJ
            results[i].detJ[j] = detJ;
            for (int k = 0; k < factors.nodes; ++k) {
                results[i].dnDx[j][k] = 1./detJ * (J[1][1] * factors.ksi_array[k][j] - J[0][1] * factors.eta_array[k][j]);  //dN/dx
                results[i].dnDy[j][k] = 1./detJ * (-J[1][0] * factors.ksi_array[k][j] + J[0][0] * factors.eta_array[k][j]); //dN/dy
            }
            cout<<endl;
            cout << "wartosc dN/dx rowna:: " << endl;
            for (int k = 0; k < factors.nodes; ++k) {
                cout << results[i].dnDx[j][k] << ", ";
            }
            cout << endl;
            cout << "wartosc dN/dy rowna: " << endl;
            for (int k = 0; k < factors.nodes; ++k) {
                cout << results[i].dnDy[j][k] << " ";
            }
            cout<< endl;

            cout << "Macierz Jakobiego dla elementu nr "<< i+1 << " i punktu "<< j+1 << endl;
            cout << J[0][0] << ",  " << J[0][1] << endl;
            cout << J[1][0] << ",  " << J[1][1] << endl;
            cout << "DetJ = " << detJ << endl;

        }
        //det[j] = dx/dksi = L/2
        results[i].detJ1D[0] = 0.5 * sqrt( (x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]) );
        results[i].detJ1D[1] = 0.5 * sqrt( (x[2] - x[1]) * (x[2] - x[1]) + (y[2] - y[1]) * (y[2] - y[1]) );
        results[i].detJ1D[2] = 0.5 * sqrt( (x[2] - x[3]) * (x[2] - x[3]) + (y[2] - y[3]) * (y[2] - y[3]) );
        results[i].detJ1D[3] = 0.5 * sqrt( (x[3] - x[0]) * (x[3] - x[0]) + (y[3] - y[0]) * (y[3] - y[0]) );
    }
    return results;
}

double*** matrixH (struct GlobalData globalData, struct shapeFactors factors,struct DnDxDnDy* dnDxDnDy_array){
    double*** matrix = new double**[globalData.elementsNumber];
    for (int i = 0; i < globalData.elementsNumber; ++i) {
        matrix[i] = new double*[factors.nodes];
        for (int j = 0; j < factors.nodes; ++j) { // alokacja

            matrix[i][j] = new double[factors.nodes];
            for (int k = 0; k < factors.nodes; ++k) {
                matrix[i][j][k] = 0;
            }
        }
        for (int j = 0; j < factors.points; ++j) {

            for (int k = 0; k < factors.nodes; ++k) {
                for (int l = 0; l < factors.nodes; ++l) {
                    matrix[i][k][l] += (dnDxDnDy_array[i].dnDx[j][k] * dnDxDnDy_array[i].dnDx[j][l] + dnDxDnDy_array[i].dnDy[j][k] * dnDxDnDy_array[i].dnDy[j][l]) * dnDxDnDy_array[i].detJ[j] * factors.w_array[j] * globalData.conductivity;
                }
            }
        }
    }
    return matrix;
}

double*** matrixHbc (struct GlobalData globalData, struct shapeFactors factors, struct DnDxDnDy* dnDxDnDy_array) {
    double*** matrix = new double**[globalData.elementsNumber];
    for (int i = 0; i < globalData.elementsNumber; ++i) {
        matrix[i] = new double*[factors.nodes];
        for (int j = 0; j < factors.nodes; ++j) { // alokacja
            matrix[i][j] = new double[factors.nodes];
            for (int k = 0; k < factors.nodes; ++k) {
                matrix[i][j][k] = 0;  //zerowanie
            }
        }
        //sprawdzenie czy punkty sa brzegowe
        for (int j = 0; j < 4; ++j) { //deklaruje krawedzie ktore sa
            if (dnDxDnDy_array[i].boundaryEdges[j] == 1){ //czy dana krawedz zawiera warunek brzegowy
                for (int m = 0; m < factors.edge_points; ++m) {
                    for (int k = 0; k < factors.nodes; ++k) { //numeracja funkcji ksztaltu
                        for (int l = 0; l < factors.nodes; ++l) { //numeracja funkcji ksztaltu
                            matrix[i][k][l] += factors.shape_func_boundary[j][m][k] * factors.shape_func_boundary[j][m][l] * factors.w_array1D[m] * globalData.alfa * dnDxDnDy_array[i].detJ1D[j]; //Hbc
                        }
                    }
                }
            }
        }
    }
    return matrix;
}

double*** C (struct GlobalData globalData, struct shapeFactors factors, struct DnDxDnDy* dnDxDnDy_array){ //macierz C
    double*** matrix = new double**[globalData.elementsNumber]; //alokacja
    for (int i = 0; i < globalData.elementsNumber; ++i) {
        matrix[i] = new double*[factors.nodes];//alokowane
        for (int j = 0; j < factors.nodes; ++j) {
            matrix[i][j] = new double[factors.nodes]; //alokowane
            for (int k = 0; k < factors.nodes; ++k) {
                matrix[i][j][k] = 0;
            }
        }
        for (int j = 0; j < factors.points; ++j) {
            for (int k = 0; k < factors.nodes; ++k) {
                for (int l = 0; l < factors.nodes; ++l) {
                    matrix[i][k][l] += factors.shape_func[k][j] * factors.shape_func[l][j] * globalData.density * globalData.specificHeat * dnDxDnDy_array[i].detJ[j] * factors.w_array[j];
                }
            }
        }
    }
    return matrix;
}

double*** HAndHbc (struct GlobalData globalData , double*** H, double***Hbc, struct shapeFactors factors){  //sumowanie H+Hbc
    double*** fullH = new double ** [globalData.elementsNumber]; //alokacja

    for (int i = 0; i < globalData.elementsNumber; ++i) {
        fullH[i] = new double* [factors.nodes]; // jak to wytlumaczyc?
        for (int j = 0; j < factors.nodes; ++j) {
            fullH[i][j] = new double [factors.nodes];//
            for (int k = 0; k < factors.nodes; ++k) {
                fullH[i][j][k] = H[i][j][k] + Hbc[i][j][k];
            }
        }
    }
    return fullH;
}

double** P (struct GlobalData globalData, struct shapeFactors factors, struct DnDxDnDy* dnDxDnDy){ // obliczanie wektora P
    double**P = new double * [globalData.elementsNumber];
    for (int i = 0; i < globalData.elementsNumber; ++i) {
        P[i] = new double [4];
        for (int j = 0; j < 4; ++j) {
            P[i][j] = 0;
        }
        for (int j = 0; j < 4; ++j) {
            if(dnDxDnDy[i].boundaryEdges[j] == 1){
                for (int k = 0; k < factors.edge_points; ++k) {
                    for (int l = 0; l < factors.nodes; ++l) {
                        P[i][l] += factors.shape_func_boundary[j][k][l] * factors.w_array1D[k] * globalData.alfa * globalData.tot * dnDxDnDy[i].detJ1D[j];
                    }
                }
            }
        }
    }
    return P;
}

double* agregateP (double** P, struct GlobalData globalData, int** elements) { //agregacja wektora P
    double* agregateVector = new double[globalData.nodesNumber];
    for (int i = 0; i < globalData.nodesNumber; ++i) {
        agregateVector[i] = 0;
    }
    for (int i = 0; i < globalData.elementsNumber; ++i) {
        for (int j = 0; j < 4; ++j) {
            agregateVector[elements[i][j]] += P[i][j];
        }
    }
    return agregateVector;
}

double** agregation (double*** H, struct GlobalData globalData, int** elements){  //agregacja dla  FullH i C
    double** agregatedMatrix = new double* [globalData.nodesNumber];
    for (int i = 0; i < globalData.nodesNumber; ++i) {
        agregatedMatrix[i] = new double [globalData.nodesNumber];
        for (int j = 0; j < globalData.nodesNumber; ++j) {
            agregatedMatrix[i][j] = 0;
        }
    }
    for (int i = 0; i < globalData.elementsNumber; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                agregatedMatrix[elements[i][j]][elements[i][k]] += H[i][j][k];
            }
        }
    }
    return agregatedMatrix;
}

double* gauss (double** Hinitial, double* Pinitial, int N) { //N - ilosc wezlow
    double** H = new double* [N];
    double* P = new double[N];
    for (int i = 0; i < N; ++i) {
        H[i] = new double [N];
        for (int j = 0; j <  N; ++j) {
            H[i][j] = Hinitial[i][j]; //przypisanie wyniku z HandC sumMatrix
        }
        P[i] = Pinitial[i]; //PandC
    }
    for (int i = 0; i < N; ++i) {
        double c = H[i][i]; //wspolczynnik w algorytmie gaussa
        for (int j = i; j < N; ++j) {
            H[i][j] = H[i][j] / c; //uzyskanie 1 na diagonali
        }
        P[i]=P[i] / c; //dzielenie c dla prawej strony rownania

        for (int j = i+1; j < N; ++j) {
            c = H[j][i];
            for (int k = i; k < N; ++k) {
                H[j][k] -= c * H[i][k];
            }
            P[j] -= P[i] * c;
        }
    }
    for (int i = N-1; i >= 0 ; --i) {
        for (int j = i-1; j >= 0; --j){
            double c = H[j][i];
            for (int k = i; k >= j; --k) {
                H[j][k] -= c * H[i][k];
            }
            P[j] -= P[i] * c;
        }
    }
    return P;
}
double** sumMatrix (double** H, double** C, struct GlobalData globalData){ //HandC
    double** matrix = new double* [globalData.nodesNumber];
    for (int i = 0; i < globalData.nodesNumber; ++i) {
        matrix[i] = new double [globalData.nodesNumber];
        for (int j = 0; j < globalData.nodesNumber; ++j) {
            matrix[i][j] = H[i][j] + C[i][j]/globalData.simulationStepTime; // H + C/delta Tau
        }
    }
    return matrix;
}

double* PandC (double** C, double* P, double* temp, struct GlobalData globalData){ //sumowanie wektora P i macierzy C, druga czesc rownania niestac ostatecznego
    double* vec = new double [globalData.nodesNumber];
    for (int i = 0; i < globalData.nodesNumber; ++i) {
        vec[i] = 0;
        for (int j = 0; j < globalData.nodesNumber; ++j) {
            vec[i] += C[i][j] * temp[j]; //C * t0
        }
        vec[i] = P[i] + vec[i] / globalData.simulationStepTime; // (C*t0)/delta Tau + P
    }
    return vec;
}

double** symulacjaCzasowa (double** H, double** C, double* P, struct GlobalData globalData){
    int t = 0;
    int N = globalData.simulationTime / globalData.simulationStepTime + 1;

    double** temp = new double* [N]; // przechowywanie temperatury w i-tym kroku czasowym w j-tym węźle
    for (int i = 0; i < N; ++i) {
        temp[i] = new double [globalData.nodesNumber];
        for (int j = 0; j < globalData.nodesNumber; ++j) {
            temp[i][j] = 0;
        }
    }
    double** HandC = sumMatrix(H, C, globalData);  //H+C lewa czesc rownania
    for (int i = 0; i < globalData.nodesNumber; ++i) {
        temp[0][i] = globalData.initialTemp; //temp w 0 kroku czasowym t0
    }

    for (int i = 0; i < N-1; ++i) {
        double* vec = PandC(C, P, temp[i], globalData); //prawa czesc rownania + P

        temp[i+1] =  gauss(HandC, vec, globalData.nodesNumber); //t1
    }
    return temp;
}

void writeToVTK (double** temp, struct GlobalData globalData, int** elements, double** nodes, string filename){
    int N = globalData.simulationTime / globalData.simulationStepTime + 1; //ilośc krokow czasowych
    for (int i = 0; i < N; ++i) {
        string name = filename + to_string(i) + ".vtk";
        ofstream file;
        file.open(name);
        file << "# vtk DataFile Version 2.0" << endl;
        file << "Unstructured Grid Example"<< endl;
        file << "ASCII"<< endl;
        file << "DATASET UNSTRUCTURED_GRID"<< endl;
        file << endl;
        file << "POINTS " << globalData.nodesNumber << " float" << endl;
        for (int j = 0; j < globalData.nodesNumber; ++j) {
            file << nodes[j][0] << " " << nodes[j][1] << " 0" << endl;
        }
        file << endl;
        file << "CELLS " << globalData.elementsNumber << " " << globalData.elementsNumber * 5 << endl;
        for (int j = 0; j < globalData.elementsNumber; ++j) {
            file << "4 " << elements[j][0] << " " << elements[j][1] << " " << elements[j][2] << " " << elements[j][3] << endl;
        }
        file << endl;
        file << "CELL_TYPES " << globalData.elementsNumber << endl;
        for (int j = 0; j < globalData.elementsNumber; ++j) {
            file << "9" << endl;
        }
        file << endl;
        file << "POINT_DATA "  << globalData.nodesNumber << endl;
        file << "SCALARS Temp float 1" << endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int j = 0; j < globalData.nodesNumber; ++j) {
            file << temp[i][j] << endl;
        }
        file.close();
    }
}

int main() {
    int** elements;
    double** nodes;
    int*BC;

    struct GlobalData globalData;
    string fileName;
    cout << "Podaj plik siatki: \n";
    cin >> fileName;
    ifstream file;
    file.open(fileName.c_str());

    if (file.good()){
        int paramValue;
        string paramName, paramName1;
        for (int i = 0; i < 10; ++i) {
            file >> paramName >> paramName1;

            if (paramName1 == "number"){
                file >> paramValue;
            }
            else {
                paramValue = stoi(paramName1);
            }
            if (paramName == "SimulationTime"){
                globalData.simulationTime = paramValue;
            }
            else if (paramName == "SimulationStepTime"){
                globalData.simulationStepTime = paramValue;
            }
            else if (paramName =="Conductivity"){
                globalData.conductivity = paramValue;
            }
            else if (paramName =="Alfa"){
                globalData.alfa = paramValue;
            }
            else if (paramName =="Tot"){
                globalData.tot = paramValue;
            }
            else if (paramName =="InitialTemp"){
                globalData.initialTemp = paramValue;
            }
            else if (paramName =="Density"){
                globalData.density = paramValue;
            }
            else if (paramName =="SpecificHeat"){
                globalData.specificHeat = paramValue;
            }
            else if (paramName =="Nodes"){
                globalData.nodesNumber = paramValue;
            }
            else if (paramName =="Elements"){
                globalData.elementsNumber = paramValue;

            }
            else{
                cout << "nieznany parametr: " << paramName << endl;
            }
        }

        nodes = new double* [globalData.nodesNumber];
        for (int i = 0; i < globalData.nodesNumber; ++i) {
            nodes[i] = new double[2];
        }
        elements = new int* [globalData.elementsNumber];
        for (int i = 0; i < globalData.elementsNumber; ++i) {
            elements[i] = new int[4];
        }

        string xx;
        file >> xx;
        for (int i = 0; i < globalData.nodesNumber; ++i) {
            int n;
            float x, y;
            char comma;
            file >> n >> comma >> x >> comma >> y;
            nodes [i][0] = x;
            nodes [i][1] = y;
        }
        file >> xx;
        file >> xx;
        for (int i = 0; i < globalData.elementsNumber; ++i) {
            char comma;
            int n, a, b, c, d;
            file >> n >> comma >> a >> comma >> b >> comma >> c >> comma >>d;
            elements [i][0] = a-1;
            elements [i][1] = b-1;
            elements [i][2] = c-1;
            elements [i][3] = d-1;
        }

        BC = new int[globalData.nodesNumber]; //wczytywanie wezlow brzegowych

        for (int i = 0; i < globalData.nodesNumber; ++i) {
            BC[i] = 0;
        }

        file >> xx;

        vector<int> boundaryNodes; //zapisujemy BC w wektorze

        while(file >> xx){
            int k = xx.find(",");
            boundaryNodes.push_back(stoi(xx.substr(0, k))-1); // od 0 do indexu po ktorym znalazlo ", ", stoi - konwersja na int, push_back dodaje do vectora
        }
        for (auto it = boundaryNodes.begin(); it!= boundaryNodes.end(); it++){ //petla do iterowania po vectorze BC
            BC[*it] = 1;
        }

    }

    else {
        cout<< "Bląd odczytu siatki" << endl;
        return 1;
    }

    cout << "SimulationTime: " << globalData.simulationTime << endl;
    cout << "SimulationStepTime: " << globalData.simulationStepTime << endl;
    cout << "Conductivity: " << globalData.conductivity << endl;
    cout << "Alfa: " << globalData.alfa << endl;
    cout << "Tot: " << globalData.tot << endl;
    cout << "InitialTemp: " << globalData.initialTemp << endl;
    cout << "Density: " << globalData.density << endl;
    cout << "SpecificHeat: " << globalData.specificHeat << endl;
    cout << "Nodes Number:" << globalData.nodesNumber << endl;
    cout << "Elements Number: " << globalData.elementsNumber << endl;

    cout << endl;

    for (int i = 0; i < globalData.nodesNumber; ++i) {
        cout << " " << "x: " << nodes[i][0] << " " << "y: " << nodes[i][1] << " status = " << BC[i] <<endl;
    }

    cout << endl;

    for (int i = 0; i < globalData.elementsNumber; ++i) {
        cout << "ID"<< i+1 << "[" << elements[i][0] + 1 << "," << elements[i][1] + 1 << "," << elements[i][2] + 1  << "," << elements[i][3] + 1 << "]"<< endl;
    }

    cout << endl;

    cout << "Wybor schematu calkowania calkowania (2,3,4) w formacie n+1, wprowadz n, (przyklad schemat 2 pkt n = 1):" << endl;
    int c = 0;
    cin >> c;

    switch (c) {
        case 1: {
            cout << "Kolejnosc punktow calkowania" << endl;
            struct shapeFactors factors2 = calculateDEtaDKsi(4, 1);
            struct DnDxDnDy* dxdy2 = calculateDnDxDnDy(nodes, elements, factors2, globalData, BC);
            double*** H = matrixH(globalData, factors2,dxdy2);
            double*** Hbc = matrixHbc(globalData, factors2, dxdy2);
            double*** fullH = HAndHbc(globalData, H, Hbc, factors2);
            double** vecP = P(globalData, factors2, dxdy2);
            double*** matrixC = C(globalData, factors2, dxdy2);
            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Macierz H dla elementu " << i+1 << endl;
                for (int j = 0; j < factors2.nodes; ++j) {
                    for (int k = 0; k < factors2.nodes; ++k) {
                        cout << H[i][j][k] << "  ";
                    }
                    cout << endl;
                }
            }

            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Macierz Hbc dla elementu " << i+1 << endl;
                for (int j = 0; j < factors2.nodes; ++j) {
                    for (int k = 0; k < factors2.nodes; ++k) {
                        cout << Hbc[i][j][k] << "  ";
                    }
                    cout << endl;
                }
            }
            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Macierz H+Hbc dla elementu " << i+1 << endl;
                for (int j = 0; j < factors2.nodes; ++j) {
                    for (int k = 0; k < factors2.nodes; ++k) {
                        cout << fullH[i][j][k] << "  ";
                    }
                    cout << endl;
                }
            }
            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Wektor P dla elementu " << i+1 << endl;
                for (int j = 0; j < factors2.nodes; ++j) {

                    cout << vecP[i][j] << "  ";
                }
                cout<<endl;
            }
            double** agregatedMatrix = agregation(fullH, globalData, elements); //agregacja FullH
            double** agregatedMatrixC = agregation(matrixC, globalData, elements); //zagregowana macierz C

            double* agregateVector = agregateP(vecP, globalData, elements);//agregacja wektora P
            cout << "Agregacja wektora P" << endl;
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                cout << agregateVector[i] << " ";
            }
            cout<<endl;
            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout << "Macierz C dla elementu - " << i+1 << endl;
                for (int j = 0; j < factors2.nodes; ++j) {
                    for (int k = 0; k < factors2.nodes; ++k) {
                        cout << matrixC[i][j][k] << "  ";
                    }
                    cout << endl;

                }
            }
            cout<<endl;
            double** HandC = sumMatrix(agregatedMatrix, agregatedMatrixC, globalData);
            cout<<"H+C" << endl;
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                for (int j = 0; j < globalData.nodesNumber; ++j) {
                    cout<< HandC[i][j] << " ";
                }
                cout << endl;
            }
            double* solution = gauss (agregatedMatrix, agregateVector, globalData.nodesNumber); //rozwiazanie stacjonarne
            cout<<"Temp ukl. stacj."<<endl;
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                cout << solution[i] << " ";
            }
            cout<<endl;
            cout<<"Agregacja Macierzy FullH"<<endl;
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                for (int j = 0; j < globalData.nodesNumber; ++j) {
                    cout<< agregatedMatrix[i][j] << " ";
                }
                cout << endl;
            }
            cout<< endl;
            //symulacja
            double** temp = symulacjaCzasowa(agregatedMatrix, agregatedMatrixC, agregateVector, globalData);

            cout<<"Symulacja: "<< endl;

            for (int i = 0; i <  globalData.simulationTime / globalData.simulationStepTime + 1; ++i) {
                double minT = numeric_limits<double>::max();
                double maxT = numeric_limits<double>::min();
                cout<<"T = " << i * globalData.simulationStepTime << endl;
                for (int j = 0; j < globalData.nodesNumber; ++j) {
                    //cout<< temp[i][j] << " ";
                    if(temp[i][j] > maxT) {
                        maxT = temp[i][j];
                    }
                    if(temp [i][j] < minT) { //i krok czasowy j element
                        minT = temp[i][j];
                    }

                }
                //cout<< endl;
                cout << "T min= " << minT << " " << "T max= " << maxT << endl;
            }
            //writeToVTK(temp, globalData, elements, nodes, "Test_4_");
            break;
        }
        case 2: {
            struct shapeFactors factors3 = calculateDEtaDKsi(4, 2);
            struct DnDxDnDy* dxdy3 = calculateDnDxDnDy(nodes, elements, factors3, globalData, BC);
            double*** H = matrixH(globalData, factors3,dxdy3);
            double*** Hbc = matrixHbc(globalData, factors3, dxdy3);
            double*** fullH = HAndHbc(globalData, H, Hbc, factors3);
            double** vecP = P(globalData, factors3, dxdy3);
            double*** matrixC = C(globalData, factors3, dxdy3);

            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Macierz H dla elementu " << i+1 << endl;
                for (int j = 0; j < factors3.nodes; ++j) {
                    for (int k = 0; k < factors3.nodes; ++k) {
                        cout << H[i][j][k] << " ";
                    }
                    cout << endl;
                }
            }
            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Macierz Hbc dla elementu " << i+1 << endl;
                for (int j = 0; j < factors3.nodes; ++j) {
                    for (int k = 0; k < factors3.nodes; ++k) {
                        cout << Hbc[i][j][k] << "  ";
                    }
                    cout << endl;
                }
            }
            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Macierz H+Hbc dla elementu " << i+1 << endl;
                for (int j = 0; j < factors3.nodes; ++j) {
                    for (int k = 0; k < factors3.nodes; ++k) {
                        cout << fullH[i][j][k] << "  ";
                    }
                    cout << endl;
                }
            }
            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Wektor P dla elementu " << i+1 << endl;
                for (int j = 0; j < factors3.nodes; ++j) {

                    cout << vecP[i][j] << "  ";
                }
                cout<<endl;
            }
            double** agregatedMatrix = agregation(fullH, globalData, elements);
            double** agregatedMatrixC = agregation(matrixC, globalData, elements);
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                for (int j = 0; j < globalData.nodesNumber; ++j) {
                    cout<< agregatedMatrix[i][j] << " ";
                }
                cout << endl;
            }
            cout<<endl;
            double* agregateVector = agregateP(vecP, globalData, elements);
            cout << "Agregacja wektora P" << endl;
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                cout << agregateVector[i] << " ";
            }
            cout<<endl;

            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout << "Macierz C dla elementu - " << i+1 << endl;
                for (int j = 0; j < factors3.nodes; ++j) {
                    for (int k = 0; k < factors3.nodes; ++k) {
                        cout << matrixC[i][j][k] << "  ";
                    }
                    cout << endl;

                }
            }
            cout<<endl;
            double** HandC = sumMatrix(agregatedMatrix, agregatedMatrixC, globalData);
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                for (int j = 0; j < globalData.nodesNumber; ++j) {
                    cout<< HandC[i][j] << " ";
                }
                cout << endl;
            }

            cout << "Gauss elim" << endl;
            double* solution =  gauss (agregatedMatrix, agregateVector, globalData.nodesNumber);
            cout<<"Temperatury ukl. stacj."<<endl;
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                cout << solution[i] << " ";
            }
            cout<<endl;
            cout<<"Agregacja Macierzy fullH"<<endl;
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                for (int j = 0; j < globalData.nodesNumber; ++j) {
                    cout<< agregatedMatrix[i][j] << " ";
                }
                cout << endl;
            }
            cout<<endl;
            double** temp = symulacjaCzasowa(agregatedMatrix, agregatedMatrixC, agregateVector, globalData); // rozklad temperatur w danej chwili czasu

            cout<<"Symulacja: "<< endl;

            for (int i = 0; i <  globalData.simulationTime / globalData.simulationStepTime + 1; ++i) {
                double minT = numeric_limits<double>::max();
                double maxT = numeric_limits<double>::min();
                cout<<"T = " << i * globalData.simulationStepTime << endl;
                for (int j = 0; j < globalData.nodesNumber; ++j) {
                    //cout<< temp[i][j] << " ";
                    if(temp[i][j] > maxT) {
                        maxT = temp[i][j];
                    }
                    if(temp [i][j] < minT) {
                        minT = temp[i][j];
                    }
                }
                cout << "T min= " << minT << " " << "T max= " << maxT << endl;
            }
            //  writeToVTK(temp, globalData, elements, nodes, "Test_2_");
            break;

        }

        case 3: {
            cout << "Kolejnosc punktow calkowania" << endl;
            struct shapeFactors factors4 = calculateDEtaDKsi(4, 3);
            struct DnDxDnDy* dxdy4 = calculateDnDxDnDy(nodes, elements, factors4, globalData, BC);
            double*** H = matrixH(globalData, factors4,dxdy4);
            double*** Hbc = matrixHbc(globalData, factors4, dxdy4);
            double*** fullH = HAndHbc(globalData, H, Hbc, factors4);
            double** vecP = P(globalData, factors4, dxdy4);
            double*** matrixC = C(globalData, factors4, dxdy4);

            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Macierz H dla elementu " << i+1 << endl;
                for (int j = 0; j < factors4.nodes; ++j) {
                    for (int k = 0; k < factors4.nodes; ++k) {
                        cout << H[i][j][k] << " ";
                    }
                    cout << endl;
                }
            }

            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Macierz Hbc dla elementu " << i+1 << endl;
                for (int j = 0; j < factors4.nodes; ++j) {
                    for (int k = 0; k < factors4.nodes; ++k) {
                        cout << Hbc[i][j][k] << "  ";
                    }
                    cout << endl;
                }
            }
            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Macierz H+Hbc dla elementu " << i+1 << endl;
                for (int j = 0; j < factors4.nodes; ++j) {
                    for (int k = 0; k < factors4.nodes; ++k) {
                        cout << fullH[i][j][k] << "  ";
                    }
                    cout << endl;
                }
            }
            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout<< "Wektor P dla elementu " << i+1 << endl;
                for (int j = 0; j < factors4.nodes; ++j) {

                    cout << vecP[i][j] << "  ";
                }
                cout<<endl;
            }
            double** agregatedMatrix = agregation(fullH, globalData, elements);
            double** agregatedMatrixC = agregation(matrixC, globalData, elements);

            double* agregateVector = agregateP(vecP, globalData, elements);
            cout << "Agregacja wektora P" << endl;
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                cout << agregateVector[i] << " ";
            }
            cout<<endl;
            for (int i = 0; i < globalData.elementsNumber; ++i) {
                cout << "Macierz C dla elementu - " << i+1 << endl;
                for (int j = 0; j < factors4.nodes; ++j) {
                    for (int k = 0; k < factors4.nodes; ++k) {
                        cout << matrixC[i][j][k] << "  ";
                    }
                    cout << endl;

                }
            }
            cout<<endl;
            double** HandC = sumMatrix(agregatedMatrix, agregatedMatrixC, globalData);
            cout<<"H+C"<<endl;
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                for (int j = 0; j < globalData.nodesNumber; ++j) {
                    cout<< HandC[i][j] << " ";
                }
                cout << endl;
            }
            cout << "Gauss elim" << endl;
            double* wynik = gauss (agregatedMatrix, agregateVector, globalData.nodesNumber);
            cout<<"Temp ukl. stacj."<<endl;
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                //cout << agregateVector[i] << " ";
                cout<<wynik[i] << " ";
            }
            cout<<endl;

            cout<<"Agregacja macierzy fullH"<<endl;
            for (int i = 0; i < globalData.nodesNumber; ++i) {
                for (int j = 0; j < globalData.nodesNumber; ++j) {
                    cout<< agregatedMatrix[i][j] << " ";
                }
                cout << endl;
            }//
            cout<<endl;
            double** temp = symulacjaCzasowa(agregatedMatrix, agregatedMatrixC, agregateVector, globalData); // rozklad temperatur w danej chwili czasu

            cout<<"Symulacja: "<< endl;

            for (int i = 0; i <  globalData.simulationTime / globalData.simulationStepTime + 1; ++i) {
                double minT = numeric_limits<double>::max();
                double maxT = numeric_limits<double>::min();
                cout<<"T = " << i * globalData.simulationStepTime << endl;
                for (int j = 0; j < globalData.nodesNumber; ++j) {
                    //cout<< temp[i][j] << " ";
                    if(temp[i][j] > maxT) {
                        maxT = temp[i][j];
                    }
                    if(temp [i][j] < minT) {
                        minT = temp[i][j];
                    }

                }
                //cout<< endl;
                cout << "temp min= " << minT << " " << "temp max= " << maxT << endl;
            }
            // writeToVTK(temp, globalData, elements, nodes, "Test_4_");
            delete[] H;
            delete[] Hbc;
            delete[] fullH;
            delete[] vecP;
            delete[] matrixC;
            delete[] agregatedMatrix;
            delete[] agregatedMatrixC;
            delete[] agregateVector;
            break;
        }
        default:
            cout<< "Wprowadzono niepoprawna wartosc ilosci punktow calkowania"<< endl;
            break;
    }
    delete[] nodes;
    delete[] elements;
    delete[] BC;

    return 0;
}
