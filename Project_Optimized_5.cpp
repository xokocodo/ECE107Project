#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <sstream>

//Distances in mm
#define SQRTN 20
#define N SQRTN*SQRTN
#define PI 3.1415926
#define Eo 8.854*pow(10,-12)
#define Vo 50


using namespace std;

float absval(float x){
if(x<=0){
    return -x;
}
else
    return x;
}

void MultiplyMatrices(vector<float>& I,float V[N*2],float Q[N*2]){
int i,j=0;

for(i=0;i<N*2;i++){
for(j=0;j<N*2;j++){
Q[i]+= I[2*N*i+j]*V[j];
}
}

return;
}

void inverse_gaussian(vector<float>& M,vector<float>& I){
   int k,l,i,j;
   float factor;
   float coefficient;

   for(k=0; k<N*2; k++){
       // cout << "Pivot = " << k << endl;
      //one in the pivot
      factor = M[2*N*k+k];
      for(i=0;i<N*2;i++){
        M[2*N*k+i] = M[2*N*k+i]/factor;
        I[2*N*k+i] = I[2*N*k+i]/factor;
            }
      //zeroing the column
      for(l=0; l<N*2; l++){
        coefficient = M[2*N*l+k];
        for(i=0;i<N*2;i++){
            if(k!=l){
            M[2*N*l+i]-=coefficient*M[2*N*k+i];
            I[2*N*l+i]-=coefficient*I[2*N*k+i];
            if(I[2*N*l+i]==0)
                I[2*N*l+i]=0;
            if(M[2*N*l+i]==0)
                M[2*N*l+i]=0;
            }
        }
      }
   }
   return;
}

float determinant(vector<float>& M, int x, int y, int Nt){
int i,j;
float det=0;

for(i=0; i<=Nt; i++){
det+=M[N*y+(x+i)]*determinant(M, x+1, y+1, Nt-1);
}


}

void inverse_analytical(vector<float>& M,vector<float>& I, int Nt){

   return;
}

float absdist(int i,int j,float k,int l,int m,float n){
float dist=0;
dist = sqrt(pow((i-l),2)+pow((j-m),2)+pow((k-n),2));
//cout  << i << " " << j << " " << k << "::" << l << " " << m << " " << n << "::" << dist << endl;
return dist;
}

void Identity(vector<float>& I){
int l,k;
for(k=0;k<N*2;k++){  //Create Identity Vector
for(l=0;l<N*2;l++){
    if(k==l){
        I[2*N*l+k]=1;
        }
    else{
        I[2*N*l+k]=0;
        }
   }
   }

}


int main(){

int i,j,k,l,m,n;  // i=x' j=y' k=z'   l=x m=y n=z
int counter1=0;   //Counter for the NxN Z array
int counter2=0;   //Counter for the NxN Z array

vector<float> Z(4*N*N);   // NxN Z Vector
float V[N*2]; // Nx1 Voltage (V) Vector
float Q[N*2]; // Nx1 Charge (Q) Vector
vector<float> I(4*N*N);  // Nx1 Identity Vector
float D;
time_t timer;
time(&timer);
string filename;
ofstream output, log;
std::ostringstream s;

log.open("log.txt");

float TotalQ;         // Sum of Q
float surface1[SQRTN][SQRTN];
float surface2[SQRTN][SQRTN];

for(int d=1; d<=40; d++)  //d is in mm
{
counter1=0;
counter2=0;

D = d*(((float)SQRTN)/(10));
cout << "D(in units)=" << D << endl;
cout << "Each unit is X mm:" << 10/((float)SQRTN) << endl;
cout << "Each mm is X units:" << ((float)SQRTN)/10 << endl;

s.str(std::string());
s <<  "Output_D_" << D << "_N_" << N << "__" << timer << ".txt";
filename=s.str();
cout << filename <<endl;
output.open(filename.c_str());  //Open Output File

Identity(I);  // Set Identity Vector

    for(n=0;n<=1;n++){
        for(l=0;l<SQRTN;l++){
            for(m=0;m<SQRTN;m++){  // rm = r
                for(k=0;k<=1;k++){
                    for(i=0;i<SQRTN;i++){   //rn = r'
                        for(j=0;j<SQRTN;j++){
                            if((i!=l)||(j!=m)||(k!=n)){
                                Z[counter1*2*N + counter2]= 1/(4*PI*Eo*absdist(i+.5,j+.5,k*D,l+.5,m+.5,n*D));  // rm=/=rn
                            }
                            else{
                                Z[counter1*2*N + counter2]= 1/(2*Eo*sqrt(PI));    // rm==rn
                            }
                            counter2++;
                            }
                        }
                    }
                if(n==1){         //If top plate, set V=Vo/2
                V[counter1]=Vo/2;
                }
                else          //If bottom plate, set V=-Vo/2
                {
                V[counter1]=-Vo/2;
                }
                Q[counter1]=0;
                counter1++;
                counter2=0;
            }
        }
    }


cout << "Initial Z Complete" << endl;

inverse_gaussian(Z,I);  // Invert Z

cout << "Inversion Complete" << endl;

MultiplyMatrices(I,V,Q);   //  Q = Z^-1 * V

cout << "Multiplication Complete" << endl;

TotalQ=0.0;         // Sum of Q

counter1=0;
counter2=0;

output << "Surface 1" << endl;

for(l=0;l<SQRTN;l++){
for(m=0;m<SQRTN;m++){
surface1[l][m]=Q[counter1];
output /*<< l <<" " << m << ":: "*/ << surface1[l][m] << ",";
counter1++;
}
output << endl;
}

output << endl << endl << endl << endl;

output << "Surface 2" << endl;
for(l=0;l<SQRTN;l++){
for(m=0;m<SQRTN;m++){
surface2[l][m]=Q[counter2+N];
output /*<< l <<" " << m << ":: "*/<< surface2[l][m] << ",";
counter2++;
}
output << endl;
}

output << endl << endl << endl << endl;


output << "Charge Layout" << endl;
for(i=0;i<N*2;i++){        //Sum Total Q
output << Q[i] << endl;
TotalQ+=absval(Q[i]);      //absolute value - to make positive
}

output << endl << endl << endl << endl;



output << "Current Run Time: " << timer << endl;
output << "Vo = " << Vo << "  D (in units)= " << D << "  N (in units)=" << N << endl;
output << "###################" << endl;
output << "Results" << endl;
output << "###################" << endl;
output << "Total Charge on One Plate: " <<TotalQ/2 << endl;       //Total Charge (One Plate)
output << "Total Capacitance: "<< TotalQ/(Vo*2)/(100*SQRTN) << endl;    //Total Capacitance
output << "Estimated Capacitance: "<< Eo*(N/D)/(100*SQRTN) << endl;       //Estimate Capacitance - For Comparison
output << "Error Margin :"<< (TotalQ/(2*Vo))/(Eo*N/D) << endl;  //Error
log << "Each Unit is (in mm):" << 10/((float) SQRTN) << " " << "D (in mm): " << D/((float) SQRTN/10) << " " <<"Error Margin :"<< (TotalQ/(2*Vo))/(Eo*N/D) << endl;  //Error

cout << "All Complete" << endl;

output.close();

}

log.close();

return 0;

}
