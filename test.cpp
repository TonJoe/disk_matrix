//#include <Polynomial.h>
//#include "Polynomial.h"
//#include "LA1.h"
#include "Monca.h"

using namespace std;
void test()  
{  
    Polynomial p1;  
	Cdouble c1(3.,0.);
	Cdouble c2(4.,0);
	Cdouble d1(3.,4);
	p1.NewTerm(c1,0);
	p1.NewTerm(c2,1);
    Polynomial p2;	 
    p2.NewTerm(d1,0);  //cout<<p2;

    cout<<"("<<p1<<") * ("<<p2<<") = "<<p1*p2<<endl;  
}  

int main()
{
	Monca experiment(6,3,1);
	qmnumber qm1[6],qm2[6];
	qm1[0].n=0; qm1[0].m=0;
	qm1[1].n=0; qm1[1].m=1;
	qm1[2].n=0; qm1[2].m=2;
	qm1[3].n=1; qm1[3].m=-1;
	qm1[4].n=1; qm1[4].m=0;
	qm1[5].n=1; qm1[5].m=1;
	
	qm2[0].n=0; qm2[0].m=0;
	qm2[1].n=0; qm2[1].m=1;
	qm2[2].n=0; qm2[2].m=2;
	qm2[3].n=0; qm2[3].m=3;
	qm2[4].n=1; qm2[4].m=-1;
	qm2[5].n=2; qm2[5].m=-2;
	//cout<<experiment.Metrop(2000000,qm2)<<endl;
	cout<<"6,3,1   "<<experiment.MetropMatrix(1000000,qm2, qm1)<<endl;
}
