#define BEGINNING 300000
#include<Polynomial.h>
#include<LA1.h>
#include<det.h>
double EllipticE(double x);
Cdouble Rmv()
{
	double r=drand48();
	if(r<0.25)
	{
		Cdouble c(1.,1.);
		return c;
	}
	else if(0.25<=r&&r<0.5)
	{
		Cdouble c(1.,-1.);
		return c;
	}	
	else if(0.5<=r&&r<0.75)
	{
		Cdouble c(-1.,1.);
		return c;
	}
	else
	{
		Cdouble c(-1.,-1.);
		return c;
	}
}
typedef struct
{
	int n;
	int m;
}qmnumber;
class Monca
{
public:
	///Polynomial *LY;			//Landaul level
	Polynomial *JAS;		//Jastrow Factor
	Polynomial *ynm;		//Projected landau level function
	Polynomial *ynm1;
	
	Polynomial *tmpJAS;		//temp polynomial, for the buffer in Monte-Carlo
	Polynomial *tmpynm;
	Polynomial *tmpynm1;	//Baobao is tired because Bao has to deal with the case with 2 wavefunctions.
	
public:
	//Monca();
	Monca(int N_p, int N, int P);
	//~Monca();
	
	void DoJas(const Polynomial *JAS, Polynomial *&ynm, qmnumber *qm);
	void Build(Polynomial *JAS, Polynomial *&ynm, const Cdouble *z, qmnumber *qm);
	Cdouble CF_Wave( Polynomial *tynm, const Cdouble *tz);
	double Metrop(int steps, qmnumber *qm);
	Cdouble MetropMatrix(int steps, qmnumber *qm1, qmnumber *qm2);
	double Vee(int n, complex<double> *z);
	double Vbb(int n, double niu);
	double Vbe(int n, double niu, double RN, complex<double> *z);

	
	//void G_Metrop(int n_p, int n, int p, int steps);	//Ground state Metropolis. particle number, filling factor n, flux number p,(v=n/(2pn+1)), Metropolis steps
	//void Ex_Metrop(int n_p, int n, int p, int steps);	//Exciton Metropolis
	
	Cdouble *z,*r;
	int n_p,n,p;
	double niu,RN;
};
Monca::Monca(int N_p, int N, int P)
{
	n_p=N_p;
	n=N;
	p=P;
	niu=double(n)/(2.*n*n_p+1.);
	RN=sqrt(2.0*n_p/niu);
	
	z=new Cdouble[n_p];
	r=new Cdouble[n_p];
	///LY=new Polynomial [n_p*n_p];
	ynm=new Polynomial [n_p*n_p];
	ynm1=new Polynomial [n_p*n_p];
	tmpynm=new Polynomial [n_p*n_p];
	tmpynm1=new Polynomial [n_p*n_p];
	
	JAS=new Polynomial [n_p];
	tmpJAS=new Polynomial [n_p];
}
void Monca::DoJas(const Polynomial *tJAS, Polynomial *&tynm, qmnumber *qm)
{
	if(tJAS==NULL)
	{
		cout<<"Void Jastrow."<<endl;
		return;
	}
	if(tynm!=NULL)
	{
		delete []tynm;
		tynm=new Polynomial[n_p*n_p];
	}
	
	
	for(int i=0;i<n_p*n_p;i++)
	{
		tynm[i].Clear0();
		//cout<<tynm[i];
	}
	int l=0,m=0;
	int i,j;
	Polynomial tmp,tmpd;
	for(i=0;i<n_p;i++)	//i'th row
	{
		for(j=0;j<n_p;j++)	//j'th column
		{
			l=qm[i].n;
			m=qm[i].m;
			/**Now start calculating wave function**/
			/**
			for(int k=0;k<=n;k++)
			{
				double c;
				c=pow(-1.,k)*C(l+m,l-k)/Fact(k);
				Cdouble coef(c,0.);
				Polynomial tmp;
				tmp.Clear0();
				tmp.NewTerm(coef,k+m);
				tynm[i*n_p+j]=tynm[i*n_p+j]+(tmp*tJAS[j]).Deriv(k);
				
			}**/
			double c;
			c=pow(-1.,l)/Fact(l);
			Cdouble coef(c,0.);
			
			tmp.Clear0();
			tmpd.Clear0();
			tmp.NewTerm(coef,m+l);
			tmpd=tJAS[j];			
			tynm[i*n_p+j]=tmp*tmpd.Deriv(l);
			//cout<<endl<<"t="<<t<<"  ynm("<<l<<","<<m<<")J"<<j<<"="<<ynm[i*n_p+j]<<endl;
		}
	}
	
}

void Monca::Build(Polynomial *tJAS, Polynomial *&tynm, const Cdouble *tz, qmnumber *qm)	//This gives single CF_wave_function.
{
	//JAS=new Polynomial[n_p];
	//ynm=new Polynomial[n_p];
	for(int i=0;i<n_p;i++)
	{
		//tJAS[i].Clear1();
		tJAS[i].Jas(tz,n_p,i); 
	}
	DoJas(tJAS, tynm, qm);
}

Cdouble Monca::CF_Wave( Polynomial *tynm, const Cdouble *tz)	//Calculate determinant of wave-funtion: many-body function
{
	Cdouble *matrix;
	
	matrix=new Cdouble[n_p*n_p];
	double sum=0.;
	for(int i=0;i<n_p;i++)	// i'th function
	{
		for(int j=0;j<n_p;j++)	//j'th coordinate
		{
			//cout<<ynm[i]<<endl;
			matrix[i*n_p+j]=tynm[i*n_p+j].Eval(tz[j]);
			//cout<<tynm[i*n_p+j]<<endl;
		}
	}
	for(int i=0;i<n_p;i++)
	{
		sum=sum+norm(tz[i]);
	}
	sum=exp(-0.25*sum);
	Cdouble ep(sum,0.);
	//return polar(1.,0.);
	return cpxdbl_det0(matrix, n_p)*ep;
}

double Monca::Metrop(int steps, qmnumber *qm)
{
	
	cout<<"The current program only applies for the case p=1, since I didn't make the exponent of Jastrow."<<endl;
	double Energy=0;
	double count=0;
	///Cdouble ranmv;
	srand48 (time(NULL));
	for(int i=0;i<n_p;i++)
	{
		z[i]=polar((double)(double(rand())/double(RAND_MAX))*RN,(double)(double(rand())/double(RAND_MAX))*2.*PI);
		r[i]=z[i];
		//cout<<z[i]<<endl;  
		
	}
	
	
	for(int st=0;st<steps;st++)
	{	
		//for(int i=0;i<n_p;i++)cout<<"&&&&&&&&&&&&&&&&&"<<z[i]-r[i]<<endl;
		
		for(int j=0;j<n_p;j++)	//Copy to a buffer
		{
			r[j]=z[j]; 
		}	
		
		///r[st%n_p]=r[st%n_p]+polar(0.1*RN,drand48()*2.*PI);
		for(int i=0;i<n_p;i++)
		{
			//ranmv.real=0.02*drand48();
			//ranmv.imag=0.02*drand48();
			Cdouble ranmv(0.2*(drand48()-0.5),0.2*(drand48()-0.5));
			r[i]=r[i]+0.425*Rmv();
			//r[i]=r[i]+polar(0.013*RN,drand48()*2.*PI);
		}
		Build(JAS, ynm, z, qm);
		Build(tmpJAS, tmpynm, r, qm);
		//cout<<RAND_MAX<<" Step:"<<st<<":: "<<norm(CF_Wave(ynm, z))<<"  "<<norm(CF_Wave(tmpynm,r))<<endl;
		if(norm(CF_Wave(tmpynm, r)/CF_Wave(ynm, z))>drand48() &&norm(r[st%n_p])<RN*RN)
		{	
			count=count+1.0;
			//cout<<"accept:"<<count/(double)st<<endl;
			
			//z[st%n_p]=r[st%n_p];
			for(int i=0;i<n_p;i++)
		{
			z[i]=r[i];
		}
			//Build(JAS, ynm, z);
		}
	
		if(st>BEGINNING)
		{	
			Energy=Energy+Vee(n_p, z);
			if(st%500==0)
				cout<<st<<":: "<<norm(CF_Wave(ynm, z))<<"  "<<norm(CF_Wave(tmpynm,r))<<endl;	
			if(st%10000==0)
			{
				cout<<st<<":"<<count/10000<<endl;
				count=0;
			}
		}
	}
	//cout<<"Now error sofar"<<endl;
	return Energy/(steps-BEGINNING);
}

Cdouble Monca::MetropMatrix(int steps, qmnumber *qm1, qmnumber *qm2)	//Calculate matrix elements between 2 CF_wavefunctions
{
	cout<<"The current program only applies for the case p=1, since I didn't make the exponent of Jastrow."<<endl;
	Cdouble Up(0.0,0.0),Down(0.0,0.0);
	double count=0;
	srand48 (time(NULL));
	for(int i=0;i<n_p;i++)
	{
		z[i]=polar((double)(double(rand())/double(RAND_MAX))*RN,(double)(double(rand())/double(RAND_MAX))*2.*PI);
		r[i]=z[i];
	}
	for(int st=0;st<steps;st++)
	{	
		for(int j=0;j<n_p;j++)	//Copy to a buffer
		{
			r[j]=z[j]; 
		}	
		for(int i=0;i<n_p;i++)
		{
			Cdouble ranmv(0.2*(drand48()-0.5),0.2*(drand48()-0.5));
			r[i]=r[i]+0.425*Rmv();
		}
		Build(JAS, ynm, z, qm1); //cout<<CF_Wave(ynm, z);
		Build(JAS, ynm1, z, qm2); //cout<<CF_Wave(ynm1, z);
		
		//Build(tmpJAS, tmpynm, r, qm1);
		Build(tmpJAS, tmpynm1, r, qm2);
		
		//cout<<RAND_MAX<<" Step:"<<st<<":: "<<norm(CF_Wave(ynm, z))<<"  "<<norm(CF_Wave(tmpynm,r))<<endl;
		if(norm(CF_Wave(tmpynm1, r)/CF_Wave(ynm1, z))>drand48() &&norm(r[st%n_p])<RN*RN)
		{	
			count=count+1.0;
			if(st>BEGINNING)
			{
				//Up=Up+conj(CF_Wave(ynm,z)/CF_Wave(ynm1, z));
				Up=Up+polar(Vee(n_p,z),0.0)*conj(CF_Wave(ynm,z)/CF_Wave(ynm1, z));
				
			}
			for(int i=0;i<n_p;i++)
			{
				z[i]=r[i];
			}
		}
		if(st%500==0)
			cout<<st<<":: "<<norm(CF_Wave(ynm1, z))<<"  "<<norm(CF_Wave(tmpynm1,r))<<endl;	
		if(st%10000==0)
		{
			cout<<st<<":"<<count/10000<<endl;
			count=0;
		}
	}
	Up=Up*polar(1./(steps-BEGINNING),0.);	cout<<"Up="<<Up<<endl;
	
	for(int st=0;st<steps;st++)
	{	
		for(int j=0;j<n_p;j++)	//Copy to a buffer
		{
			r[j]=z[j]; 
		}	
		for(int i=0;i<n_p;i++)
		{
			Cdouble ranmv(0.2*(drand48()-0.5),0.2*(drand48()-0.5));
			r[i]=r[i]+0.425*Rmv();
			//r[i]=r[i]+polar(0.013*RN,drand48()*2.*PI);
		}
		Build(JAS, ynm, z, qm1);
		Build(JAS, ynm1, z, qm2);
		
		Build(tmpJAS, tmpynm1, r, qm2);
		
		if(norm(CF_Wave(tmpynm1, r)/CF_Wave(ynm1, z))>drand48() &&norm(r[st%n_p])<RN*RN)
		{	
			count=count+1.0;
			if(st>BEGINNING)
			{	
				Down=Down+polar(norm(CF_Wave(ynm,z)/CF_Wave(ynm1, z)),0.0);
				
			}
			for(int i=0;i<n_p;i++)
			{
				z[i]=r[i];
			}
		}
		if(st%500==0)
			cout<<st<<":: "<<norm(CF_Wave(ynm1, z))<<"  "<<norm(CF_Wave(tmpynm1,r))<<endl;	
		if(st%10000==0)
		{
			cout<<st<<":"<<count/10000<<endl;
			count=0;
		}
	
		
	}
	Down=Down*polar(1./(steps-BEGINNING),0.);	cout<<"Down="<<Down<<endl;
	return Up/sqrt(Down);
}
/***********************************/
/**Three potentials: e-e, e-b, b-b**/
/***********************************/
double Monca::Vee(int n, complex<double> *z)
{
	double V=0;
	for(int i=0;i<n;i++)
		for(int j=i+1;j<n;j++)
		{
			V=V+1./abs(z[i]-z[j]);
		}
	return V;
}
double Monca::Vbb(int n, double niu)
{
	return n*8./3./PI*sqrt(niu*n/2);
}
double Monca::Vbe(int n, double niu, double RN, complex<double> *z)
{
	double sum=0.0;
	for(int i=0;i<n;i++)
	{
		sum=sum+EllipticE((abs(z[i])/RN)*(abs(z[i])/RN));
	}
	return -sqrt(2.*niu*n)*sum*2/PI;
}

double EllipticE(double x)	//Elliptic function, for the use of the evaluation for some wave functions.
{
	return PI/2.-PI/8.*x-3.*PI/128.*x*x-5.*PI/512.*pow(x,3)-175.*PI/32768.*pow(x,4)-441.*PI/131072.*pow(x,5)-4851*PI/2097152.*pow(x,6)-14157.*PI/8388608.*pow(x,7)-2760615.*PI/2147483648.*pow(x,8)-8690825.*PI/8589934592*pow(x,9)-112285459.*PI/137438953472.*pow(x,10);

}
