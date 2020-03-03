
#include "MyLinearAlgebra.h"




namespace MyLinearAlgebra {


TVector::TVector() noexcept {

	this->data=new long double[1];
	n=1;
}

TVector::TVector(int n) noexcept
{
	this->data=new long double[n];
	this->n=n;
	for (int i=0;i<n;i++)
	{
	   data[i]=0;
	}
}


TVector::TVector(long double *a,int n)
{
	this->n=n;
	this->data=new long double[n];
			for(int j=0;j<n;j++)
			{
				data[j]=a[j];
			}
}

TVector::TVector(const TVector& value):n(1),data(nullptr)
{
	(*this)=value;
}

TVector& TVector::operator=(const TVector& value)
{
	if (this!=&value)
	{
		if(n!=value.n)
		{
			if (data) {delete[] data;}
			data=new long double[value.n];
			n=value.n;
		}
       memcpy(data,value.data,sizeof(long double)*n);
	}
	return (*this);
}

TVector TVector::Clone()
{
	TVector vect(this->GetSize());

					for(int i=0;i<this->GetSize();i++)
					{
						vect.SetItem(i,this->GetItem(i));
					}
					return vect;
}

void TVector::Resize(int n)
{
	long double *a=new long double[this->GetSize()];

					  for(int i=0;i<this->GetSize();i++)
					  {
						  a[i]=data[i];
					  }
	                delete[] data;
	               data=nullptr;
	                this->n = n;

			        this->data=new long double[n];
			        for(int i=0;i<n;i++)
			        {
			         if (i<this->GetSize())
			  		   {
			  		   this->SetItem(i,  a[i]);
			  		   }
			        }
                    delete[] a;
}


TVector operator *(long double value,const TVector& a)
{
	   TVector V(a.n);
		for (int i=0;i<a.n;i++)
		{
            V[i]=value*a[i];
		}
		return V;
}
TVector TVector::operator -(const TVector& a) const
{
#ifdef _DEBUG
	if(n!=a.n) throw 1;
#endif
TVector V(n);
for (int i=0;i<n;i++)
{
	V[i]=data[i]-a[i];
}
return V;
}

TVector TVector::operator +(const TVector& a) const
{
#ifdef _DEBUG
	if(n!=a.n) throw 1;
#endif

TVector V(n);
for (int i=0;i<n;i++)
{
	V[i]=data[i]+a[i];
}
return V;
}

TVector TVector::operator *(long double a) const
{
	TVector V(n);
	for (int i=0;i<n;i++)
	{
		V[i]=data[i]*a;
	}
	return V;
}



long double TVector::operator *(const TVector& a) const
{
#ifdef _DEBUG
	if(n!=a.n) throw 1;
#endif
long double b=0;
for (int i=0;i<n;i++)
{
	b+=data[i]*a[i];
}
return b;

}

TVector TVector::operator*(const TMatrix& a) const
{
#ifdef _DEBUG
	if(n!=a.n) throw 1;
#endif
    TVector V(n);

    for (int i = 0; i < n; i++)
    {
        V[i] = 0;
        for (int j = 0; j < a.GetColCount(); j++)
            V[i] += a(i,j) * data[j];
    }

    return V;
}

TVector TVector::operator ^ (const TVector& a) const
{
#ifdef _DEBUG
	if(n!=a.n) throw 1;
#endif
	TVector Cross(n);
	Cross[0]= data[1]*a[2]-data[2]*a[1];
	Cross[1]= data[2]*a[0]-data[0]*a[2];
    Cross[2]= data[0]*a[1]-data[1]*a[0];
			return Cross;
}


long double TVector::length() const noexcept
{
           double b=0;
			for(int i=0;i<n;i++)
	        {
	       	 b+=data[i]*data[i];
	        }
			return pow(b, 0.5);
}

TVector TVector::norm() noexcept
{
	       TVector vect(this->GetSize());
			for(int i=0; i<this->GetSize();i++)
			   {
				  vect[i]=(data[i]/this->length());
			   }
			return vect;
}

void TVector::Print()
{
	std::cout << "( ";
			for (int i=0; i<GetSize(); i++) {
				std::cout << std::to_string(GetItem(i))+" ";
			}
			std::cout << ")" << std::endl;
}

TVector TVector::operator-()const
{
	TVector Neg = (*this);
			for (int i=0; i<Neg.GetSize(); i++) {

					Neg.SetItem(i, Neg.GetItem(i)*(-1));

			}
			return Neg;
}

TMatrix TVector::inMatrix() const
{
	try
			{
				if ((int)(sqrt(this->GetSize()))!=sqrt(this->GetSize())) throw 1;
			}
			catch (std::exception& ex)
			{
			//	this->Resize((int)pow((int)round(sqrt(this->GetSize())),2));
			}
			TMatrix A(sqrt(this->GetSize()),sqrt(this->GetSize()));
			for(int j =0;j<(int)sqrt(this->GetSize());j++)
					{
						for(int i =0;i<(int)sqrt(this->GetSize());i++)
							{
					            A(i, j)=data[A.GetColCount()*j+i];
							}
					}
			return A;
}

TVector TVector::Concat(int i, int j) const
{
	TVector vect(j-i+1);
			for(int l=i;l<=j;l++)
			{
                vect[l-i]= data[l];
			}
			return vect;
}

TVector TVector::Clip(TVector a)
{
	TVector vect((int)(this->GetSize()+a.GetSize()));
			for(int i=0;i<this->GetSize();i++)
			{
				vect[i]= this->GetItem(i);
			}
			for(int j=this->GetSize();j<vect.GetSize();j++)
			{
				vect[j]= a[j-this->GetSize()];
			}
			return vect;
}

TVector TVector::rotateByRodrigFormule(TVector& e, long double phi) const
{
	        TVector p1(this->GetSize());
			TQuaternion Q(phi,e);
			TVector theta(3);
			 theta[0] = 2*Q.vect()[0]/Q.scal();
			 theta[1] = 2*Q.vect()[1]/Q.scal();
			 theta[2] = 2*Q.vect()[2]/Q.scal();
			 p1=(*this)+theta*(1/(1+pow(tan(phi/2), 2)));
			 p1=p1 ^ ((*this)+theta ^ (*this)* 0.5);
			    return p1;
}

TVector TVector::rotateByQuaternion(TQuaternion& Q) const
{
	        TVector p1(3);
			TQuaternion t;
			t=Q*(*this)*(!Q);
			p1=t.vect();
		    return p1;
}

TVector TVector::toreturnKrilov(TQuaternion& Q) const
{

			long double yaw=atan((2*Q.scal()*Q.vect()[1]-2*Q.vect()[0]*Q.vect()[2])/(2*Q.scal()*Q.scal()+2*Q.vect()[0]*Q.vect()[0]-1));
			long double pitch=asin((2*Q.vect()[0]*Q.vect()[1]+2*Q.scal()*Q.vect()[2]));
			long double roll=atan((2*Q.scal()*Q.vect()[0]-2*Q.vect()[1]*Q.vect()[2])/(2*Q.scal()*Q.scal()+2*Q.vect()[1]*Q.vect()[1]-1));
			long double *A=new long double[3] {yaw,pitch,roll};
			TVector p(A,3);
			return p;
}

TVector::~TVector() {

	if (data){
		delete[] data;
		n=0;
		data=nullptr;
	}
}



TMatrix::TMatrix() noexcept  : n(1), m(1), data(nullptr) {}


TMatrix::TMatrix(int n, int m) noexcept  : n(n),m(m), data(nullptr)
{
	data=new long double*[n];
	for(int i=0;i<n;i++)
	{
		data[i]=new long double[m];
	}
	for(int i=0;i<n;i++)
		{
		for(int j=0;j<m;j++)
			data[i][j]=0;
		}
}


TMatrix::TMatrix(const TMatrix& rvalue) : n(1), m(1), data(nullptr)  {
(*this) = rvalue;
}

TMatrix::TMatrix(long double ** A, int n, int m)
{
    this->n=n;
    this->m=m;

	data=new long double*[n];
		for(int i=0;i<n;i++)
		{
			data[i]=new long double[m];
		}

			for (int i=0;i<n;i++)
			{
				for(int j=0;j<m;j++)
				{
					data[i][j]=A[i][j];
				}
			}
}



TMatrix& TMatrix::operator = (const TMatrix& rvalue) {

if (this != &rvalue) {

this->~TMatrix();

this->Resize(rvalue.n, rvalue.m);

for (int i = 0; i < n; i++)
memcpy(this->data[i], rvalue.data[i], sizeof(long double)*m);
}

return (*this);
}

void TMatrix::Resize(int n, int m)
{

if ((this->n==n) & (this->m==m)) return;

if ((this->n==1) && (this->m==1))
{
	data=new long double*[1];
		for(int j=0;j<1;j++)
		{
			data[j]=new long double[1];
		}
}

long double arg[this->GetRowCount()][this->GetColCount()];



for(int i=0;i<this->GetRowCount();i++)
	       {
	    	   for(int j=0;j<this->GetColCount();j++)
	    	   {
	    		   arg[i][j]=data[i][j];
	    	   }
	       }


for (int k = 0; k < this->n; k++)
delete[] data[k];
delete[] data;
data = NULL;

data=new long double*[n];
	for(int j=0;j<n;j++)
	{
		data[j]=new long double[m];
	}
	       for(int i=0;i<n;i++)
	       {
	    	   for(int j=0;j<m;j++)
	    	   {
	    		   if (i<this->GetRowCount() && j<this->GetColCount())
	    		   {
	    		   data[i][j]=arg[i][j];
	    		   }
	    		   else data[i][j]=0;
	    	   }
	       }
	       this->n = n;
	       this->m = m;
}



TMatrix TMatrix::E(int n) noexcept {
TMatrix E(n,n);
for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++)
		{
		    if(i==j)
		    { E(i,j) = 1; }
		}
}
return E;
}



TMatrix TMatrix::operator + (const TMatrix& arg) const {
#ifdef _DEBUG
if ((n != arg.n) || (m != arg.m))
throw 1;
#endif
TMatrix M(n, m);
for (int i = 0; i < M.GetRowCount(); i++)
{
	for (int j = 0; j < M.GetColCount(); j++)
	{M(i,j) = data[i][j] + arg(i,j);}
}
	return M;
	}

TMatrix TMatrix::operator - (const TMatrix& arg) const {
#ifdef _DEBUG
if ((n != arg.n) || (m != arg.m))
throw 1;
#endif
TMatrix M(n, m);
for (int i = 0; i < M.GetRowCount(); i++)
{
	for (int j = 0; j < M.GetColCount(); j++)
	{M(i,j) = data[i][j] - arg(i,j);}
}
	return M;
	}

TVector TMatrix::operator * (const TVector& arg) const
		{
#ifdef _DEBUG
if ((n != arg.n))
throw 1;
#endif
TVector V(n);

for (int i = 0; i < n; i++)
{
    V[i] = 0;
    for (int j = 0; j < m; j++)
        V[i] += data[i][j] * arg[j];
}

return V;
		}

TMatrix TMatrix::operator * (const TMatrix& arg) const
		{
#ifdef _DEBUG
if ((n != arg.n) || (m != arg.m))
throw 1;
#endif
TMatrix M(n, arg.m);

for (int i = 0; i < M.n; i++)
    for (int j = 0; j < M.m; j++)
    {
        M(i, j) = 0;
        for (int k = 0; k < m; k++)
            M(i,j) += data[i][k] * arg(k, j);
    }

return M;
		}

TMatrix TMatrix::operator * (long double arg) const
		{
	TMatrix matrix(this->GetRowCount(),this->GetColCount());
			for(int i=0;i<this->GetRowCount();i++)
			{
				for(int j=0;j<this->GetColCount();j++)
				{
						matrix(i, j)= (data[i][j]*arg);
				}
			}
			return matrix;
		}

TMatrix operator * (long double lvalue, const TMatrix& rvalue)
		{
	TMatrix matrix(rvalue.GetRowCount(),rvalue.GetColCount());
				for(int i=0;i<rvalue.GetRowCount();i++)
				{
					for(int j=0;j<rvalue.GetColCount();j++)
					{
							matrix(i, j)= (rvalue(i,j)*lvalue);
					}
				}
				return matrix;
		}


TMatrix TMatrix::operator - () const
{
	TMatrix Neg = (*this);
			for (int i=0; i<Neg.GetRowCount(); i++) {
				for (int j=0; j<Neg.GetColCount(); j++) {
						Neg(i,j) = Neg(i,j)*(-1);
				}
			}
			return Neg;
}


void TMatrix::Print()
{
	for (int i=0; i<GetRowCount(); i++) {
				for (int j=0; j<GetColCount(); j++) {
					std::cout <<std::to_string(data[i][j])+" ";
				}
				std::cout << (" ") <<std::endl;
			}
}

TMatrix TMatrix::t() const noexcept
{
	TMatrix mt2(this->GetColCount(),this->GetRowCount());
			for(int i=0;i<this->GetRowCount();i++)
			{
				for(int j=0;j<this->GetColCount();j++)
				{
					mt2(j, i)= data[i][j];
				}
			}
			return mt2;
}


TMatrix TMatrix::swapRows(int i, int j)
{
	TMatrix A=(*this);
	            long double vrem;
				for(int d=0;d<this->GetColCount();d++)
				{
					vrem = A(i,d);
					A(i,d) = A(j,d);
					A(j,d)=vrem;
				}
				return A;
}

TMatrix TMatrix::swapCols(int i, int j)
{
	TMatrix A=(*this);
	long double vrem;

				for(int d=0;d<this->GetColCount();d++)
				{
					vrem = A(d, i);
					A(d, i)= A(d, j);
					A(d, j)= vrem;
				}
				return A;
}

long double TMatrix::det() const
{
	        int count;
			long double result=1;

			TMatrix mt1=(*this);

			TMatrix mt2=(*this);


			for(int i=0;i<this->GetColCount();i++)
			{
				if (mt1(i, i)==0 && i<this->GetRowCount()-1)
				{
					count=i;
					do {
						count++;
						if (mt1(i, i)==0 && i<this->GetColCount()-1)
						{
					    count=i;
						do {
						count++;
						} while (mt1(count, 0)==0);
						mt1=mt1.swapCols(i, count);

						result=-result;
						}
                    } while (mt1(count, 0)==0);
					mt1=mt1.swapRows(i, count);

					result=-result;
				}

				for(int x=0;x<this->GetRowCount();x++)
				{
					for(int y=0;y<this->GetColCount();y++)
					{
						mt2(x, y)= mt1(x, y);
					}
				}

				long double tmp;
				for(int j=i+1;j<this->GetRowCount();j++)
				{
					for(int k=0;k<this->GetColCount();k++)
					{
						tmp=mt1(j, k)-((mt1(i, k)*mt2(j, i))/mt1(i, i));
						mt1(j, k)= tmp;
					}
				}

			}
			for(int x=0;x<this->GetRowCount();x++)
			{
				result=result*mt1(x, x);
			}
			return result;
}


TMatrix TMatrix::operator !() const
		{
	double const1 = 0.00001;
	try{
			if (this->GetRowCount()!=this->GetColCount()) {
				 throw 1;
			}
	}
	catch (const std::exception& err)
	{
         std::cout << "The end"<< std::endl;
	}
			TMatrix x(this->GetRowCount(),this->GetColCount());
			TMatrix A = (*this);
			std::vector<int> list1;
			std::vector<int> list2;
			for (int i=0; i<x.GetRowCount(); i++) {
				if (A(i, i)<const1 && A(i, i)>-const1) {
					int k=1;
					long double max = A(i, i);
					int numRows = i;
					for (k = 1; i+k<A.GetRowCount(); k++) {
						if (!(A(i+k, i)<const1 && A(i+k, i)>-const1) && abs(A(i+k, i))>abs(max)) {
							max = A(i+k, i);
							numRows = i+k;
						}
					}
					A=A.swapRows(i, numRows);
					list1.push_back(i);
					list2.push_back(numRows);
					}
				}
			for (int i=0; i<x.GetRowCount(); i++) {
				x(i, i)= 1;
			}

			for (int i=0; i<x.GetRowCount(); i++) {
				TMatrix r = A;
				for (int j=0; j<x.GetColCount(); j++) {
					x(i, j)= x(i, j)/r(i,i);
					A(i, j)= A(i, j)/r(i,i);
				}

				r=A;
				for (int j=0; j<x.GetColCount(); j++) {
					for (int k=0; k<x.GetColCount(); k++) {
						if (k!=i) {
							A(k, j)= A(k, j)-r(k, i)*A(i, j);
							x(k, j)= x(k, j)-r(k, i)*x(i, j);
						}
					}
				}
			}


			if (list1.size()>0) {
				for (int i=0; i<(int)list1.size(); i++) {
					x=x.swapCols(list1[list1.size()-1-i], list2[list2.size()-1-i]);
				}
			}

			return x;
		}

TVector TMatrix::inVector() const
{
	TVector vect(this->GetRowCount()*this->GetColCount());
			for(int j=0;j<this->GetRowCount();j++)
			{
				for(int i=0;i<this->GetRowCount();i++)
				{
					vect[this->GetColCount()*j+i]=data[i][j];
				}
			}
			return vect;
}

void TMatrix::toText(std::string name)
{
     std::ofstream out;
     out.open(name+".txt");
     if (out.is_open())
     {
            out.precision(16);
    		for (int i=0;i<this->GetRowCount()-1;i++)
    				{
    						for(int j=0;j<this->GetColCount();j++)
    						{
                                out<< this->data[i][j];
    							out<<(" ");
    						}
    						out <<std::endl;

    				}
    				out.close();
     }
}



TMatrix::~TMatrix()
{
if (data) {
for (int i = 0; i < n; i++)
delete[] data[i];
delete[] data;
data = nullptr;
n = m = 0;
}
}

TQuaternion::TQuaternion(): q0(0),Q(3){}


TQuaternion::TQuaternion(long double l0, long double l1, long double l2, long double l3): q0(l0),Q(3)
{
	Q[0]= l1;
	Q[1]= l2;
	Q[2]= l3;
}

TQuaternion::TQuaternion(long double phi, const TVector& e)
{
	            try {
                if (e.GetSize()!=3)  throw "Bad size";
				this->q0=cos(phi/2);
		    	Q=e;
		    	Q=Q.norm();
		    	Q=Q*(sin(phi/2));
				}
				catch(std::exception& ex)
				{
					ex.what();
				}
}

TQuaternion::TQuaternion(long double yaw,long double pitch,long double roll)

    {
    	this->q0=cos(roll/2)*cos(yaw/2)*cos(pitch/2)-sin(roll/2)*sin(yaw/2)*sin(pitch/2);
    	TVector Q(3);
    	Q.SetItem(0, sin(roll/2)*cos(yaw/2)*cos(pitch/2)+cos(roll/2)*sin(yaw/2)*sin(pitch/2));
    	Q.SetItem(1, cos(roll/2)*sin(yaw/2)*cos(pitch/2)+sin(roll/2)*cos(yaw/2)*sin(pitch/2));
    	Q.SetItem(2, cos(roll/2)*cos(yaw/2)*sin(pitch/2)-sin(roll/2)*sin(yaw/2)*cos(pitch/2));
    }


TQuaternion::TQuaternion(const TQuaternion& rvalue) {*this=rvalue;}


TQuaternion& TQuaternion::operator = (const TQuaternion& rvalue)
{
    if(this!=&rvalue)
    {
    	q0=rvalue.q0;
    	Q=rvalue.Q;
    }
    return *this;
}

TQuaternion TQuaternion::operator - (const TQuaternion& arg) const
		{
	        TQuaternion Q;
			Q.q0=q0-arg.q0;
			Q.Q=this->Q-arg.Q;
			return Q;
		}

TQuaternion TQuaternion::operator + (const TQuaternion& arg) const
		{
	  	  	  	TQuaternion Q;
				Q.q0=this->q0+arg.q0;
				Q.Q=this->Q+arg.Q;
				return Q;
		}

TQuaternion TQuaternion::operator * (const TQuaternion& arg) const
		{
	        TQuaternion C;
			C.q0=q0*arg.q0-Q*arg.Q;
			C.Q=Q*arg.q0+arg.Q*q0+Q^arg.Q;
			return C;
		}

TQuaternion TQuaternion::operator * (const TVector& arg) const
		{
	        TQuaternion Q=*this;
			TQuaternion Q1(0,arg[0],arg[1],arg[2]);
			Q=Q*Q1;
			return Q;
		}

TQuaternion TQuaternion::operator * (long double arg) const
		{
	        TQuaternion Q=*this;
            Q.q0=Q.q0*arg;
            Q.Q=Q.Q*arg;
			return Q;
		}

TQuaternion operator * (long double lvalue, const TQuaternion & rvalue)
		{
	            TQuaternion Q;
				Q.q0=rvalue.q0*lvalue;
				Q.Q=rvalue.Q*lvalue;
				return Q;
		}

TQuaternion TQuaternion::norm() noexcept
{
	  TQuaternion Q;
	  long double b=this->Length();
			Q.q0=q0*1/b;
			Q.Q=this->Q*(1/b);
			return Q;
}

long double TQuaternion::Length() const noexcept
{
	        long double a=q0*q0;
			for(int i=0;i<3;i++)
			a+=Q[i]*Q[i];
			return pow(a, 0.5);
}

TQuaternion TQuaternion::conj() const
{
	        TQuaternion Q;
			Q.q0=this->scal();
			Q.Q=this->vect()*(-1);
			return Q;
}

TQuaternion TQuaternion::operator ! () const
		{
	        return this->conj()*(1/(this->Length() * this->Length()));
		}

TMatrix TQuaternion::toRotateMatrix()
{
	TMatrix A(3,3);

			TQuaternion Q=*this;
			Q=Q.norm();
			A(0, 0)= Q.scal()*Q.scal()+Q.Q[1]*Q.Q[1]-Q.Q[2]*Q.Q[2]-Q.Q[3]*Q.Q[3];
			A(1, 1)= Q.scal()*Q.scal()-Q.Q[1]*Q.Q[1]+Q.Q[2]*Q.Q[2]-Q.Q[3]*Q.Q[3];
			A(2, 2)= Q.scal()*Q.scal()-Q.Q[1]*Q.Q[1]-Q.Q[2]*Q.Q[2]+Q.Q[3]*Q.Q[3];


			A(0, 1)= 2*(Q.Q[1]*Q.Q[2]-Q.scal()*Q.Q[3]);
			A(0, 2)= 2*(Q.scal()*Q.Q[2]+Q.Q[1]*Q.Q[3]);

			A(1, 0)= 2*(Q.Q[1]*Q.Q[2]+Q.scal()*Q.Q[3]);

			A(2, 0)= 2*(Q.Q[1]*Q.Q[3]-Q.scal()*Q.Q[2]);

			A(1, 2)= 2*(Q.Q[2]*Q.Q[3]-Q.scal()*Q.Q[1]);
			A(2, 1)= 2*(Q.scal()*Q.Q[1]+Q.Q[2]*Q.Q[3]);
			return A;
}

void TQuaternion::Print()
{
   std::cout << "(";
   	   std::cout << scal()<< " ";
   	   	   for(int i=0;i<3;i++)
   	   	   {
   	   		   std::cout<< vect()[i] <<" ";
   	   	   }
   	   	   	   std::cout<<")"<<std::endl;
}

TQuaternion TQuaternion::toReturnQ(const TMatrix& A) const
{
	try {
		long double q0=pow((1+A(0, 0)+A(1, 1)+A(2, 2)),0.5)/2;
		long double q1=(A(2, 1)-A(1, 2))/(4*q0);
		long double q2=(A(0, 2)-A(2, 0))/(4*q0);
		long double q3=(A(1, 0)-A(0, 1))/(4*q0);
		TQuaternion Q(q0,q1,q2,q3);
		return Q;
		}
		catch (std::exception& ex)
		{
            std::cout << "Can't" <<std::endl;
			TQuaternion Q;
			return Q;
		}
}



}


