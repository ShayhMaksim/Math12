

#ifndef MYLINEARALGEBRA_H_
#define MYLINEARALGEBRA_H_
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <string.h>

namespace MyLinearAlgebra {
class TMatrix;
class TQuaternion;
class TVector {
protected:
int n;
long double *data;
public:
	TVector() noexcept;
	TVector(int n) noexcept;

	TVector(long double a[],int n);

	TVector(const TVector& value);
	virtual ~TVector();
	inline long double GetItem(int i)  {return data[i];}
	inline void SetItem(int i,long double value)  {data[i]=value;}
    inline int GetSize() const noexcept  {return n;}
	inline int GetHight() const noexcept {return n-1;}
	inline long double& operator[](int i){return data[i];}
	inline const long double& operator[](int i) const {return data[i];}
	TVector& operator=(const TVector& value);
	TVector Clone();
	void Resize(int n);
	TVector operator-()const;
	TVector operator-(const TVector& a) const;
	TVector operator+(const TVector& a) const;
	TVector operator *(long double a) const;
	long double operator *(const TVector& a) const;
	TVector operator *(const TMatrix& a) const;
	TQuaternion operator *(const TQuaternion& a) const;
    TVector operator ^(const TVector& a) const;
	friend TVector operator *(long double value,const TVector& a);
    long double length() const noexcept ;
    TVector norm() noexcept ;
	void Print();
	TVector rotateByRodrigFormule(TVector& e, long double phi) const;
	TVector rotateByQuaternion(TQuaternion& Q) const;
	TVector toreturnKrilov(TQuaternion& Q) const;
    TMatrix inMatrix() const;
    TVector Concat(int i, int j) const;
    TVector Clip(TVector a);
};

class TMatrix
{
protected:



int n, m;


long double **data;

public:



TMatrix() noexcept ;



TMatrix(int n, int m) noexcept ;



TMatrix(const TMatrix& rvalue);

TMatrix(long double ** A, int n, int m);

void Print() ;



TMatrix& operator = (const TMatrix& rvalue);



virtual ~TMatrix();



inline int GetRowCount() const noexcept { return n; }



inline int GetColCount() const noexcept  { return m; }



inline int GetRowHigh() const noexcept { return n-1; }


inline int GetColHigh() const noexcept { return m-1; }



void Resize(int n, int m);



inline long double& operator()(int i, int j) { return data[i][j]; }



inline const long double& operator()(int i, int j) const { return data[i][j]; }



TMatrix operator - () const;



TMatrix operator - (const TMatrix& arg) const;


TMatrix operator + (const TMatrix& arg) const;


TMatrix operator * (long double arg) const;


TMatrix operator * (const TMatrix& arg) const;



TVector operator * (const TVector& arg) const;



friend TMatrix operator * (long double lvalue, const TMatrix& rvalue);



TMatrix operator ! () const;



long double det() const;



TMatrix t() const noexcept ;



static TMatrix E(int n) noexcept ;



TMatrix swapRows(int i, int j);
TMatrix swapCols(int i, int j);
void toText(std::string name);
TVector inVector() const;

};
class TQuaternion
{


	long double q0;


	TVector Q;

	public:


	TQuaternion();



	TQuaternion(long double l0, long double l1, long double l2, long double l3);


	TQuaternion(long double phi, const TVector& e);


	TQuaternion(const TQuaternion& rvalue);


	TQuaternion(long double yaw, long double pitch,long double roll);



	TQuaternion& operator = (const TQuaternion& rvalue);



	TQuaternion operator - (const TQuaternion& arg) const;


	TQuaternion operator + (const TQuaternion& arg) const;



	TQuaternion operator * (const TQuaternion& arg) const;



	TQuaternion operator * (const TVector& arg) const;


	TQuaternion operator * (long double arg) const;


	friend TQuaternion operator * (long double lvalue, const TQuaternion & rvalue);


	TQuaternion operator ! () const;



	inline long double scal() const noexcept { return q0; }


	inline TVector vect() const noexcept { return Q; }



	TQuaternion norm() noexcept ;



	long double Length() const noexcept ;



	TQuaternion conj() const;


	TMatrix toRotateMatrix() ;



	TQuaternion toReturnQ(const TMatrix& A) const;

    void Print();
};
} /* namespace MyLinearAlgebra */

#endif /* TVECTOR_H_ */
