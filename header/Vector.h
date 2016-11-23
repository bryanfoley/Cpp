#ifndef VECTOR_H_
#define VECTOR_H_

#include<iostream>
#include<math.h>

using namespace std;

class Vector{
	//Friend overloaded binary functions
	friend istream & operator >> (istream & is, Vector & v){
		is >> v._x >> v._y >> v._phi;
		return is;
	}

	friend ostream & operator << (ostream & os, const Vector & v){
		os << v._x << " " << v._y << " " << v._phi;
		return os;
	}
	//Vector addition
	friend Vector operator + (const Vector & v1, const Vector & v2){
		Vector res(v1);
		res+=v2;
		return res;
	}

	//Vector subtraction
	friend Vector operator - (const Vector & v1, const Vector & v2){
		Vector res(v1);
		res-=v2;
		return res;
	}

	//Scalar multiplication
	friend Vector operator * (double c, const Vector & p){
		Vector res=p;
		res *= c;
		return res;
	}

	friend Vector operator * (const Vector & p, double c){
		return c*p;
	}
public:
	//Default Constructor
	Vector();
	//Explicit constructor
	explicit Vector(double x, double y, double phi);

	//Access private members
	double & x();
	double x() const;
	double & y();
	double y() const;
	double & phi();
	double phi() const;

	//Overloaded binary operators
	//Addition
	const Vector & operator += (const Vector & p);

	//Subtraction
	const Vector & operator -= (const Vector & p);

	//Multiplication by a scalar
	const Vector & operator *= (double c);

	//Assign a vector new values
	void operator = (double c);

	//Private--------------------------------------------
private:
	double _x;		//x-component
	double _y;		//y-component
	double _phi;	//angular-component
};


#endif /* VECTOR_H_ */
