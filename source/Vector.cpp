#ifndef _Vector_cpp
#define _Vector_cpp

#include<iostream>
#include<math.h>
using namespace std;

class Vector {
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

	//Public---------------------------------------------
	public:
		//Explicit constructor
		explicit Vector(double x=0, double y=0, double phi=0): _x(x), _y(y), _phi(phi){};
		//explicit Vector(const Vector & ivec): _x(ivec._x),_y(ivec._y),_phi(ivec._phi){};
		
		//Access private members
		double & x(){return _x;}
		double x() const {return _x;}
		double & y(){return _y;}
		double y() const {return _y;}
		double & phi(){return _phi;}
		double phi() const {return _phi;}

		//Overloaded binary operators
		//Addition
		const Vector & operator += (const Vector & p){
			_x += p._x; _y += p._y; _phi += p._phi;
			return *this;
		}

		//Subtraction
		const Vector & operator -= (const Vector & p){
			_x -= p._x; _y -= p._y; _phi -= p._phi;
			return *this;
		}

		//Multiplication by a scalar
		const Vector & operator *= (double c){
			_x *= c; _y *= c; _phi *= c;
			return *this;
		}

		//Private--------------------------------------------
	private:
		double _x;		//x-component
		double _y;		//y-component
		double _phi;	//angular-component
};
const Vector null(0,0,0);
#endif
