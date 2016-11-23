#ifndef _Vector_cpp
#define _Vector_cpp

#include<iostream>
#include"Vector.h"
using namespace std;

Vector::Vector(){
	_x=0.0;
	_y=0.0;
	_phi=0.0;
}

Vector::Vector(double x, double y, double phi){
	_x=x;
	_y=y;
	_phi=phi;
}
		
//Access private members
double & Vector::x() {return _x;}
double Vector::x() const {return _x;}
double & Vector::y() {return _y;}
double Vector::y() const {return _y;}
double & Vector::phi() {return _phi;}
double Vector::phi() const {return _phi;}

//Overloaded binary operators
//Addition
const Vector & Vector::operator += (const Vector & p){
	_x += p._x; _y += p._y; _phi += p._phi;
	return *this;
}

//Subtraction
const Vector & Vector::operator -= (const Vector & p){
	_x -= p._x; _y -= p._y; _phi -= p._phi;
	return *this;
}

//Multiplication by a scalar
const Vector & Vector::operator *= (double c){
	_x *= c; _y *= c; _phi *= c;
	return *this;
}

//Assign a vector new values entirely
void Vector::operator = (double c){
	_x = c; _y = c; _phi = c;
}
#endif
