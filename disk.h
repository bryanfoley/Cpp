#ifndef _disk_h
#define _disk_h
#include<iostream>
#include<fstream>
#include "Vector.h"
using namespace std;

extern Vector G;

inline double normalise(double dx, double L){
			while(dx< -L/2){
				dx+=L;
			}
			while(dx>=L/2){
				dx-=L;
			}
			return dx;
		}

class disk{
		//Calculate the force between two particles
		friend void force(disk & p1, disk & p2, double lx, double ly, int tag);
		//Calculate the packing density of the system
		friend double density(disk iArray[], double lx, double ly);
		friend double avg_coordination_number(disk iArray[]);
		friend double Distance(const disk & p1, const disk & p2, 
			 double lx, double ly){
				double dx=normalise(p1.rtd0.x()-p2.rtd0.x(),lx);
				double dy=normalise(p1.rtd0.y()-p2.rtd0.y(),ly);
				return sqrt(dx*dx+dy*dy);
		}
		friend istream & operator >> (istream & is, disk & p1);
		friend ostream & operator >> (ostream & os, const disk & p1);

	public:
		//Constructor
		disk():rtd0(null),rtd1(null),rtd2(null),rtd3(null),rtd4(null){}

		//Access private members
		//Position vector
		Vector & pos() {return rtd0;}	//Access a dynamic vector
		Vector pos() const {return rtd0;}	//Access a static vector
		//x-position
		double & x() {return rtd0.x();}
		double x() const {return rtd0.x();}
		//y-position
		double & y() {return rtd0.y();}
		double y() const {return rtd0.y();}
		//Angle
		double & phi() {return rtd0.phi();}
		double phi() const {return rtd0.phi();}
		
		//Velocity vector
		const Vector & velocity() const {return rtd1;}
		//x-velocity
		double & vx() {return rtd1.x();}
		double vx() const {return rtd1.x();}
		//y-velocity
		double & vy() {return rtd1.y();}
		double vy() const {return rtd1.y();}
		//Angular velocity
		double & omega() {return rtd1.phi();}
		double omega() const {return rtd1.phi();}

		//Particle properties
		//Moment of Inertia
		double & J() {return _J;}
		double J() const {return _J;}
		//Radius
		double & r() {return _r;}
		double r() const {return _r;}
		//Mass
		double & m() {return _m;}
		double m() const {return _m;}
		//Particle type
		int & type() {return _type;}
		int type() const {return _type;}
		//Material Properties
		double & mu() {return _mu;}
		double & gamma() {return _gamma;}
		double & Y() {return _Y;}
		double & A() {return _A;}
		//Coordination number
		int & z() {return _z;}
		int z() const {return _z;}

		//Functions that operate on the _force vector of a particle
		void predict(double dt);
		void add_force(const Vector & f) {_force+=f;}
		void correct(double dt);
		void set_force_to_zero() {_force=null;}
		void boundary_conditions(int n, double dt, double t);

		//Other functions
		void inc_coordination_number() {_z++;}
		void reset_coordination_number() {_z=0;}
		double kinetic_energy() const;
		double linear_momentum() const;
		double angular_momentum() const;
			
	private:
		double _r;	//Radius
		double _m;	//Mass
		double _J;	//Moment of Inertia
		int _z;		//Coordination Number
		int _type;	//Type of particle
		double _mu,_gamma,_Y,_A;//Material properties
		Vector rtd0,rtd1,rtd2,rtd3,rtd4;	//Position vector and its higher order time derivatives
		Vector _force;	//Force vector, has x, y and angular components
};
#endif
