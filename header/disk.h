#ifndef _disk_h
#define _disk_h
#include<iostream>
#include<fstream>
#include "Vector.h"
using namespace std;


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
		//Calculate the average coordination number of the system
		friend double avg_coordination_number(disk iArray[]);
		//Calculate the distance between the centre of two particles
		friend double Distance(const disk & p1, const disk & p2, double lx, double ly){
				double dx=normalise(p1.rtd0.x()-p2.rtd0.x(),lx);
				double dy=normalise(p1.rtd0.y()-p2.rtd0.y(),ly);
				return sqrt(dx*dx+dy*dy);
		}
		//Input data from a file to an object
		friend istream & operator >> (istream & is, disk & p1);
		//Output data from an object to a file
		friend ostream & operator << (ostream & os, const disk & p1);
		friend ostream & operator << (ostream & os, double val);

	public:
		//Constructor
		disk();

		//Access private members
		//Position vector
		Vector & pos();	//Access a dynamic vector
		Vector pos() const;	//Access a static vector
		//x-position
		double & x();
		double x() const;
		//y-position
		double & y();
		double y() const;
		//Angle
		double & phi();
		double phi() const;
		
		//Velocity vector
		const Vector & velocity() const;
		//x-velocity
		double & vx();
		double vx() const;
		//y-velocity
		double & vy();
		double vy() const;
		//Angular velocity
		double & omega();
		double omega() const;

		//Particle properties
		//Moment of Inertia
		double & J();
		double J() const;
		//Radius
		double & r();
		double r() const;
		//Mass
		double & m();
		double m() const;
		//Particle type
		int & type();
		int type() const;
		//Material Properties
		double & mu();
		double & gamma();
		double & Y();
		double & A();
		//Coordination number
		int & z();
		int z() const;

		//Functions that operate on the _force vector of a particle
		void set_force_to_zero();
		void add_force(const Vector & f);
		void predict(double dt);
		void correct(double dt);
		void boundary_conditions(int n, double dt, double t);

		//Other functions
		void inc_coordination_number();
		void reset_coordination_number();
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
