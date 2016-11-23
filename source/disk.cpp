#include "disk.h"
#include <stdlib.h>

Vector G;

//Constructor
disk::disk(){
	_r = 0.0;	//Radius
	_m = 0.0;	//Mass
	_J = 0.0;	//Moment of Inertia
	_z = 0;		//Coordination Number
	_type = 0;	//Type of particle
	_mu = 0.0;	//Material properties
	_gamma = 0.0;
	_Y = 0.0;
	_A = 0.0;
}

//Access private members
//Position vector
Vector & disk::pos() {return rtd0;}	//Access a dynamic vector
Vector disk::pos() const {return rtd0;}	//Access a static vector
//x-position
double & disk::x() {return rtd0.x();}
double disk::x() const {return rtd0.x();}
//y-position
double & disk::y() {return rtd0.y();}
double disk::y() const {return rtd0.y();}
//Angle
double & disk::phi() {return rtd0.phi();}
double disk::phi() const {return rtd0.phi();}

//Velocity vector
const Vector & disk::velocity() const {return rtd1;}
//x-velocity
double & disk::vx() {return rtd1.x();}
double disk::vx() const {return rtd1.x();}
//y-velocity
double & disk::vy() {return rtd1.y();}
double disk::vy() const {return rtd1.y();}
//Angular velocity
double & disk::omega() {return rtd1.phi();}
double disk::omega() const {return rtd1.phi();}

//Particle properties
//Moment of Inertia
double & disk::J() {return _J;}
double disk::J() const {return _J;}
//Radius
double & disk::r() {return _r;}
double disk::r() const {return _r;}
//Mass
double & disk::m() {return _m;}
double disk::m() const {return _m;}
//Particle type
int & disk::type() {return _type;}
int disk::type() const {return _type;}
//Material Properties
double & disk::mu() {return _mu;}
double & disk::gamma() {return _gamma;}
double & disk::Y() {return _Y;}
double & disk::A() {return _A;}
//Coordination number
int & disk::z() {return _z;}
int disk::z() const {return _z;}

//Functions that operate on the _force vector of a particle
void disk::set_force_to_zero() {_force=0.0;}
void disk::add_force(const Vector & f) {_force+=f;}

//Other functions
void disk::inc_coordination_number() {_z++;}
void disk::reset_coordination_number() {_z=0;}

//Gears Integration Scheme - Predict
void disk::predict(double dt){
			double a1 = dt;
			double a2 = a1 * dt/2;
			double a3 = a2 * dt/3;
			double a4 = a3 * dt/4;

			rtd0 += a1*rtd1 + a2*rtd2 + a3*rtd3 + a4*rtd4;
			rtd1 += a1*rtd2 + a2*rtd3 + a3*rtd4;
			rtd2 += a1*rtd3 + a2*rtd4;
			rtd3 += a1*rtd4;
		}

//Gears Integration Scheme - Correct
void disk::correct(double dt){
			static Vector accel, corr;

			double dt_recip = 1/dt;

			const double coeff0 = double(19)/double(90)*(dt*dt/double(2));
			const double coeff1 = double(3)/double(4)*(dt/double(2));
			const double coeff3 = double(1)/double(2)*(double(3)*dt_recip);
			const double coeff4 = double(1)/double(12)*(double(12)*(dt_recip*dt_recip));

			accel=Vector((((1/_m)*_force.x())+G.x()),(((1/_m)*_force.y())+G.y()),(((1/_J)*_force.phi())+G.phi()));

			corr=accel-rtd2;

			rtd0 += coeff0*corr;
			rtd1 += coeff1*corr;
			rtd2 = accel;
			rtd3 += coeff3*corr;
			rtd4 += coeff4*corr;
		}

//Apply boundary conditions based on Particle type
void disk::boundary_conditions(int n, double dt, double t)
{
  switch(type()){
  case(0): break; //Free particle
  case(1): break; //Wall particle
  case(2): {
      x()=0.5-0.4*cos(10*t);
      y()=0.1;
      vx()=10*0.4*sin(t);
      vy()=0;
    } break;
  case(3): {
      double xx=x()-0.5;
      double yy=y()-0.5;
      double xp=xx*cos(dt)-yy*sin(dt);
      double yp=xx*sin(dt)+yy*cos(dt);
      x()=0.5+xp;
      y()=0.5+yp;
      vx()=-yp;
      vy()= xp;
      omega()=1;
    } break;
  case(4): {
      x()=0.5+0.1*cos(t) + 0.4*cos(t+2*n*3.14159265/128);
      y()=0.5+0.1*sin(t) + 0.4*sin(t+2*n*3.14159265/128);
      vx()=-0.1*sin(t) - 0.4*sin(t+2*n*3.14159265/128);
      vy()= 0.1*cos(t) - 0.4*cos(t+2*n*3.14159265/128);
      omega()=1;
    } break;
  case(5): {//Vertical Vibrating Particle
	  y()=y()+((0.02*sin(30*t))*0.1);
      vx()=0;
      vy()=0.02*30*cos(30*t);
    } break;
  case(6): {
      int i=n/2;
      y()=i*0.02+0.1+0.02*sin(30*t);
      vx()=0;
      vy()=0.02*30*cos(30*t);
    } break;
  case(7): {//Expander
	  r()+=0.001;
	}break;
  case(8): {//Contractor
	  r()-=0.001;
	}break;
  default: {
      cout << "ptype: " << type() << " not implemented\n";
      abort();
    }
  }
}

//The kinetic energy of a particle
double disk::kinetic_energy() const{
	return _m*(rtd1.x()*rtd1.x()/2 + rtd1.y()*rtd1.y()/2) + _J*rtd1.phi()*rtd1.phi()/2;
}

//The Angular Momentum of a particle
double disk::angular_momentum() const{
	return _J*rtd1.phi();
}

//The Linear Momentum of a particle
double disk::linear_momentum() const{
	return _m*(rtd1.x() + rtd1.y()) + _J*rtd1.phi();
}
