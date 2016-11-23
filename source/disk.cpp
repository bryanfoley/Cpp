#include "disk.hpp"
#include <stdlib.h>

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
