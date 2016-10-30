// freeze_thaw.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "disk.h"
#include "Vector.h"
//#include "init_hopper.h"
//#include "init_silo.h"
//#include "init_window_box.h"
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<assert.h>
#include<time.h>
#include<stdlib.h>
#include<set>
#include<vector>
#include<cmath>
#include<iterator>
#include<list>
#include<stdio.h>
#include<limits>

#define length(a) ( sizeof ( a ) / sizeof ( *a ) )//Length of an array
//#define DEBUG
using namespace std;

//Function Headers
//Dump particle properties to the console
void dump_grain_to_console(int i, disk p);
//Dump particle properties to file
void dump_grain_to_file(ostream & os, int i, disk p);

//Calculate the packing density of the system
double density();
//Calculate the Total Kinetic Energy of the System
double total_kinetic_energy();
//Calculate the Total Linear Momentum of the System
double total_linear_momentum();
//Calculate the Total Angular Momentum of the System
double total_angular_momentum();
//Calculate the Average Coordination Number of the System
double avg_coordination_number();
//The Verlet List functions
bool make_verlet();
bool do_touch(int i, int k);
bool verlet_needs_update();
void predict(double dt);
void correct(double dt);

int particles;						//The number of particles in the system
int tag;							//Used to select the force law
int output;							//Used to dump the phase properties
Vector G;
vector<disk> particle;//The Gravity vector
double Z = 0.0;						//The System Coordination number
double potential_sum = 0.0;			//The sum of the system potential energy
double kinetic_sum = 0.0;			//The sum of the system kinetic energy
double linear_momentum_sum = 0.0;	//The sum of the system linear momentum
double angular_momentum_sum = 0.0;	//The sum of the system angular momentum
//Verlet variables
vector<set<int> > verlet;			//The Verlet List
vector<vector<vector<int> > > celllist;//The Grid List
vector<disk> safe;
double verlet_distance = 0.00005;
double verlet_ratio = 0.6;
double verlet_grid = 0.05;
double verlet_increase = 1.1;
double dx,dy;
int vnx,vny;
double lx = 2.0;
double ly = 2.0;
double x_0 = 0.0;
double y_0 = 0.0;
double Timesafe;

int main(int argc, char* argv[])//Main program
{	
	double dt,tmax, totalSteps;		//The timestep and the max time
	time_t start, end;	//Timer variables
	double diff;		//The difference in time
	int j;			//Counters for loops
	int counter = 0;	//Print out counter


	ofstream out1 ("Output/run1.out");	//File writer object, writes phase information at predetermined timesteps
	ofstream out2 ("Input/run2.in");	//File writer object, writes the final phase information at the end
	//ifstream in ("input.dat");//File reader object
	ifstream in ("Input/window_box.in");//File reader object, change this to run any tests
    #ifdef DEBUG
	    //Testing output file, this is for any information
	    ofstream vel ("Output/vel2.out");

	    //Output file for the energy of the system
	    ofstream eng ("Output/eng2.out");

	    //Output file for the Linear Momentum of the system
	    ofstream lmn ("Output/lmn2.out");

	    //Output file for the Angular Momentum of the system
	    ofstream amn ("Output/amn2.out");

	    //Output file for the Coordination Numberof the system
	    ofstream cdn ("Output/cdn.out");
    #endif

	if(!out1){//If the file is not accessible
		cout << "Cannot open the first file to be written to!\n";
		return 1;
	}

	if(!out2){//If the file is not accessible
		cout << "Cannot open the seccond file to be written to!\n";
		return 1;
	}

	if(!in){//If the file is not present
		cout << "Cannot access the input file!\n";
		return 1;
	}
	
	cout << "The input file is being read:\n";
	//Read data from file, Gravity
	in >> particles;
	in >> tag;
	in >> tmax;
	in >> dt;
	in >> G.x();
	in >> G.y();
	in >> G.phi();

	output = 100000;

	//Print Pre-processing details to the screen
    totalSteps = tmax/dt;
	cout << "The Simulation is of " << particles << " particles" << ".\n" 
		<< "It simulates "<< tmax << " second(s) at a time step of " << dt << ".\n"
		<< "Using force law " << tag << ".\n"
	    << "This equates to " << totalSteps << " steps. \n";
	
	particle.resize(particles);
	safe=particle;
	Timesafe=0.0;
	vnx = int(lx/verlet_grid);
	vny = int(ly/verlet_grid);
	if(vnx == 0){
		vnx = 1;
	}
	if(vny == 0){
		vny = 1;
	}
	dx = lx/vnx;
	dy = ly/vny;

	celllist.resize(vnx);
	for(int i = 0; i < vnx; i++){
		celllist[i].resize(vny);
	}
	make_verlet();

	cout << "The particles are being assigned their properties:\n";
	//Input the particle properties from a file
	for(unsigned int i = 0; i < particle.size(); i++){
		in >> particle[i].x();
		in >> particle[i].y();
		in >> particle[i].phi();
		in >> particle[i].vx();
		in >> particle[i].vy();
		in >> particle[i].omega();
		in >> particle[i].m();
		in >> particle[i].r();
		in >> particle[i].type();
		in >> particle[i].mu();
		in >> particle[i].gamma();
		in >> particle[i].Y();
		in >> particle[i].A();
		particle[i].J() = 0.0000000001;
	}
	make_verlet();

	kinetic_sum=total_kinetic_energy();
	linear_momentum_sum = total_linear_momentum();
	angular_momentum_sum = total_angular_momentum();

	
	//Particle Properties, goes to run1.out
	for(unsigned int j = 0; j < particles; j++){
		dump_grain_to_file(out1,j,particle[j]);
	}
    #ifdef DEBUG
	//Coordination Number
	    cdn << Z << "\n";
	
	//Testing output file
	    vel << particle[0].vx() << "\t" << particle[0].vy() << "\t" << particle[0].omega() << "\n";
	
	//Print the system energy
	    eng << kinetic_sum << "\t" << potential_sum << "\t" << kinetic_sum + potential_sum << "\n";

	//Print the system linear momentum
	    lmn << particle[0].linear_momentum() << "\t" << particle[1].linear_momentum() << "\t" << linear_momentum_sum << "\n";

	//Print the system linear momentum
	    amn << particle[0].angular_momentum() << "\t" << particle[1].angular_momentum() << "\t" << angular_momentum_sum << "\n";
    #endif
	//make verlet
	//Start time
	time(&start);
	//Print it to screen
	cout << "The Main loop has started:\n";
	for(double t = 0; t <= tmax; t+=dt){
        #ifdef DEBUG
		    cout << "The time is: " << t << "\n";
        #endif
		counter++;
        #ifdef DEBUG
		    cout << "Counter increased to..." << counter << "\n";
        #endif
		bool ok = true;
		bool newverlet = false;
		if(verlet_needs_update()){
            #ifdef DEBUG
			    cout << "(1)Verlet needs update...\n";
            #endif
			ok = make_verlet();
			#ifdef DEBUG
			    cout << "(1a)...\n";
			#endif
			newverlet = true;
			#ifdef DEBUG
			    cout << "(1b)...\n";
			#endif
		}
		if(!ok){
			#ifdef DEBUG
			    cout << "(2)Verlet fail! Reneging...\n";
			#endif
			particle = safe;
			#ifdef DEBUG
			    cout << "(2a)...\n";
			#endif
			t = Timesafe;
			#ifdef DEBUG
			    cout << "(2b)...\n";
			#endif
			verlet_distance *= verlet_increase;
			#ifdef DEBUG
			    cout << "(2c)...\n";
			#endif
			make_verlet();
			#ifdef DEBUG
			    cout << "(2d)...\n";
			#endif
		}
		if(newverlet && ok){
			#ifdef DEBUG
			    cout << "(3) Everything is ok...\n";
			#endif
			safe = particle;
			#ifdef DEBUG
			    cout << "(3a)...\n";
			#endif
			Timesafe = t;
			#ifdef DEBUG
			    cout << "(3b)...\n";
			#endif
		}
		#ifdef DEBUG
		    cout <<"Predicting...\n";
		#endif
		for(unsigned int i = 0; i < particle.size(); i++){
			particle[i].reset_coordination_number();
			if(particle[i].type()==0){
			    particle[i].set_force_to_zero();
				#ifdef DEBUG
			        cout << "Predict...\n";
			        dump_grain_to_console(i,particle[i]);
				#endif
			    particle[i].predict(dt);
				#ifdef DEBUG
			        dump_grain_to_console(i,particle[i]);
				#endif
			}
			else {
					particle[i].boundary_conditions(i,dt,t);
				 }
			#ifdef DEBUG
			    dump_grain_to_file(out1,i,particle[i]);
			#endif
		}

		//cout << "Cycling...\n";
		//Cycle through the particles and look for collisions
		for(unsigned int i = 0; i < particle.size(); i++){
			set<int>::iterator iter;
			for(iter = verlet[i].begin(); iter != verlet[i].end(); iter++){
				force(particle[i],particle[*iter],lx,ly,tag);
				}
			}
		//Calculate the average Coordination number of the system
		Z = avg_coordination_number();

		//Correct the positions and reset the Coordination Numbers
		//cout << "Correcting...\n";
		for(unsigned int i = 0; i < particle.size(); i++) {
			if(particle[i].type()==0){
				#ifdef DEBUG
				    cout << "Correct...\n";
				    dump_grain_to_console(i,particle[i]);
				#endif
			    particle[i].correct(dt);
				#ifdef DEBUG
			        dump_grain_to_console(i,particle[i]);
			        dump_grain_to_file(out1,i,particle[i]);
				#endif
			}
			particle[i].reset_coordination_number();
			//The system kinetic energy
			kinetic_sum = total_kinetic_energy();
			//The system linear momentum
			linear_momentum_sum = total_linear_momentum();
			//The system angular momentum
			angular_momentum_sum = total_angular_momentum();

		}

		//Output to File, goes to run1.out
		if(counter%output==0){
			cout << float(counter/totalSteps)*100.0 << " percent complete \n";
			//cout << "The time is: " << t << "\n";
		for(j = 0; j < particles; j++){
			dump_grain_to_file(out1,j,particle[j]);
		}
        #ifdef DEBUG
			//Coordination Number
			cdn << Z << "\n";

			//Testing output file
			vel << particle[0].vx() << "\t" << particle[0].vy() << "\t" << particle[0].omega() << "\n";

			//Print the system energy
			eng << kinetic_sum << "\t" << potential_sum << "\t" << kinetic_sum + potential_sum << "\n";

			//Print the System linear momentum
			lmn << particle[0].linear_momentum() << "\t" << particle[1].linear_momentum() << "\t" << linear_momentum_sum << "\n";

			//Print the System angular momentum
			amn << particle[0].angular_momentum() << "\t" << particle[1].angular_momentum() << "\t" << angular_momentum_sum << "\n";
        #endif
		}
		kinetic_sum=potential_sum=0.0;
	}
	//End Time
	time(&end);
	
	//Final Phase Plot, goes to run2.out
	/*for(j = 0; j < particles; j++){
		dump_grain_to_file(out2,j,particle[j]);
	}*/
	//Calculate the time to run the loop
	diff = difftime(end,start);
	cout << "It has taken " << diff/60 << " minutes to run the main loop.\n";

	//Calculate the system density
	//double rho = density();
	//print it to the output file
	//out << "The system density is: " << rho << "\n";

	//Close the input and output files
	in.close();
	out1.close();
	out2.close();
    #ifdef DEBUG
	    vel.close();
	    eng.close();
	    lmn.close();
	    amn.close();
	    cdn.close();
    #endif

	return 0;
}

//--------------------------------------------Functions-----------------------------------------
void dump_grain_to_console(int i, disk p){
	cout << i << "\t" <<
			p.x() <<
			"\t" << p.y() <<
			"\t" << p.phi() <<
			"\t" << p.vx() <<
			"\t" << p.vy() <<
			"\t" << p.omega() <<
            "\t" << p.type() <<
			"\t" << p.A() << "\n\n\n";
}

void dump_grain_to_file(ostream & os, int i, disk p){
						os << i << "\t" <<
								p.x() <<
								"\t" << p.y() <<
								"\t" << p.phi() <<
								"\t" << p.vx() <<
								"\t" << p.vy() <<
								"\t" << p.omega() <<
								"\t" << p.m() <<
								"\t" << p.r() <<
								"\t" << p.type() <<
								"\t" << p.mu() <<
								"\t" << p.gamma() <<
								"\t" << p.Y() <<
								"\t" << p.A() << "\n";
}

double density (){
	double sum;	//The sum of the disk areas, the system width and height
	double result;
	
			
	sum = 0.0;

	//Calculate the sum of the disk areas, from r^2 * PI
	for(unsigned int j = 0; j < particle.size(); j++){
		sum += pow((particle[j].r()),2.0)*3.14159265;
	}
	result = sum/(lx*ly);

	return result;
}

//The kinetic energy of a particle
double disk::kinetic_energy() const{
	return _m*(rtd1.x()*rtd1.x()/2 + rtd1.y()*rtd1.y()/2) + _J*rtd1.phi()*rtd1.phi()/2;
}
//The total kinetic energy of the system
double total_kinetic_energy()
{
  double sum=0.0;
  for(unsigned int i=0;i<particle.size();i++){
    if(particle[i].type()==0){
      sum+=particle[i].kinetic_energy();
    }
  }
  return sum;
}

//The Linear Momentum of a particle
double disk::linear_momentum() const{
	return _m*(rtd1.x() + rtd1.y()) + _J*rtd1.phi();
}

//The System Linear Momentum
double total_linear_momentum(){
	double sum= 0.0;

	for(unsigned int i = 0; i < particle.size(); i++){
		if(particle[i].type()==0){
			sum += particle[i].linear_momentum();
		}
	}
	return sum;
}

//The Angular Momentum of a particle
double disk::angular_momentum() const{
	return _J*rtd1.phi();
}

//The System Angular Momentum
double total_angular_momentum(){
	double sum= 0.0;

	for(unsigned int i = 0; i < particle.size(); i++){
		if(particle[i].type()==0){
			sum += particle[i].angular_momentum();
		}
	}
	return sum;
}

//The avergae Coordination number
double avg_coordination_number(){
	double sum=0.0;
	double result;

	for(unsigned int a = 0; a < particle.size(); a++){
		sum+=particle[a].z();
	}

	result = sum/particles;
	return result;
}

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
			//cout <<"_J: " << _J << "\n";
			//cout << "1/_J: " << double(1.0/double(_J)) << "\n";
			//cout <<"_force.phi(): " <<_force.phi() << "\n";

			accel=Vector((((1/_m)*_force.x())+G.x()),(((1/_m)*_force.y())+G.y()),(((1/_J)*_force.phi())+G.phi()));
			//cout << accel;
			//cout << "\n";
			
			corr=accel-rtd2;

			rtd0 += coeff0*corr;
			rtd1 += coeff1*corr;
			rtd2 = accel;
			rtd3 += coeff3*corr;
			rtd4 += coeff4*corr;
		}

void force(disk & p1, disk & p2, double lx, double ly, int tag){
			double dx=normalise(p1.x()-p2.x(),lx);
			double dy=normalise(p1.y()-p2.y(),ly);
			double rr=sqrt(dx*dx+dy*dy);
			double r1=p1.r();
			double r2=p2.r();
			double m1 = p1.m();
			double m2 = p2.m();
			double xi=r1+r2-rr;

			if(xi>0){
				p1.inc_coordination_number();
				p2.inc_coordination_number();
				potential_sum += (1*(xi*xi))/2;
				double Y=p1._Y*p2._Y/(p1._Y+p2._Y);
				double A=0.5*(p1._A+p2._A);
				double mu = (p1._mu<p2._mu ? p1._mu : p2._mu);
				double gamma = (p1._gamma<p2._gamma ? p1._gamma : p2._gamma);
				gamma=100.0;
				double gamma_t = 100.0;
				double kn = 300000.0;
				double reff = (r1*r2)/(r1+r2);
				double meff = (m1*m2)/(m1+m2);
				double dvx=p1.vx()-p2.vx();
				double dvy=p1.vy()-p2.vy();
				double rr_rez=1/rr;
				double ex=dx*rr_rez;
				double ey=dy*rr_rez;
				double xidot=-(ex*dvx+ey*dvy);
				double vtrel=-dvx*ey + dvy*ex + p1.omega()*p1.r()-p2.omega()*p2.r();
				double fn,fx,fy;
				if(tag==0){
						fn = kn*xi + gamma*meff*xidot;			//1st force law
					}
				else if(tag==1){
						fn=sqrt(xi)*Y*sqrt(reff)*(xi+A*xidot);	//2nd force law
					}
				else if(tag==3){
						fx = kn*xi - (gamma*meff*(dvx*ex));		//3rd force law
						fy = kn*xi - (gamma*meff*(dvy*ex));
					}
				double ft=-gamma_t*vtrel;

			if(fn<0) fn=0;
			if(ft<-mu*fn) ft=-mu*fn;
			if(ft>mu*fn) ft=mu*fn;
			
			if((tag==0)||(tag==1)){
				if(p1.type()==0) {
				    p1.add_force(Vector(fn*ex-ft*ey, fn*ey+ft*ex, r1*ft));
				}
			    if(p2.type()==0) {
				    p2.add_force(Vector(-fn*ex+ft*ey, -fn*ey-ft*ex, -r2*ft));
				}
			}
			else{
			    if(p1.type()==0) {
				    p1.add_force(Vector(fx*ex-ft*ey, fy*ey+ft*ex, r1*ft));
				}
			    if(p2.type()==0) {
				    p2.add_force(Vector(-fx*ex+ft*ey, -fy*ey-ft*ex, -r2*ft));
				}
			}
		}
}

void disk::boundary_conditions(int n, double dt, double t)
{
  switch(type()){
  case(0): break;
  case(1): break;
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

bool make_verlet(){
	bool ok = true;
    #ifdef DEBUG
	    cout << "Making the new Verlet list...\n";
    #endif

	verlet.resize(particles);
    #ifdef DEBUG
	    cout << "Resized...\n";
    #endif
	for(int ix = 0;ix < vnx; ix++){
		for(int iy = 0; iy < vny; iy++){
			celllist[ix][iy].clear();
		}
	}
    #ifdef DEBUG
	    cout << "Celllist cleared...\n";
    #endif
	for(unsigned int i = 0; i < particle.size(); i++){
        #ifdef DEBUG
		    cout << i << " " << particle[i].x() << " " << particle[i].y() << "...\n";
        #endif
		int ix = int((particle[i].x()-x_0)/dx);
		int iy = int((particle[i].y()-y_0)/dy);
        #ifdef DEBUG
		    cout << ix << "..." << iy << "...\n";
        #endif
		celllist[ix][iy].push_back(i);
        #ifdef DEBUG
		    cout << "Remaking...\n";
        #endif
	}
    #ifdef DEBUG
	    cout << "Cells remade...\n";
	#endif
	for(unsigned int i = 0; i < particles; i++){
		set<int> oldverlet=verlet[i];
		verlet[i].clear();
		int ix = int((particle[i].x()-x_0)/dx);
		int iy = int((particle[i].y()-y_0)/dy);
		for(int iix = ix - 1; iix <= ix + 1; iix++){
			for(int iiy = iy - 1; iiy <= iy + 1; iiy++){
				int wx = (iix+vnx)%vnx;
				int wy = (iiy+vny)%vny;
				for(unsigned int k = 0; k < celllist[wx][wy].size(); k++){
					int pk = celllist[wx][wy][k];
					if(pk < (int)i){
						if(Distance(particle[i],particle[pk],lx,ly) < particle[i].r() + particle[pk].r() + verlet_distance){
							if((particle[i].type()==0) || (particle[pk].type()==0)){
								verlet[i].insert(pk);

								if(oldverlet.find(pk)==oldverlet.end()){
									if(do_touch(i,pk)){
										ok=false;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return ok;
}
bool do_touch(int i, int k){
	return(Distance(particle[i],particle[k],lx,ly) < particle[i].r() + particle[k].r());
}

bool verlet_needs_update(){
	for(unsigned int i = 0; i < particle.size(); i++){
		if(Distance(particle[i],safe[i],lx,ly) >= verlet_ratio*verlet_distance){
			return true;
		}
	}
	return false;
}
