// freeze_thaw.cpp : Defines the entry point for the console application.
//************************************************************************//
//*******************************Inclusions*******************************//
//************************************************************************//
#include "stdafx.h"
#include "disk.h"
#include "Vector.h"
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

//************************************************************************//
//*********************************Defines********************************//
//************************************************************************//
#define length(a) ( sizeof ( a ) / sizeof ( *a ) )//Length of an array
//#define DEBUG
using namespace std;

//************************************************************************//
//****************************Function Headers****************************//
//************************************************************************//
//Assign the system properties from a file
void assign_system_properties(istream & in);
//Assign the system properties from a file
void assign_particle_properties(istream & in);
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
//Calculate the time taken to run the simulation
void calc_execution_time(time_t startTime, time_t endTime);
//The Verlet List functions
bool make_verlet();
bool do_touch(int i, int k);
bool verlet_needs_update();
void resize_cells(vector<disk> &particle, vector<vector<vector<int> > > &celllist);

//************************************************************************//
//****************************System Variables****************************//
//************************************************************************//
int particles;						//The number of particles in the system
int tag;							//Used to select the force law
int output;							//Used to dump the phase properties
double dt,tmax, totalSteps;			//The timestep,max time and total steps in the simulation
time_t start, end;					//Timer variables
double diff;						//The difference in time
int counter = 0;					//Print out counter

Vector G;							//The Gravity vector
vector<disk> particle;				//The vector that holds all particle objects
double Z = 0.0;						//The System Coordination number
double potential_sum = 0.0;			//The sum of the system potential energy
double kinetic_sum = 0.0;			//The sum of the system kinetic energy
double linear_momentum_sum = 0.0;	//The sum of the system linear momentum
double angular_momentum_sum = 0.0;	//The sum of the system angular momentum

//************************************************************************//
//****************************Verlet Variables****************************//
//************************************************************************//
vector<set<int> > verlet;			//The Verlet List
vector<vector<vector<int> > > celllist;//The Grid List
vector<disk> safe;					//Backup of the vector that holds all particle objects

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

//************************************************************************//
//**********************************Main**********************************//
//************************************************************************//

int main(int argc, char* argv[])//Main program
{	
	//Open files for input/output
	ofstream out1;					//File writer object, writes phase information at predetermined timesteps
	out1.open("Output/run1.bin",ios::out | ios::binary);

	ifstream in;//File reader object, change this to run any tests
	in.open("Input/grv_test.in",ios::in);

	if(!out1){//If the file is not accessible
		cerr << "Cannot open the first file to be written to!\n";
		return 1;
	}

	if(!in){//If the file is not present
		cerr << "Cannot access the input file!\n";
		return 1;
	}
	int j;								//Counters for loops

	assign_system_properties(in);

	//Print Pre-processing details to the screen
    totalSteps = tmax/dt;
	output = int(totalSteps/100);
	cout << "The Simulation is of " << particles << " particles" << ".\n" 
		<< "It simulates "<< tmax << " second(s) at a time step of " << dt << ".\n"
		<< "Using force law " << tag << ".\n"
	    << "This equates to " << totalSteps << " steps. \n"
	    << "Every " << output << " steps there will be a file write operation.\n";
	
	resize_cells(particle,celllist);

	assign_particle_properties(in);

	//Close the input files
	in.close();
	make_verlet();

	kinetic_sum=total_kinetic_energy();
	linear_momentum_sum = total_linear_momentum();
	angular_momentum_sum = total_angular_momentum();

	
	//Particle Properties, goes to run1.out
	for(unsigned int j = 0; j < particles; j++){
		dump_grain_to_file(out1,j,particle[j]);
	}
	
	//Start time
	time(&start);
	//Print it to screen
	cout << "The Main loop has started:\n";
	for(double t = 0; t <= tmax; t+=dt){
		counter++;
		bool ok = true;
		bool newverlet = false;
		//Check if the verlet list must be updated
		if(verlet_needs_update()){
			ok = make_verlet();
			newverlet = true;
		}
		//If error detected, roll back to the last known correct positions and remake the verlet list.
		if(!ok){
			particle = safe;
			t = Timesafe;
			verlet_distance *= verlet_increase;
			make_verlet();
		}
		//Roll the simulation back to the last correct state if needed
		if(newverlet && ok){
			safe = particle;
			Timesafe = t;
		}
		//Predict the new properties of the free particles, or apply boundary conditions
		for(unsigned int i = 0; i < particle.size(); i++){
			particle[i].reset_coordination_number();
			if(particle[i].type()==0){
			    particle[i].set_force_to_zero();
			    particle[i].predict(dt);
			}
			else {
					particle[i].boundary_conditions(i,dt,t);
				 }
		}
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
		for(unsigned int i = 0; i < particle.size(); i++) {
			if(particle[i].type()==0){
			    particle[i].correct(dt);
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
			for(j = 0; j < particles; j++){
			    dump_grain_to_file(out1,j,particle[j]);
		    }
		}
		kinetic_sum=potential_sum=0.0;
	}
	//End Time
	time(&end);
	
	//Calculate the time to run the loop
	calc_execution_time(start,end);

	//Close the output files
	out1.close();

	return 0;
}

//************************************************************************//
//********************************Functions*******************************//
//************************************************************************//

void assign_system_properties(istream & in){
	cout << "Assigning the system properties...";
	in >> particles;
	in >> tag;
	in >> tmax;
	in >> dt;
	in >> G.x();
	in >> G.y();
	in >> G.phi();
	cout << "OK\n";
}

void assign_particle_properties(istream & in){
	cout << "The particles are being assigned their properties...";
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
	cout<< "OK\n";
}

void dump_grain_to_console(int i, disk p){
	cout << i << ", " <<
			p.x() <<
			", " << p.y() <<
			", " << p.phi() <<
			", " << p.vx() <<
			", " << p.vy() <<
			", " << p.omega() <<
			", " << p.m() <<
			", " << p.r() <<
			", " << p.type() <<
			", " << p.mu() <<
			", " << p.gamma() <<
			", " << p.Y() <<
			", " << p.A() << "\n";
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

	result = double(sum/particles);
	return result;
}

void calc_execution_time(time_t startTime,time_t endTime){
	//Calculate the time to run the loop
	diff = difftime(endTime,startTime);
	cout << "It has taken " << diff/60 << " minutes to run the main loop.\n";
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
				//gamma=100.0;
				double kn = 3.0;
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
				double ft=-gamma*vtrel;

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

bool make_verlet(){
	bool ok = true;
	verlet.resize(particles);

	for(int ix = 0;ix < vnx; ix++){
		for(int iy = 0; iy < vny; iy++){
			celllist[ix][iy].clear();
		}
	}
	for(unsigned int i = 0; i < particle.size(); i++){
		int ix = int((particle[i].x()-x_0)/dx);
		int iy = int((particle[i].y()-y_0)/dy);
		celllist[ix][iy].push_back(i);
	}
	for(unsigned int i = 0; i < particles; i++){
		set<int> oldverlet=verlet[i];
		verlet[i].clear();
		int ix = int((particle[i].x()-x_0)/dx);
		int iy = int((particle[i].y()-y_0)/dy);
		if (ix < 0 || iy < 0){
		}
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
void resize_cells(vector<disk> &particle, vector<vector<vector<int> > > &celllist){
	cout << "Resizing the particle list and the cell lists...";
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
	//make_verlet();
	cout << "OK\n";
}
