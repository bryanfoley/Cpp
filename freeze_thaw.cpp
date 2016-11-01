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
//Display the progress bar
void display_progress(double progress);
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
	in.open("Input/hopper.in",ios::in);

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
		//If no roll back is needed, update the safe values
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
			//cout << float(counter/totalSteps)*100.0 << " percent complete \n";
			display_progress(double(counter/totalSteps));
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
#include "helper.hpp"
