#include<fstream>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>

void dump_particle(ostream & os, double x, double y, double phi, double vx, double vy,
					double omega, double mass, double radius, int type, double mu,
					double gamma, double Y, double A){
	os << x << "\t" << y << "\t" << phi << "\t" << vx << "\t" << vy << "\t" << omega
		 << "\t" << mass << "\t" << radius << "\t" << type << "\t" << mu << "\t" << gamma
		 << "\t" << Y << "\t" << A << "\n"; 
}
int main(){
	ofstream fout ("Input/large_hopper.in");
	fout << 1711 << "\n" << 0 << "\n" << 5 << "\n" << 0.000001 << "\n"
		<< 0 << "\t" << -9.81 << "\t" << 0 << "\n";

	//Base of the Hopper
	for(int i = 0; i < 11; i++){
		dump_particle(fout,0.45+i*0.01,0.19,0,0,0,0,1,0.005,1,0.5,10,0.0003,0.01);
	}

	//Walls of the Hopper
	for(int i = 0; i < 50; i++){
		dump_particle(fout,0.55+(i+0.5)*0.005,i*0.01+0.2,0,0,0,0,1,0.005,1,0.5,10,100000,0.01);
		dump_particle(fout,0.45-(i+0.5)*0.005,i*0.01+0.2,0,0,0,0,1,0.005,1,0.5,10,100000,0.01);
	}
	double Rmax = 0.006;
	double Rmin = 0.004;

	//Particles above the hopper
	for(int i = 0; i < 20; i++){
		for(int k = 0; k < 20; k++){
			double centerx = 0.35+0.013*i;
			double centery = 0.65+0.013*k;
			double r = Rmin + rand() * (Rmax - Rmin) / RAND_MAX;
			//double r = Rmin*Rmax/(Rmax - z*(Rmax-Rmin));
			dump_particle(fout,centerx,centery,0,0,0,0,(r*r/(Rmax*Rmax)),r,0,0.5,10,100000,0.01);
		}
	}
	fout.close();
}