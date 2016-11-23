#ifndef HELPER_HPP_
#define HELPER_HPP_

#include<iostream>

extern Vector G;

using namespace std;

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

void dump_grain_to_console(){
	for (unsigned int i=0; i < particle.size(); i++){
	    cout << i << ", " <<
			particle[i].x() <<
			", " << particle[i].y() <<
			", " << particle[i].phi() <<
			", " << particle[i].vx() <<
			", " << particle[i].vy() <<
			", " << particle[i].omega() <<
			", " << particle[i].m() <<
			", " << particle[i].r() <<
			", " << particle[i].type() <<
			", " << particle[i].mu() <<
			", " << particle[i].gamma() <<
			", " << particle[i].Y() <<
			", " << particle[i].A() << "\n";
	}
}

void dump_grain_to_file(ostream & os){
	for (unsigned int i=0; i< particle.size(); i++){
						os << i << "\t" <<
								particle[i].x() <<
								"\t" << particle[i].y() <<
								"\t" << particle[i].phi() <<
								"\t" << particle[i].vx() <<
								"\t" << particle[i].vy() <<
								"\t" << particle[i].omega() <<
								"\t" << particle[i].m() <<
								"\t" << particle[i].r() <<
								"\t" << particle[i].type() <<
								"\t" << particle[i].mu() <<
								"\t" << particle[i].gamma() <<
								"\t" << particle[i].Y() <<
								"\t" << particle[i].A() << "\n";
	}
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

//The average Coordination number
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
	cout << "\nIt has taken " << diff/60 << " minutes to run the main loop.\n";
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
	for(int i = 0; i < particles; i++){
		int ix = int((particle[i].x()-x_0)/dx);
		int iy = int((particle[i].y()-y_0)/dy);
		celllist[ix][iy].push_back(i);
	}
	for(int i = 0; i < particles; i++){
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

double calc_average_radius(){
	double sum_radii;
	double avg_radii;
	for (unsigned int i=0; i < particle.size(); i++){
		sum_radii +=particle[i].r();
	}
	avg_radii = ceil(sum_radii/particles);
	return avg_radii;
}

void update_verlet_variables(){
	verlet_distance = calc_average_radius();
	verlet_ratio = 0.6;
	verlet_grid = verlet_distance*1.1;
	verlet_increase = 1.1;
}

void resize_cells(vector<disk> &particle, vector<vector<vector<int> > > &celllist){
	cout << "";
	particle.resize(particles);
	cout << "";
	update_verlet_variables();
	cout << "";
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
	cout << endl;
}

void display_progress(double progress){
	if (progress < 1.0){
	    int barWidth = 70;

	    cout << "[";
	    int pos = barWidth * progress;
	    for (int i = 0; i < barWidth; ++i){
	        if (i < pos) cout << "=";
	        else if (i == pos) cout << ">";
	        else cout << " ";
	    }
	    cout << "] " << int(progress * 100.0) << " %\r";
	    cout.flush();
	}
}

void header(){
	int barWidth = 80;
	for (int i=0; i<5; i++){
			for (int j=0; j< barWidth; j++){
				cout << "#";
			}
			cout << "\n";
	}
}

string get_date(void){
   time_t now;
   int MAX_DATE = 20;
   char the_date[MAX_DATE];

   the_date[0] = '\0';

   now = time(NULL);

   if (now != -1)
   {
      strftime(the_date, MAX_DATE, "%d_%m_%Y_%H_%M_%S", gmtime(&now));
   }

   return string(the_date);
}
#endif /* HELPER_HPP_ */
