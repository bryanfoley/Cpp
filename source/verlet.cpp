#include "defines.h"
#include "verlet.h"
using namespace std;

void verlet::verlet():
		m_verlet_distance(0.00005),
		m_verlet_ratio(0.6),
		m_verlet_grid(0.05),
		m_verlet_increase(1.1),
		m_lx(2.0),
		m_ly(3.0),
		m_x_0(0.0),
		m_y_0(-1.5),
		m_vnx(int(m_lx/m_verlet_grid)),
		/*if(vnx == 0){
			vnx = 1;
		}
		if(vny == 0){
			vny = 1;
		}*/
		m_vny(int(m_ly/m_verlet_grid)),
		m_dx(m_lx/m_vnx),
		m_dy(m_ly/m_vny),
		m_Timesafe(0.0)
{
	cout << "Created a Verlet object\n";
}

bool verlet::make_verlet(vector<set<int> > *p_verletList, vector<vector<vector<int> > > *p_cellList){
	bool ok = true;
	//cout << "Making the new Verlet list...\n";

	*p_verletList->resize(particles);
	//cout << "Resized...\n";
	for(int ix = 0;ix < m_vnx; ix++){
		for(int iy = 0; iy < m_vny; iy++){
			*p_cellList[ix][iy].clear();
		}
	}
	//cout << "*p_cellList cleared...\n";
	for(unsigned int i = 0; i < particle.size(); i++){
		int ix = int((particle[i].x()-m_x_0)/m_dx);
		int iy = int((particle[i].y()-m_y_0)/m_dy);
		*p_cellList[ix][iy]->push_back(i);
		//cout << "Remaking...\n";
	}
	//cout << "Cells remade...\n";
	for(unsigned int i = 0; i < particles; i++){
		set<int> oldverletList=*p_verletList[i];
		*p_verletList[i].clear();
		int ix = int((particle[i].x()-m_x_0)/m_dx);
		int iy = int((particle[i].y()-m_y_0)/m_dy);
		for(int iix = ix - 1; iix <= ix + 1; iix++){
			for(int iiy = iy - 1; iiy <= iy + 1; iiy++){
				int wx = (iix+m_vnx)%m_vnx;
				int wy = (iiy+m_vny)%m_vny;
				for(unsigned int k = 0; k < *p_cellList[wx][wy].size(); k++){
					int pk = *p_cellList[wx][wy][k];
					if(pk < (int)i){
						if(Distance(particle[i],particle[pk],m_lx,m_ly) < particle[i].r() + particle[pk].r() + m_verlet_distance){
							if((particle[i].type()==0) || (particle[pk].type()==0)){
								*p_verletList[i]->insert(pk);

								if(oldverletList.find(pk)==oldverletList.end()){
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

bool verlet::do_touch(int i, int k){
	return(Distance(particle[i],particle[k],lx,ly) < particle[i].r() + particle[k].r());
}

bool verlet::verlet_needs_update(){
	for(unsigned int i = 0; i < particle.size(); i++){
		if(Distance(particle[i],safe[i],m_lx,m_ly) >= m_verlet_ratio*m_verlet_distance){
			return true;
		}
	}
	return false;
}