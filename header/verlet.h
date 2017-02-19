#ifndef _verlet_h
#define _verlet_h
#include<iostream>
#include<fstream>
#include<vector>
#include<set>
#include "vector.h"

class verlet{
	//Constructor
	//verlet();

	//Member functions
	public:
		verlet();	
		bool make_verlet(vector<set<int> > *p_verletList, vector<vector<vector<int> > > *p_cellList);
		bool do_touch(int i, int k);
		bool verlet_needs_update();
	
	//Member variables
		double m_verlet_distance;
		double m_verlet_ratio;
		double m_verlet_grid;
		double m_verlet_increase;
		double m_dx,m_dy;
		int m_vnx,m_vny;
		double m_lx;
		double m_ly;
		double m_x_0;
		double m_y_0;
		double m_Timesafe;
};
#endif