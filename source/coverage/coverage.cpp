/***************************************************************************+
+																			+
+  Sensor Coverage Calculation   											+
+																			+		
+  Copyright (c) 2005  by  Kuo-Chuan Lee	CS, NTHU, TAIWAN				+
+  2005/10/07  Modified by Chun-Hsien Wu	CS, NTHU, TAIWAN				+
+  2018/10/28  Arranged by Chun-Hsien Wu	ITRI, TAIWAN				    +
+																			+	
****************************************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cstdlib>

using namespace std;

#define OBSTACLE 111
#define CANDIDATE 222
#define SENSOR 333

int NUM; // # of grid points
int init_sensor_num;

struct MyPoint
{
	double x;
	double y;
};
vector<vector <MyPoint> > Obstacle;
vector<MyPoint> Sensor;

struct GridPoint
{
	double x;
	double y;
	int type;
	double miss; // missing probabilty
};
vector<GridPoint> Grid;

// check if point c between point a and b
bool between(MyPoint a, MyPoint b, MyPoint c)
{
	if (((a.x - c.x) * (b.x - c.x) <= 0) && ((a.y - c.y) * (b.y - c.y) <= 0))
		return true;
	else
		return false;
}

// check if line segment (a,b) intersects line segment (c,d)
bool intersect(MyPoint a, MyPoint b, MyPoint c, MyPoint d)
{
	double cross1 = (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x); // (c-a) x (b-a)
	double cross2 = (d.x - a.x) * (b.y - a.y) - (d.y - a.y) * (b.x - a.x); // (d-a) x (b-a)
	double cross3 = (a.x - c.x) * (d.y - c.y) - (a.y - c.y) * (d.x - c.x); // (a-c) x (d-c)
	double cross4 = (b.x - c.x) * (d.y - c.y) - (b.y - c.y) * (d.x - c.x); // (b-c) x (d-c)

	if ((cross1 == 0) && between(a, b, c)) // c is on the line segment (a,b)
		return true;
	if ((cross2 == 0) && between(a, b, d)) // d is on the line segment (a,b)
		return true;
	if ((cross3 == 0) && between(c, d, a)) // a is on the line segment (c,d)
		return true;
	if ((cross4 == 0) && between(c, d, b)) // b is on the line segment (c,d)
		return true;
	if ((cross1 * cross2 < 0) && (cross3 * cross4 < 0)) // line segment (a,b) and line segment (c,d) separate each other
		return true;
	else
		return false;
}

// check if point c within line segment (a,b)
bool online(MyPoint a, MyPoint b, MyPoint c)
{
	double cross = (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x); // (c-a) x (b-a)

	if ((cross == 0) && between(a, b, c)) // check if c is on the line segment (a,b)
		return true;
	else
		return false;
}

bool test_within(MyPoint pt, int edge_num, vector<MyPoint> polygon)
{
	int count = 0;

	if (pt.x < 0 || pt.y < 0) // negative point is not allowed
		return false;

	for (int i = 0; i < edge_num; i++)
	{
		double ax = polygon[i].x;
		double ay = polygon[i].y;
		double bx = polygon[(i + 1) % edge_num].x;
		double by = polygon[(i + 1) % edge_num].y;

		// check if pt is on the edge(a,b)
		if (online(polygon[i], polygon[(i + 1) % edge_num], pt))
			return true;
		if (ay != by) // edge(a,b) is not horizontal
		{
			// check if a or b with larger y-coordinate is on L(pt,(0,pt.y))
			if ((ay == pt.y && ax < pt.x) || (by == pt.y && bx < pt.x))
			{
				if (pt.y == max(ay, by))
					count++;
			}
			else // check if edge(a,b) intersects L
			{
				MyPoint o = {0, pt.y};

				if (intersect(polygon[i], polygon[(i + 1) % edge_num], pt, o))
					count++;
			}
		}
	}
	if (count % 2 == 1)
		return true;
	else
		return false;
}

double Coverage()
{
	double max = 0.0;
	double min = 1.0;
	double sigma = 0.0;
	double avg = 0.0;
	int count = 0;
	int miss_grid = 0;

	for (int i = 0; i < NUM; i++)
	{
		if (Grid[i].type != OBSTACLE)
		{
			if (Grid[i].miss > max)
				max = Grid[i].miss;
			if (Grid[i].miss < min)
				min = Grid[i].miss;
			sigma = sigma + Grid[i].miss;
			count++;
		}
		if (Grid[i].miss == 1)
			miss_grid++;
	}
	avg = sigma / count;

	cout << endl
		 << "Total grids = " << NUM << endl;
	cout << "Obstacles = " << (int)Obstacle.size() << endl;
	cout << "Sensors = " << (int)Sensor.size() << endl;
	cout << "max miss= " << max << endl;
	cout << "min miss= " << min << endl;
	cout << "avg miss= " << avg << endl;
	cout << "miss gird= " << miss_grid << endl;

	return (1.0 - avg) * 100.0;
}

void draw_graph(string filename)
{
	fstream psfile;
	psfile.open(filename.c_str(), ios_base::out | ios_base::in | ios_base::trunc);
	psfile << "%!" << endl;

	for (int i = 0; i < NUM; i++)
	{
		if (Grid[i].miss == 1.0)
		{
			psfile << "0.5 0.0 0.5 setrgbcolor" << endl;
			psfile << Grid[i].x << " " << Grid[i].y << " " << 0.6 << " " << 0 << " " << 360 << " "
				   << "arc"
				   << " "
				   << "fill"
				   << " "
				   << "stroke" << endl;
		}
		else
		{
			psfile << (1.0 - Grid[i].miss) << " 0.0 0.0 setrgbcolor" << endl;
			psfile << Grid[i].x << " " << Grid[i].y << " " << 0.6 << " " << 0 << " " << 360 << " "
				   << "arc"
				   << " "
				   << "fill"
				   << " "
				   << "stroke" << endl;
		}
	}

	psfile << "0.0 0.0 1.0 setrgbcolor" << endl;
	for (int i = 0; i < init_sensor_num; i++)
	{
		psfile << Sensor[i].x << " " << Sensor[i].y << " " << 1 << " " << 0 << " " << 360 << " "
			   << "arc"
			   << " "
			   << "fill"
			   << " "
			   << "stroke" << endl;
	}

	psfile << "0.0 1.0 0.0 setrgbcolor" << endl;
	for (int i = init_sensor_num; i < (int)Sensor.size(); i++)
	{
		psfile << Sensor[i].x << " " << Sensor[i].y << " " << 1 << " " << 0 << " " << 360 << " "
			   << "arc"
			   << " "
			   << "fill"
			   << " "
			   << "stroke" << endl;
	}

	psfile << "0.4 0.4 0.4 setrgbcolor" << endl;
	for (int i = 0; i < (int)Obstacle.size(); i++)
	{
		psfile << "newpath" << endl;
		psfile << Obstacle[i][0].x << " " << Obstacle[i][0].y << " moveto" << endl;
		for (int j = 1; j < (int)Obstacle[i].size(); j++)
			psfile << Obstacle[i][j].x << " " << Obstacle[i][j].y << " lineto" << endl;
		psfile << "closepath" << endl
			   << "fill" << endl;
	}
	psfile << "showpage" << endl;
	psfile.close();
}

/*******************************************************************************/

int main(int argc, char *argv[])
{
	MyPoint mp;

	string filename_in = "input-out"; // default input filename
	fstream inFile;
	if (argc > 1)
		filename_in = argv[1];
	cout << "Input filename: " << filename_in.c_str() << endl
		 << endl;

	// open input file
	try
	{
		inFile.open(filename_in.c_str());
		if (!inFile)
		{
			cout << "File open error!" << endl;
			return -1;
		}
	}
	catch (string filename_in)
	{
		cout << "Fatal error: File " << filename_in.c_str() << " not found." << endl
			 << "Abort process." << endl;
		exit(1);
	}

	/************************** test parameters ****************************/
	int boundary_corner_num; // # of the sensing area's corners
	double SRange;			 // 100% sensing range
	double PRange;			 // probable sensing range
	double alpha, beta;		 // parameters of probabilistic sensor detection model
	int obstacle_num;		 // # of obstacles
	int corner_num;			 // # of the obstacle's corners
	bool output_coverage;	// display coverage?
	bool output_graph;		 // output ps coverage graph?

	// read test file =============================================================

	vector<MyPoint> Boundary;

	inFile >> boundary_corner_num; // info. of sensing area
	for (int i = 0; i < boundary_corner_num; i++)
	{
		inFile >> mp.x;
		inFile >> mp.y;
		Boundary.push_back(mp);
	}

	inFile >> SRange;
	inFile >> PRange;
	inFile >> alpha;
	inFile >> beta;

	inFile >> obstacle_num;
	for (int i = 0; i < obstacle_num; i++)
	{
		vector<MyPoint> obs; // vector used to save the position of each corner

		inFile >> corner_num;
		for (int j = 0; j < corner_num; j++)
		{
			inFile >> mp.x;
			inFile >> mp.y;
			obs.push_back(mp);
		}
		Obstacle.push_back(obs);
	}

	inFile >> output_coverage;
	inFile >> output_graph;

	inFile >> init_sensor_num; // deployed sensors
	while (inFile >> mp.x >> mp.y)
		Sensor.push_back(mp);

	// create grid points ===============================================================
	GridPoint gp;
	int X_low = (int)Boundary[0].x;
	int X_up = (int)Boundary[2].x;
	int Y_low = (int)Boundary[0].y;
	int Y_up = (int)Boundary[2].y;
	int Xdist = X_up - X_low;
	int Ydist = Y_up - Y_low;
	int Xnum = Xdist + 1;
	int Ynum = Ydist + 1;
	NUM = Xnum * Ynum; // total # of grid points

	for (int i = 0; i < Ynum; i++)
	{
		for (int j = 0; j < Xnum; j++)
		{
			gp.x = (double)j;
			gp.y = (double)i;
			gp.miss = 1.0;

			bool o_free = true;

			for (int k = 0; k < obstacle_num && o_free; k++)
			{
				MyPoint test = {gp.x, gp.y};

				if (test_within(test, (int)Obstacle[k].size(), Obstacle[k]))
					o_free = false;
			}
			if (o_free)
				gp.type = CANDIDATE;
			else
				gp.type = OBSTACLE;
			Grid.push_back(gp);
		}
	}

	// calculate sensing probabilty ==========================================================

	for (int s = 0; s < (int)Sensor.size(); s++)
	{
		int startx = (int)Sensor[s].x - (int)(SRange + PRange);
		int endx = (int)Sensor[s].x + (int)(SRange + PRange);
		int starty = (int)Sensor[s].y - (int)(SRange + PRange);
		int endy = (int)Sensor[s].y + (int)(SRange + PRange);

		if (startx < X_low)
			startx = X_low;
		if (endx > X_up)
			endx = X_up;
		if (starty < Y_low)
			starty = Y_low;
		if (endy > Y_up)
			endy = Y_up;

		for (int i = starty; i <= endy; i++)
		{
			for (int j = startx; j <= endx; j++)
			{
				if (Grid[i * Xnum + j].type != OBSTACLE)
				{
					double miss_rate;
					double d = sqrt((Sensor[s].x - j) * (Sensor[s].x - j) + (Sensor[s].y - i) * (Sensor[s].y - i));

					if (d >= SRange + PRange)
						miss_rate = 1.0;
					else
					{
						bool o_free = true;
						for (int r = 0; r < obstacle_num && o_free; r++)
						{
							for (int t = 0; t < (int)Obstacle[r].size() && o_free; t++)
							{
								MyPoint p = {(double)j, (double)i};
								MyPoint q = Sensor[s];
								MyPoint m = Obstacle[r][t];
								MyPoint n = Obstacle[r][(t + 1) % (int)Obstacle[r].size()];

								if (intersect(p, q, m, n) == true)
									o_free = false;
							}
						}
						if (o_free != true)
							miss_rate = 1.0;
						else
						{
							if (d <= SRange - PRange)
								miss_rate = 0.0;
							else
								miss_rate = 1.0 - exp(-alpha * pow((d - (SRange - PRange)) / (2 * PRange), beta));
						}
					}
					Grid[i * Xnum + j].miss = Grid[i * Xnum + j].miss * miss_rate;
				}
			}
		}
	}

	// output coverage rate
	if (output_coverage)
		cout << "coverage= " << Coverage() << endl;

	// output coverage graph
	if (output_graph)
		draw_graph(filename_in.append("_coverage.ps"));
}