/***************************************************************************+
+  DT-Score Algorithm                                                       +
+																			+		
+  Copyright (c) 2005  by  Kuo-Chuan Lee	CS, NTHU, TAIWAN				+
+  2005/10/07  Modified by Chun-Hsien Wu	CS, NTHU, TAIWAN				+	
+  2006/08/21  Modified by Chun-Hsien Wu	CS, NTHU, TAIWAN				+
+  2018/10/28  Arranged by Chun-Hsien Wu	ITRI, TAIWAN				    +
+                                                                           +
+  The Delaunay triangulation algorithm (Delaunay-tree) is based on the:    +
+                                                                           +
+  O. Devillers, Improved incremental randomized Delaunay triangula-        +
+  tion, in: Proceedings of the Fourteenth Annual Symposium on              +
+  Computational Geometry, 1998, pp. 106–115.                               +
+                                                                           +
/***************************************************************************+
+                                                                           +
+  Delaunay-tree                                                            +
+                                                                           +
+                                                                           +
+  Copyright (c) 1993  by  INRIA Prisme Project                             +
+  2004 route des Lucioles BP 93 06902 Sophia Antipolis Cedex               +
+  All rights reserved.                                                     +
+                                                                           +
****************************************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <cstdlib>
#include <sys/time.h>

using namespace std;

// measure execution time
#define MICRO_PER_SECOND 1000000
#define MICRO_PER_MILLISECOND 1000

// output file
fstream psfile;
const float OFFSET = 50.0;							   // offset for drawing ps graph
unsigned int init_sensor_num = 0, grid_sensor_num = 0; // initial sensors

struct MyPoint
{
	double x;
	double y;
};
vector<vector <MyPoint> > Obstacle;
vector<MyPoint> Boundary;
vector<MyPoint> Sensor;
struct MyNode
{
	struct MyPoint v0;
	struct MyPoint v1;
	struct MyPoint v2;
	struct MyPoint center;
	double radius;
};
vector<MyNode> Node; // center of empty circles

/*******************************************************************************/

class Delaunay_tree;
class DT_node;
class DT_list;
class point;
typedef unsigned char idx; // used for flag and index in array

vector<DT_node *> DT_node_Collection;
vector<DT_list *> DT_list_Collection;
vector<point *> point_Collection;

struct MyPoint mp;
struct MyNode mn;

class Delaunay_tree
{
  private:
	int nb;		   // number of the current operation
	DT_node *root; // the root of the delaunay_tree
  public:
	Delaunay_tree();					// initialization as empty
	~Delaunay_tree();					// destructor
	Delaunay_tree &operator+=(point *); // insertion
	void output();						// the Delaunay triangulation
	void find_nodes();
};

class point
{
  private:
	double x;
	double y;

  public:
	point(double xx, double yy)
	{
		x = xx;
		y = yy;
	}
	double X() { return x; }
	double Y() { return y; }
	void lineto(point *);
	friend point operator+(point, point);
	friend point operator-(point, point);
	friend double operator*(point, point);
	friend double operator^(point, point);
};

class DT_flag
{
  private:
	idx f;
	DT_flag() { f = (idx)0; }
	void infinite(int i) { f |= (idx)i; }
	void last_finite() { f |= 8; }
	void kill() { f |= 16; }
	idx is_infinite() { return f & 7; }
	idx is_last_finite() { return f & 8; }
	idx is_dead() { return f & 16; }
	idx is_degenerate() { return f & 32; }

  public:
	friend class DT_node;
	friend class Delaunay_tree;
};

class DT_list
{
  private:
	DT_list *next;
	DT_node *key;
	DT_list(DT_list *l, DT_node *k)
	{
		next = l;
		key = k;
	}
	~DT_list();

  public:
	friend class DT_node;
	friend class Delaunay_tree;
};

class DT_node
{
  private:
	// the first veretx is the creator, that is finite
	// except for the root and its neighbors
	DT_flag flag;
	unsigned int nb;
	point *vertices[3];
	DT_node *neighbors[3];
	DT_list *sons;
	DT_node();						  // initialize the root
	DT_node(DT_node *, idx);		  // initialize nowhere
	DT_node(DT_node *, point *, idx); // father, creator, direction of stepfather
	idx conflict(point *);			  // true if the point is inside the (closed) circumdisk
	DT_node *find_conflict(point *p); // return an alive node in conflict
	void output();
	void find_nodes();
	idx cw_neighbor_index(point *p)
	{
		return ((p == vertices[0]) ? 2 : ((p == vertices[1]) ? 0 : 1));
	}
	idx neighbor_index(DT_node *n)
	{
		return ((neighbors[0] == n) ? 0 : ((neighbors[1] == n) ? 1 : 2));
	}

  public:
	friend class Delaunay_tree;
	double getV0X() { return vertices[0]->X(); }
	double getV0Y() { return vertices[0]->Y(); }
	double getV1X() { return vertices[1]->X(); }
	double getV1Y() { return vertices[1]->Y(); }
	double getV2X() { return vertices[2]->X(); }
	double getV2Y() { return vertices[2]->Y(); }
};

inline point operator+(point a, point b) { return point(a.x + b.x, a.y + b.y); }
inline point operator-(point a, point b) { return point(a.x - b.x, a.y - b.y); }
inline double operator*(point a, point b) { return a.x * b.x + a.y * b.y; }
inline double operator^(point a, point b) { return a.x * b.y - a.y * b.x; }

DT_list::~DT_list()
{
}

idx DT_node::conflict(point *p)
{
	switch (flag.is_infinite())
	{
	case 4:
		return 0;
	case 3:
		return 1;
	case 2:
		return ((*p - *vertices[0]) * (*vertices[1] + *vertices[2]) >= 0);
	case 1:
		return ((flag.is_last_finite())
					? (((*p - *vertices[2]) ^ (*vertices[2] - *vertices[0])) >= 0)
					: (((*p - *vertices[0]) ^ (*vertices[0] - *vertices[1])) >= 0));
	case 0:
		// compute the det 4*4 column: x,y,x**2+y**2,1 for p and vertices [0,1,2]
		double x, y;
		double x0, y0;
		double x1, y1;
		double x2, y2;
		double z1, z2;
		double alpha, beta, gamma;

		x = p->X();
		y = p->Y();
		x0 = vertices[0]->X();
		y0 = vertices[0]->Y();
		x1 = vertices[1]->X();
		y1 = vertices[1]->Y();
		x2 = vertices[2]->X();
		y2 = vertices[2]->Y();
		x1 -= x0;
		y1 -= y0;
		x2 -= x0;
		y2 -= y0;
		x -= x0;
		y -= y0;
		z1 = (x1 * x1) + (y1 * y1);
		z2 = (x2 * x2) + (y2 * y2);
		alpha = (y1 * z2) - (z1 * y2);
		beta = (x2 * z1) - (x1 * z2);
		gamma = (x1 * y2) - (y1 * x2);
		return ((alpha * x) + (beta * y) + gamma * ((x * x) + (y * y)) <= 0);
	}
}

void point::lineto(point *p)
{
	psfile << this->x + OFFSET << " " << this->y + OFFSET << " moveto " << p->x + OFFSET << " " << p->y + OFFSET << " lineto stroke\n";
}

DT_node::DT_node()
{
	vertices[0] = new point(1, 0);
	point_Collection.push_back(vertices[0]);
	vertices[1] = new point(-0.5, 0.8660254);
	point_Collection.push_back(vertices[1]);
	vertices[2] = new point(-0.5, -0.8660254);
	point_Collection.push_back(vertices[2]);
	flag.infinite(3);
	;
	nb = 0;
	sons = NULL;
}

DT_node::DT_node(DT_node *root, idx i)
{
	vertices[0] = root->vertices[0];
	vertices[1] = root->vertices[1];
	vertices[2] = root->vertices[2];
	flag.infinite(4);
	nb = 0;
	sons = NULL;
	neighbors[i] = root;
	root->neighbors[i] = this;
}

DT_node::DT_node(DT_node *f, point *c, idx i)
{
	// the triangle is created in ccw order
	// circumdisk and flatness are not computed
	switch (f->flag.is_infinite())
	{
	case 0:
		flag.infinite(0);
		break;
	case 1:
		if (f->flag.is_last_finite())
			flag.infinite((i == 1) ? 0 : 1);
		else
			flag.infinite((i == 2) ? 0 : 1);
		if (flag.is_infinite())
			if (f->flag.is_last_finite())
			{
				if (i == 0)
					flag.last_finite();
			}
			else
			{
				if (i == 1)
					flag.last_finite();
			}
		break;
	case 2:
		flag.infinite((i == 0) ? 2 : 1);
		if (i == 1)
			flag.last_finite();
		break;
	case 3:
		flag.infinite(2);
		break;
	}
	nb = 0;
	sons = NULL;
	f->sons = new DT_list(f->sons, this);
	DT_list_Collection.push_back(f->sons);
	f->neighbors[i]->sons = new DT_list(f->neighbors[i]->sons, this);
	DT_list_Collection.push_back(f->neighbors[i]->sons);
	f->neighbors[i]->neighbors[f->neighbors[i]->neighbor_index(f)] = this;
	vertices[0] = c;
	neighbors[0] = f->neighbors[i];
	switch (i)
	{
	case 0:
		vertices[1] = f->vertices[1];
		vertices[2] = f->vertices[2];
		break;
	case 1:
		vertices[1] = f->vertices[2];
		vertices[2] = f->vertices[0];
		break;
	case 2:
		vertices[1] = f->vertices[0];
		vertices[2] = f->vertices[1];
		break;
	}
}

DT_node *DT_node::find_conflict(point *p)
{
	DT_list *l;
	DT_node *n;
	if (!conflict(p))
		return NULL;
	if (!flag.is_dead())
		return this;
	for (l = sons; l; l = l->next)
		if (l->key->nb != nb)
		{
			l->key->nb = nb;
			n = l->key->find_conflict(p);
			if (n)
				return n;
		}
	return NULL;
}

void DT_node::output()
{
	DT_list *l;
	if (flag.is_dead())
	{
		for (l = sons; l; l = l->next)
			if (l->key->nb != nb)
			{
				l->key->nb = nb;
				l->key->output();
			}
		return;
	}
	if (neighbors[0]->nb != nb)
		if (!flag.is_infinite())
			vertices[1]->lineto(vertices[2]);
	if (neighbors[1]->nb != nb)
		if (!flag.is_infinite())
			vertices[2]->lineto(vertices[0]);
		else if ((flag.is_infinite() == 1) && (flag.is_last_finite()))
			vertices[2]->lineto(vertices[0]);
	if (neighbors[2]->nb != nb)
		if (!flag.is_infinite())
			vertices[0]->lineto(vertices[1]);
		else if ((flag.is_infinite() == 1) && (!flag.is_last_finite()))
			vertices[0]->lineto(vertices[1]);
}

// find empty circle from Delaunay triangle
void DT_node::find_nodes()
{
	DT_list *l;
	if (flag.is_dead())
	{
		for (l = sons; l; l = l->next)
		{
			if (l->key->nb != nb)
			{
				l->key->nb = nb;
				l->key->find_nodes();
			}
		}
		return;
	}
	if (!flag.is_infinite())
	{
		double x0 = vertices[0]->X();
		double y0 = vertices[0]->Y();
		double x1 = vertices[1]->X();
		double y1 = vertices[1]->Y();
		double x2 = vertices[2]->X();
		double y2 = vertices[2]->Y();

		double alpha;
		double beta;
		double gamma;
		double alpha_div;
		double beta_div;

		// calculate the center & radius of the circumcircle of a triangle
		alpha_div = ((y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2));
		beta_div = ((x1 - x2) * (y0 - y1) - (x0 - x1) * (y1 - y2));

		alpha = (((y0 - y1) * ((x1 * x1) - (x2 * x2) + (y1 * y1) - (y2 * y2)) - (y1 - y2) * ((x0 * x0) - (x1 * x1) + (y0 * y0) - (y1 * y1)))) / alpha_div;
		beta = (((x0 - x1) * ((x1 * x1) - (x2 * x2) + (y1 * y1) - (y2 * y2)) - (x1 - x2) * ((x0 * x0) - (x1 * x1) + (y0 * y0) - (y1 * y1)))) / beta_div;
		gamma = -(alpha * x0) - (beta * y0) - (x0 * x0) - (y0 * y0);

		mp.x = vertices[0]->X();
		mp.y = vertices[0]->Y();
		mn.v0 = mp;
		mp.x = vertices[1]->X();
		mp.y = vertices[1]->Y();
		mn.v1 = mp;
		mp.x = vertices[2]->X();
		mp.y = vertices[2]->Y();
		mn.v2 = mp;
		mn.center.x = -alpha / 2;
		mn.center.y = -beta / 2;
		mn.radius = sqrt(((alpha * alpha) + (beta * beta)) / 4 - gamma);
		Node.push_back(mn);
	}
}

// output the Delaunay triangulation
void Delaunay_tree::output()
{

	root->nb = ++nb;
	root->output();
}

void Delaunay_tree::find_nodes()
{
	root->nb = ++nb;
	root->find_nodes();
}

Delaunay_tree::Delaunay_tree() // initialization as empty
{
	nb = 0;
	root = new DT_node();
	DT_node_Collection.push_back(root);
	DT_node *a, *b, *c;
	a = new DT_node(root, 0);
	DT_node_Collection.push_back(a);
	b = new DT_node(root, 1);
	DT_node_Collection.push_back(b);
	c = new DT_node(root, 2);
	DT_node_Collection.push_back(c);
	root->neighbors[0]->neighbors[1] = root->neighbors[1];
	root->neighbors[0]->neighbors[2] = root->neighbors[2];
	root->neighbors[1]->neighbors[0] = root->neighbors[0];
	root->neighbors[1]->neighbors[2] = root->neighbors[2];
	root->neighbors[2]->neighbors[0] = root->neighbors[0];
	root->neighbors[2]->neighbors[1] = root->neighbors[1];
}

Delaunay_tree::~Delaunay_tree() // destructor
{
	nb++;
	for (int i = 0; i < (int)DT_node_Collection.size(); i++)
		delete DT_node_Collection[i];
	for (int i = 0; i < (int)point_Collection.size(); i++)
		delete point_Collection[i];
	for (int i = 0; i < (int)DT_list_Collection.size(); i++)
		delete DT_list_Collection[i];
	DT_node_Collection.clear();
	point_Collection.clear();
	DT_list_Collection.clear();
}

Delaunay_tree &Delaunay_tree::operator+=(point *p) // insertion
{
	DT_node *n;
	DT_node *created;
	DT_node *last;
	DT_node *first;
	point *q;
	point *r;
	idx i;

	point_Collection.push_back(p);

	root->nb = ++nb;
	if (!(n = root->find_conflict(p)))
		return *this;
	// test if p is already inserted
	for (i = 0; (int)i < 3 - (int)n->flag.is_infinite(); i++)
		if ((p->X() == n->vertices[i]->X()) && (p->Y() == n->vertices[i]->Y()))
			return *this;
	n->flag.kill();
	// we will turn cw around first vertex of n, till next triangle
	// is not in conflict
	q = n->vertices[0];
	while (n->neighbors[i = n->cw_neighbor_index(q)]->conflict(p))
	{
		n = n->neighbors[i];
		n->flag.kill();
	}

	first = last = new DT_node(n, p, i);
	DT_node_Collection.push_back(first);
	// we will turn cw around r, till next triangle is not in conflict
	r = n->vertices[(((int)i) + 2) % 3];
	while (1)
	{
		i = n->cw_neighbor_index(r);
		if (n->neighbors[i]->flag.is_dead())
		{
			n = n->neighbors[i];
			continue;
		}
		if (n->neighbors[i]->conflict(p))
		{
			n = n->neighbors[i];
			n->flag.kill();
			continue;
		}
		break;
	}

	while (1)
	{
		// n is killed by p
		// n->neighbors[i] is not in conflict with p
		// r is vertex i+1 of n
		created = new DT_node(n, p, i);
		DT_node_Collection.push_back(created);
		created->neighbors[2] = last;
		last->neighbors[1] = created;
		last = created;
		r = n->vertices[(((int)i) + 2) % 3]; // r turn in n ccw
		if (r == q)
			break;
		// we will turn cw around r, till next triangle is not in conflict
		while (1)
		{
			i = n->cw_neighbor_index(r);
			if (n->neighbors[i]->flag.is_dead())
			{
				n = n->neighbors[i];
				continue;
			}
			if (n->neighbors[i]->conflict(p))
			{
				n = n->neighbors[i];
				n->flag.kill();
				continue;
			}
			break;
		}
	}
	first->neighbors[2] = last;
	last->neighbors[1] = first;
	return *this;
}

// check if line segment (a,b) intersects line segment (c,d)
bool intersect(MyPoint a, MyPoint b, MyPoint c, MyPoint d)
{
	double cross1 = (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x); // (c-a) x (b-a)
	double cross2 = (d.x - a.x) * (b.y - a.y) - (d.y - a.y) * (b.x - a.x); // (d-a) x (b-a)
	double cross3 = (a.x - c.x) * (d.y - c.y) - (a.y - c.y) * (d.x - c.x); // (a-c) x (d-c)
	double cross4 = (b.x - c.x) * (d.y - c.y) - (b.y - c.y) * (d.x - c.x); // (b-c) x (d-c)

	// c is on the line segment (a,b)
	if ((cross1 == 0) && ((a.x - c.x) * (b.x - c.x) <= 0) && ((a.y - c.y) * (b.y - c.y) <= 0))
		return true;
	// d is on the line segment (a,b)
	if ((cross2 == 0) && ((a.x - d.x) * (b.x - d.x) <= 0) && ((a.y - d.y) * (b.y - d.y) <= 0))
		return true;
	// a is on the line segment (c,d)
	if ((cross3 == 0) && ((c.x - a.x) * (d.x - a.x) <= 0) && ((c.y - a.y) * (d.y - a.y) <= 0))
		return true;
	// b is on the line segment (c,d)
	if ((cross4 == 0) && ((c.x - b.x) * (d.x - b.x) <= 0) && ((c.y - b.y) * (d.y - b.y) <= 0))
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

	// check if c is on the line segment (a,b)
	if ((cross == 0) && ((a.x - c.x) * (b.x - c.x) <= 0) && ((a.y - c.y) * (b.y - c.y) <= 0))
		return true;
	else
		return false;
}

bool test_within(MyPoint pt, int edge_num, vector<MyPoint> polygon)
{
	int count = 0;
	double ax, ay, bx, by;

	if (pt.x < 0 || pt.y < 0) // negative point is not allowed
		return false;

	for (int i = 0; i < edge_num; i++)
	{
		ax = polygon[i].x;
		ay = polygon[i].y;
		bx = polygon[(i + 1) % edge_num].x;
		by = polygon[(i + 1) % edge_num].y;

		if (online(polygon[i], polygon[(i + 1) % edge_num], pt)) // check if pt is on the edge(a,b)
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

void draw_graph(string filename)
{
	psfile.open(filename.c_str(), ios_base::out | ios_base::in | ios_base::trunc);
	psfile << "%!\n";

	// sense area
	psfile << "0.95 0.95 0.95 setrgbcolor" << endl;
	psfile << "newpath" << endl;
	psfile << Boundary[0].x + OFFSET << " " << Boundary[0].y + OFFSET << " moveto" << endl;
	for (int i = 1; i < (int)Boundary.size(); i++)
	{
		psfile << Boundary[i].x + OFFSET << " " << Boundary[i].y + OFFSET << " lineto" << endl;
	}
	psfile << "closepath" << endl
		   << "fill" << endl;

	// draw obstacles
	psfile << "0.2 0.5 0.5 setrgbcolor" << endl;
	for (int i = 0; i < (int)Obstacle.size(); i++)
	{
		psfile << "newpath" << endl;
		psfile << Obstacle[i][0].x + OFFSET << " " << Obstacle[i][0].y + OFFSET << " moveto" << endl;
		for (int j = 1; j < (int)Obstacle[i].size(); j++)
			psfile << Obstacle[i][j].x + OFFSET << " " << Obstacle[i][j].y + OFFSET << " lineto" << endl;
		psfile << "closepath" << endl
			   << "fill" << endl;
	}

	// create and draw DT graph
	psfile << "0.8 0.8 0.8 setrgbcolor\n";
	Delaunay_tree DTfinal;

	for (unsigned int i = 0; i < (unsigned int)Sensor.size(); i++)
		DTfinal += new point(Sensor[i].x, Sensor[i].y);

	DTfinal.output();

	// draw initial sensor point
	psfile << "0.0 0.0 0.5 setrgbcolor\n";
	for (unsigned int i = 0; i < init_sensor_num; i++)
		psfile << Sensor[i].x + OFFSET << " " << Sensor[i].y + OFFSET << " 2 0 360 arc fill stroke" << endl;

	// draw grid sensor point
	psfile << "0.0 0.0 0.0 setrgbcolor\n";
	for (unsigned int i = init_sensor_num; i < init_sensor_num + grid_sensor_num; i++)
		psfile << Sensor[i].x + OFFSET << " " << Sensor[i].y + OFFSET << " 2 0 360 arc fill stroke" << endl;

	// draw new sensor point
	psfile << "0.0 1.0 0.0 setrgbcolor\n";
	for (unsigned int i = init_sensor_num + grid_sensor_num; i < (unsigned int)Sensor.size(); i++)
		psfile << Sensor[i].x + OFFSET << " " << Sensor[i].y + OFFSET << " 2 0 360 arc fill stroke" << endl;

	// finish drawing
	psfile << "showpage\n";
	psfile.close();
}

double pt_dist(MyPoint a, MyPoint b)
{
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

bool on_obstacle(MyPoint mp)
{
	for (unsigned int i = 0; i < (int)Obstacle.size(); i++)
	{
		if (test_within(mp, (int)Obstacle[i].size(), Obstacle[i]))
		{
			return true;
		}
	}
	return false;
}

//=============================================================================
int main(int argc, char *argv[])
{
	double delta;
	struct timeval gettime[3];

	bool stop_deploy = false, o_free = true, deploy_ok = true;
	MyPoint mp;
	string filename_in = "input"; // default input filename
	fstream inFile;

	if (argc > 1)
	{
		filename_in = argv[1];
	}
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

	string filename_out = filename_in.append("_result_DT-Score");

	/************************** test parameters ****************************/
	int boundary_corner_num;	// # of the sensing area's corners
	double dist;				// (UNUSED) distance of grid points
	double SRange;				// 100% sensing range
	double PRange;				// probable sensing range
	double alpha, beta;			// parameters of probabilistic sensor detection model
	double threshold;			// (UNUSED) miss detection threshold
	int obstacle_num;			// # of obstacles
	int corner_num;				// # of the obstacle's corners
	bool output_coverage;		// display coverage?
	bool output_graph;			// output ps coverage graph?
	unsigned int limit_num;		// maximal # of deployable sensors
	unsigned int candidate_num; // # of candidates at refined deployment step

	// timer start
	gettimeofday(&gettime[0], NULL);

	// read test file =============================================================
	inFile >> boundary_corner_num; // info. of sensing area
	for (int i = 0; i < boundary_corner_num; i++)
	{
		inFile >> mp.x;
		inFile >> mp.y;
		Boundary.push_back(mp);
	}

	inFile >> dist; // (UNUSED) distance of grid points
	inFile >> SRange;
	inFile >> PRange;
	inFile >> alpha;
	inFile >> beta;
	inFile >> threshold; // (UNUSED) miss detection threshold
	inFile >> limit_num;

	inFile >> obstacle_num; // # of obstacles
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

	inFile >> candidate_num; // # of candidates

	inFile >> init_sensor_num; // pre-deployed sensors
	for (unsigned int i = 0; i < init_sensor_num; i++)
	{
		inFile >> mp.x;
		inFile >> mp.y;
		Sensor.push_back(mp);
	}

	inFile >> output_coverage; // display coverage?
	inFile >> output_graph;	// output ps coverage graph?

	// create grid points ===============================================================
	int X_low = (int)Boundary[0].x;
	int X_up = (int)Boundary[2].x;
	int Y_low = (int)Boundary[0].y;
	int Y_up = (int)Boundary[2].y;

	double BOUNDARY_GRID_OFFSET = (SRange - PRange) / sqrt(2.0);
	double GRID_OFFSET = BOUNDARY_GRID_OFFSET;
	double GRID_SENSOR_STEP = GRID_OFFSET;
	double GRID_SENSOR_DIST = (SRange - PRange);

	// add grid points cover the obstacles
	MyPoint a, b, tmp, mp1, mp2;
	double slope, c, lower_bound_x, upper_bound_x, lower_bound_y, upper_bound_y;
	bool mp1_ok, mp2_ok, slope_ok;

	for (unsigned int i = 0; i < (unsigned int)obstacle_num; i++)
	{
		vector<MyPoint> grid_points;

		// add grid points along the edges of obstacle
		for (unsigned int j = 0; j < (unsigned int)Obstacle[i].size(); j++)
		{
			// get the end points of each edge
			a.x = Obstacle[i][j].x;
			a.y = Obstacle[i][j].y;

			if (j + 1 == Obstacle[i].size())
			{ // return to the starting point
				b.x = Obstacle[i][0].x;
				b.y = Obstacle[i][0].y;
			}
			else
			{
				b.x = Obstacle[i][j + 1].x;
				b.y = Obstacle[i][j + 1].y;
			}

			if (a.x != b.x)
			{ // horizontal boundary
				slope_ok = true;

				if (a.x > b.x)
				{ // swap a&b if a.x is larger
					tmp.x = a.x;
					tmp.y = a.y;
					a.x = b.x;
					a.y = b.y;
					b.x = tmp.x;
					b.y = tmp.y;
				}

				slope = (b.y - a.y) / (b.x - a.x); // calculate the slope of edge(ab)
				c = a.y - slope * a.x;			   // get the coefficient of line equation
			}
			else
				slope_ok = false; // vertical edge

			// Type 1: treated as horizontal line (how to select a proper value?)
			if (slope_ok && fabs(slope) <= 1)
			{

				// decide the range of x-axis
				if (a.x <= X_low + BOUNDARY_GRID_OFFSET * 3)
					lower_bound_x = X_low + BOUNDARY_GRID_OFFSET * 3; // Heuristic: deal with edge closes to boundary
				else
					lower_bound_x = a.x;

				if (b.x >= X_up - BOUNDARY_GRID_OFFSET * 3)
					upper_bound_x = X_up - BOUNDARY_GRID_OFFSET * 3;
				else
					upper_bound_x = b.x;

				// decide the range of y-axis
				if (a.y >= b.y)
				{
					upper_bound_y = a.y + GRID_OFFSET;
					lower_bound_y = b.y - GRID_OFFSET;
				}
				else
				{
					upper_bound_y = b.y + GRID_OFFSET;
					lower_bound_y = a.y - GRID_OFFSET;
				}
				// Heuristic: deal with edge closes to boundary
				if (upper_bound_y >= Y_up - BOUNDARY_GRID_OFFSET * 3)
					upper_bound_y = Y_up - BOUNDARY_GRID_OFFSET * 3;
				if (lower_bound_y <= Y_low + BOUNDARY_GRID_OFFSET * 3)
					lower_bound_y = Y_low + BOUNDARY_GRID_OFFSET * 3;

				// only place points outside of obstacle
				mp1_ok = mp2_ok = false;
				mp.x = (a.x + b.x) / 2;
				mp.y = (a.y + b.y) / 2 + 0.1;
				if (test_within(mp, (unsigned int)Obstacle[i].size(), Obstacle[i]))
					mp2_ok = true;
				mp.y = (a.y + b.y) / 2 - 0.1;
				if (test_within(mp, (unsigned int)Obstacle[i].size(), Obstacle[i]))
					mp1_ok = true;

				for (mp.x = lower_bound_x; mp.x <= upper_bound_x; mp.x = mp.x + GRID_SENSOR_STEP)
				{
					mp.y = slope * mp.x + c;

					mp1.x = mp2.x = mp.x;
					mp1.y = mp.y + GRID_OFFSET; // new point with distance GRID_OFFSET to edge
					mp2.y = mp.y - GRID_OFFSET;

					if ((mp1.y >= lower_bound_y) && (mp1.y <= upper_bound_y) && mp1_ok)
					{
						// add mp1
						grid_points.push_back(mp1);
					}
					if ((mp2.y >= lower_bound_y) && (mp2.y <= upper_bound_y) && mp2_ok)
					{
						// add mp2
						grid_points.push_back(mp2);
					}
				}
			}
			// Type 2:  treated as vertical line (fabs(slope) >= 10.0)
			if (slope_ok && fabs(slope) > 1 && fabs(slope) < 10.0)
			{

				// decide the range of x-axis
				if (a.x <= X_low + BOUNDARY_GRID_OFFSET * 3)
					lower_bound_x = X_low + BOUNDARY_GRID_OFFSET * 3; // Heuristic: deal with edge closes to boundary
				else
					lower_bound_x = a.x - GRID_OFFSET;

				if (b.x >= X_up - BOUNDARY_GRID_OFFSET * 3)
					upper_bound_x = X_up - BOUNDARY_GRID_OFFSET * 3;
				else
					upper_bound_x = b.x + GRID_OFFSET;

				// decide the range of y-axis
				if (a.y >= b.y)
				{
					upper_bound_y = a.y;
					lower_bound_y = b.y;
				}
				else
				{
					upper_bound_y = b.y;
					lower_bound_y = a.y;
				}
				// Heuristic: deal with edge closes to boundary
				if (upper_bound_y >= Y_up - BOUNDARY_GRID_OFFSET * 3)
					upper_bound_y = Y_up - BOUNDARY_GRID_OFFSET * 3;
				if (lower_bound_y <= Y_low + BOUNDARY_GRID_OFFSET * 3)
					lower_bound_y = Y_low + BOUNDARY_GRID_OFFSET * 3;

				mp1_ok = mp2_ok = false;
				mp.y = (a.y + b.y) / 2;
				mp.x = (a.x + b.x) / 2 + 0.1;
				if (test_within(mp, (unsigned int)Obstacle[i].size(), Obstacle[i]))
					mp2_ok = true;
				mp.x = (a.x + b.x) / 2 - 0.1;
				if (test_within(mp, (unsigned int)Obstacle[i].size(), Obstacle[i]))
					mp1_ok = true;

				for (mp.y = lower_bound_y; mp.y <= upper_bound_y; mp.y = mp.y + GRID_SENSOR_STEP)
				{
					mp.x = (mp.y - c) / slope;

					mp1.x = mp.x + GRID_OFFSET;
					mp2.x = mp.x - GRID_OFFSET;
					mp1.y = mp2.y = mp.y;

					if ((mp1.x >= lower_bound_x) && (mp1.x <= upper_bound_x) && mp1_ok)
					{
						// add mp1
						grid_points.push_back(mp1);
					}
					if ((mp2.x >= lower_bound_x) && (mp2.x <= upper_bound_x) && mp2_ok)
					{
						// add mp2
						grid_points.push_back(mp2);
					}
				}
			}
			// Type 3:  vertical line
			if (!slope_ok || fabs(slope) >= 10.0)
			{

				// decide the range of y-axis
				if (a.y >= b.y)
				{
					upper_bound_y = a.y;
					lower_bound_y = b.y;
				}
				else
				{
					upper_bound_y = b.y;
					lower_bound_y = a.y;
				}
				// Heuristic: deal with edge closes to boundary
				if (upper_bound_y >= Y_up - BOUNDARY_GRID_OFFSET * 3)
					upper_bound_y = Y_up - BOUNDARY_GRID_OFFSET * 3;
				if (lower_bound_y <= Y_low + BOUNDARY_GRID_OFFSET * 3)
					lower_bound_y = Y_low + BOUNDARY_GRID_OFFSET * 3;

				mp1_ok = mp2_ok = false;
				mp.y = (a.y + b.y) / 2;
				mp.x = (a.x + b.x) / 2 + 0.1;
				if (test_within(mp, (unsigned int)Obstacle[i].size(), Obstacle[i]))
					mp2_ok = true;
				mp.x = (a.x + b.x) / 2 - 0.1;
				if (test_within(mp, (unsigned int)Obstacle[i].size(), Obstacle[i]))
					mp1_ok = true;

				mp.x = a.x;

				for (mp.y = lower_bound_y; mp.y <= upper_bound_y; mp.y = mp.y + GRID_SENSOR_STEP)
				{
					mp1.x = mp.x + GRID_OFFSET;
					mp2.x = mp.x - GRID_OFFSET;
					mp1.y = mp2.y = mp.y;

					if ((mp1.x > X_low + BOUNDARY_GRID_OFFSET) && (mp1.x < X_up - BOUNDARY_GRID_OFFSET) && mp1_ok)
					{
						// add mp1
						grid_points.push_back(mp1);
					}
					if ((mp2.x > X_low + BOUNDARY_GRID_OFFSET) && (mp2.x < X_up - BOUNDARY_GRID_OFFSET) && mp2_ok)
					{
						// add mp2
						grid_points.push_back(mp2);
					}
				}
			} //else
		}	 //for

		// test if the grid points can be added to Sensor
		for (unsigned int k = 0; k < (unsigned int)grid_points.size(); k++)
		{
			mp.x = grid_points[k].x;
			mp.y = grid_points[k].y;

			deploy_ok = true;
			for (unsigned int n = 0; n < (int)Sensor.size() && deploy_ok; n++)
			{
				if (pt_dist(mp, Sensor[n]) < GRID_SENSOR_DIST)
				{ // Heuristic: add sensor if dist. to other sensors >= GRID_SENSOR_DIST
					deploy_ok = false;
				}
			}
			if ((unsigned int)Sensor.size() >= limit_num)
				deploy_ok = false;
			if (deploy_ok)
			{
				if (!on_obstacle(mp))
				{
					Sensor.push_back(mp);
					grid_sensor_num++;
				}
			}
		}
	} //for(unsigned int i=0; i<(int)Obstacle.size(); i++)

	// add grid points along the boundary
	double upper_bound, lower_bound;

	for (unsigned int i = 0; i < (unsigned int)boundary_corner_num; i++)
	{
		vector<MyPoint> grid_points;

		a.x = Boundary[i].x;
		a.y = Boundary[i].y;

		if (i + 1 == boundary_corner_num)
		{
			b.x = Boundary[0].x;
			b.y = Boundary[0].y;
		}
		else
		{
			b.x = Boundary[i + 1].x;
			b.y = Boundary[i + 1].y;
		}

		if (a.x != b.x)
		{ // horizontal boundary
			if (a.x > b.x)
			{ // swap a&b if a.x is larger
				tmp.x = a.x;
				tmp.y = a.y;
				a.x = b.x;
				a.y = b.y;
				b.x = tmp.x;
				b.y = tmp.y;
			}

			lower_bound = a.x + BOUNDARY_GRID_OFFSET; // starting at the left
			upper_bound = b.x - BOUNDARY_GRID_OFFSET; // ending at the right

			mp.y = a.y;

			for (mp.x = lower_bound; mp.x <= upper_bound; mp.x = mp.x + BOUNDARY_GRID_OFFSET * 2)
			{ // each point with step = BOUNDARY_GRID_OFFSET*2
				if (mp.y == Y_low)
				{ // add points cover lower horizontal boundary
					// add mp1 only
					mp1.x = mp.x;
					mp1.y = mp.y + BOUNDARY_GRID_OFFSET;
					grid_points.push_back(mp1);
					if (mp.x + BOUNDARY_GRID_OFFSET * 2 > upper_bound)
					{ // Heuristic: add point to increase coverage
						mp1.x = upper_bound;
						grid_points.push_back(mp1);
					}
				}
				if (mp.y == Y_up)
				{
					// add mp2 only
					mp2.x = mp.x;
					mp2.y = mp.y - BOUNDARY_GRID_OFFSET;
					grid_points.push_back(mp2);
					if (mp.x + BOUNDARY_GRID_OFFSET * 2 > upper_bound)
					{
						mp2.x = upper_bound;
						grid_points.push_back(mp2);
					}
				}
			}
		}
		else
		{ // a.x == b.x
			if (a.y > b.y)
			{
				tmp.x = a.x;
				tmp.y = a.y;
				a.x = b.x;
				a.y = b.y;
				b.x = tmp.x;
				b.y = tmp.y;
			}

			lower_bound = a.y + BOUNDARY_GRID_OFFSET;
			upper_bound = b.y - BOUNDARY_GRID_OFFSET;

			mp.x = a.x;

			for (mp.y = lower_bound; mp.y <= upper_bound; mp.y = mp.y + BOUNDARY_GRID_OFFSET * 2)
			{
				if (mp.x == X_low)
				{
					// add mp1 only
					mp1.x = mp.x + BOUNDARY_GRID_OFFSET;
					mp1.y = mp.y;
					grid_points.push_back(mp1);
					if (mp.y + BOUNDARY_GRID_OFFSET * 2 > upper_bound)
					{
						mp1.y = upper_bound;
						grid_points.push_back(mp1);
					}
				}
				if (mp.x == X_up)
				{
					// add mp2 only
					mp2.x = mp.x - BOUNDARY_GRID_OFFSET;
					mp2.y = mp.y;
					grid_points.push_back(mp2);
					if (mp.y + BOUNDARY_GRID_OFFSET * 2 > upper_bound)
					{
						mp2.y = upper_bound;
						grid_points.push_back(mp2);
					}
				}
			}
		}
		// Heuristic: test if the grid points can be added to Sensor
		for (unsigned int k = 0; k < (unsigned int)grid_points.size(); k++)
		{
			mp.x = grid_points[k].x;
			mp.y = grid_points[k].y;

			deploy_ok = true;

			for (unsigned int n = 0; n < (int)Sensor.size() && deploy_ok; n++)
			{
				if (pt_dist(mp, Sensor[n]) < PRange)
				{ // Heuristic: add sensor if dist. to other sensors >= PRange
					deploy_ok = false;
				}
			}
			if ((unsigned int)Sensor.size() >= limit_num)
				deploy_ok = false;
			if (deploy_ok)
			{
				//if(!on_obstacle(mp) && !on_semi_obstacle(mp)){
				if (!on_obstacle(mp))
				{
					Sensor.push_back(mp);
					grid_sensor_num++;
				}
			}
		}
	} //for(unsigned int i=0; i < (unsigned int)boundary_corner_num; i++)

	// draw initial graph
	if (output_graph)
		draw_graph(filename_in.append("_init.ps"));

	// get init time =======================================================================
	gettimeofday(&gettime[1], NULL);

	// calculate the seconds part of the delta between count and count + 1
	delta = (double)(gettime[1].tv_sec - gettime[0].tv_sec) * (double)MICRO_PER_SECOND;

	// calculate the microseconds part of the delta between count and count + 1
	delta += (double)(gettime[1].tv_usec - gettime[0].tv_usec);

	printf("\nInit time: %.3f seconds.\n", delta / (double)MICRO_PER_SECOND);

	// start ==========================================================================
	bool keep_deploy = true, add_pt = true;
	int max_node; // index for the node with max radius
	double max_radius;

	while ((unsigned int)Sensor.size() < limit_num && keep_deploy == true)
	{
		Delaunay_tree DT;

		// add current sensor points
		for (unsigned int i = 0; i < (unsigned int)Sensor.size(); i++)
			DT += new point(Sensor[i].x, Sensor[i].y);

		// find empty circle from Delaunay triangle
		DT.find_nodes();

		// copy radius vector
		vector<double> Radius_tmp;
		for (unsigned int i = 0; i < (unsigned int)Node.size(); i++)
		{

			Radius_tmp.push_back(Node[i].radius);
		}

		// choose cadidate_num of index, if not enough, choose all
		vector<int> Candidate_index;
		unsigned int dead_num = 0;

		for (unsigned int k = 0; (unsigned int)Candidate_index.size() < candidate_num && (unsigned int)Candidate_index.size() < ((unsigned int)Node.size() - dead_num); k++)
		{
			max_radius = 0.0; // max radius for current availiable nodes

			for (unsigned int i = 0; i < (unsigned int)Radius_tmp.size(); i++)
			{
				if (Radius_tmp[i] > max_radius)
				{
					max_radius = Radius_tmp[i];
					max_node = i;
				}
			}

			if (Node[max_node].center.x >= 0 && Node[max_node].center.x <= X_up && Node[max_node].center.y >= 0 && Node[max_node].center.y <= Y_up)
			{
				deploy_ok = true;

				if (on_obstacle(Node[max_node].center))
				{
					deploy_ok = false;
					Radius_tmp[max_node] = 0.0;
				}

				if (deploy_ok)
				{ // not on obstacles
					{
						Radius_tmp[max_node] = 0.0;
						Candidate_index.push_back(max_node); // add to candidate
					}
				}
				else
				{
					Radius_tmp[max_node] = 0.0; // on obstacles
					dead_num++;
				}
			}
			else
			{
				Radius_tmp[max_node] = -1.0; // outside of sensing area
				dead_num++;
			}
		}

		// create scoring vector
		if (stop_deploy)
			break;

		vector<double> Score;
		double score, point_dist, prob, overlap, effective_radius;
		int startx, endx, starty, endy;
		MyPoint m, n;

		if (Candidate_index.size() != 0)
		{
			// scoring
			for (int k = 0; k < (int)Candidate_index.size(); k++)
			{
				score = 0;

				mp.x = Node[Candidate_index[k]].center.x;
				mp.y = Node[Candidate_index[k]].center.y;
				startx = (int)mp.x - (int)(SRange + PRange);
				endx = (int)mp.x + (int)(SRange + PRange);
				starty = (int)mp.y - (int)(SRange + PRange);
				endy = (int)mp.y + (int)(SRange + PRange);

				if (startx < X_low)
					startx = X_low;
				if (endx > X_up)
					endx = X_up;
				if (starty < Y_low)
					starty = Y_low;
				if (endy > Y_up)
					endy = Y_up;

				for (int i = startx; i <= endx; i = i + (int)PRange) //trace all grids
				{
					for (int j = starty; j <= endy; j = j + (int)PRange)
					{
						MyPoint grid = {(double)i, (double)j};

						point_dist = pt_dist(mp, grid);

						if (point_dist >= (SRange + PRange))
							prob = 0.0;
						else
						{
							// check if blocked by obstacles
							o_free = true;
							for (int r = 0; r < obstacle_num && o_free; r++)
							{
								for (int t = 0; t < (int)Obstacle[r].size() && o_free; t++)
								{
									m = Obstacle[r][t];
									n = Obstacle[r][(t + 1) % (int)Obstacle[r].size()];

									if (intersect(grid, Node[Candidate_index[k]].center, m, n) == true)
										o_free = false;
								}
							}
							if (o_free != true)
								prob = 0.0;
							else // score if not blocked
							{
								if (point_dist <= (SRange - PRange))
									prob = 1.0;
								else
									prob = exp(-alpha * pow((point_dist - (SRange - PRange)) / (2 * PRange), beta));
								score = score + prob;
							}
						}
					} //for j
				}	 //for i
				Score.push_back(score);
			} //for

			for (int i = 0; i < (int)Score.size(); i++)
			{
				overlap = 1.0;
				effective_radius = Node[Candidate_index[i]].radius - (SRange - PRange);

				if (effective_radius < (SRange - PRange))
				{
					overlap = pow(effective_radius / (SRange - PRange), 2); // ratio of non-overlapped sensing region
				}
				Score[i] = overlap * Score[i];
			}

			// add new sensor node
			double max_score = 0.0;
			int win = 0;
			for (int i = 0; i < (int)Score.size(); i++) // find high score
			{
				if (Score[i] > max_score)
				{
					max_score = Score[i];
					win = i;
				}
			}

			Sensor.push_back(Node[Candidate_index[win]].center);

			if ((unsigned int)Node.size() != 0)
				Node.clear();
			if ((unsigned int)Radius_tmp.size() != 0)
				Radius_tmp.clear();
			if ((int)Candidate_index.size() != 0)
				Candidate_index.clear();
			if ((int)Score.size() != 0)
				Score.clear();

		} //if Candidate != 0
	}	 //end of while

	// get execution time =======================================================================
	gettimeofday(&gettime[2], NULL);

	// Calculate the seconds part of the delta between count and count + 1
	delta = (double)(gettime[2].tv_sec - gettime[1].tv_sec) * (double)MICRO_PER_SECOND;

	// Calculate the microseconds part of the delta between count and count + 1
	delta += (double)(gettime[2].tv_usec - gettime[1].tv_usec);

	printf("\nExecution time: %.3f seconds.\n", delta / (double)MICRO_PER_SECOND);

	// output result file =========================================================================
	fstream outFile;
	outFile.open(filename_out.c_str(), ios_base::out | ios_base::in | ios_base::trunc);

	outFile << boundary_corner_num;
	for (int i = 0; i < boundary_corner_num; i++)
	{
		outFile << " " << Boundary[i].x << " " << Boundary[i].y;
	}
	outFile << endl;

	outFile << SRange << endl;
	outFile << PRange << endl;
	outFile << alpha << endl;
	outFile << beta << endl;

	outFile << obstacle_num << endl;
	for (int i = 0; i < obstacle_num; i++)
	{
		outFile << (int)Obstacle[i].size() << " ";
		for (int j = 0; j < (int)Obstacle[i].size(); j++)
		{
			outFile << Obstacle[i][j].x << " " << Obstacle[i][j].y << " ";
		}
		outFile << endl;
	}

	outFile << output_coverage << endl;
	outFile << output_graph << endl;

	outFile << init_sensor_num << endl;

	for (unsigned int i = 0; i < (unsigned int)Sensor.size(); i++)
		outFile << Sensor[i].x << " " << Sensor[i].y << endl;
	outFile.close();

	// output deployment graph
	if (output_graph)
		draw_graph(filename_out.append("_deploy.ps"));

	cout << endl
		 << "init sensors = " << init_sensor_num << endl;
	cout << "grid sensors = " << grid_sensor_num << endl;
	cout << "new sensors = " << (unsigned int)Sensor.size() - init_sensor_num - grid_sensor_num << endl;
}