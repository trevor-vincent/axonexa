class Polygon_XY {

private:

vector<Line_XY> lines;
int region;

public:

Polygon_XY

void addVertice(double, double);

//determines the point of intersection between the polygon and a ray
bool inside(double [] , int &);
bool inside(Point &);
bool inside(Polygon_XY &);
bool intersect(Line &, double &);
void getnormal(double [], double []);
int getregion(double []);
double getT2(double []);
double getD(double []);

};

void Polygon_XY::addLine(double xi, double yi, double xf, double yf){

	Line_XY line(xi,yi,xf,yf);
	lines.push_back(line);

}

bool Polygon_XY::inside(double x[], int & reg){

	int polySides = lines.size();
	bool oddNodes = false;
	
	double x = x[0];
	double y = x[1];
	
	for (int i = 0; i < lines.size(); i++) {
		
		if (( lines[i].yi < y && doub_gtoe(lines[i].yf, y) ||   lines[i].yf < y && doub_gtoe(lines[i].yi,y) ) &&  ( doub_stoe(lines[i].xi,x) || doub_stoe(lines[i].xf,x))) {
			oddNodes^=(lines[i].xi +(y-lines[i].yi)/(lines[i].yf-lines[i].yi)*(lines[i].xf-lines[i].xi)<x); 
		}
	}

	if (oddNodes == true){reg = region; return true;}
	return false;
}

bool Polygon_XY::inside(Point & test){

	int polySides = lines.size();
	bool oddNodes = false;
	
	double x = test.x;
	double y = test.y;
	
	for (int i = 0; i < lines.size(); i++) {
		
		if (( lines[i].yi < y && doub_gtoe(lines[i].yf, y) ||   lines[i].yf < y && doub_gtoe(lines[i].yi,y) ) &&  ( doub_stoe(lines[i].xi,x) || doub_stoe(lines[i].xf,x))) {
			oddNodes^=(lines[i].xi +(y-lines[i].yi)/(lines[i].yf-lines[i].yi)*(lines[i].xf-lines[i].xi)<x); 
		}
	}

	return oddNodes;


}

bool Polygon_XY::inside(Polygon_XY test){

	for (int i = 0; i < test.size(); i++){
	
		Point p1(test.lines[i].xi, test.lines[i].yi);
		Point p2(test.lines[i].xf, test.lines[i].yf);
		
		if (!inside(p1)) { return false;}
		if (!inside(p2)) { return false;}
	
	}
	
	return true;
}

bool Polygon_XY::intersect(Line & line, double & v){

	double v_temp;
	vector<double> pos_v;

	for (int i = 0; i < lines.size(); i++){
		
		if (line.intersect(lines[i], v_temp) ){
			pos_v.push_back(v_temp);
		}
	}
	
	if (pos_v.isempty()){
		return false;
	}
	
	else {
		std::sort(pos_v.begin(), pos_v.end());
		v = pos_v[0];
	}
	
}

