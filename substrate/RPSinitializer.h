// /************************************************************
// *
// *	
// *	    RANDOMLY PACKED SQUARE LATTICE Initializer
// *
// *
// *               ( Basis must have a radius )
// *************************************************************/

template <class Basis>	
class RPSLatticeInitializer {
	
 public:

  std::vector<Basis> basis;
  Lattice lat;
  int basisType; //0 = cylinder, 1 = sphere
  bool verbose;

 public:
	
  
  RPSLatticeInitializer(Lattice & _lat, int _basisType, bool _verbose = false){
    basis = std::vector<Basis>(_lat.getBasisSize());
    lat = _lat;
	basisType = _basisType;
	verbose = _verbose;
  }

  ~RPSLatticeInitializer(){}
			

  void uniformRadialDist(int seed, real max_radius, real min_radius){
			
    Ran randomnumber(seed);
    std::vector<real> radii( basis.size());
    real deviate;

    for (int i = 0; i < basis.size() ; i++) {
      deviate = (max_radius-min_radius)*randomnumber.doub() + min_radius;
      if (deviate > 0.0) {radii[i] = deviate;}
      else {i--;}
    }
					
    std::sort(radii.begin(), radii.end());
					
    for (int i = 0; i < basis.size(); i++){
      basis[i].setRadius(radii[i]);
    }					

  }
			

  void gaussianRadialDist(int seed, double mean, double stddev){
			
    Normaldev nd(mean, stddev, seed);
					
    std::vector<real> radii(basis.size());
					
    double deviate;

    for (int i = 0; i < basis.size(); i++) {
      deviate = nd.dev();
      if (deviate > 0.0) {radii[i] = deviate;}
      else {i--;}
    }
					
    std::sort(radii.begin(), radii.end());
					
    for (int i = 0; i < basis.size(); i++){
      basis[i].setRadius(radii[i]);
    }
			
  }
			

  void gammaRadialDist(int seed, double alpha, double beta, real minimum, real maximum){
			
    Gammadev gd(alpha, beta, seed);
					
    std::vector<real> radii(basis.size());
					
    double deviate;

    for (int i = 0; i < basis.size(); i++) {
      deviate = gd.dev();
      if (deviate > minimum && deviate < maximum) {radii[i] = deviate;}
      else {i--;}
    }
		
    std::sort(radii.begin(), radii.end());
					
    for (int i = 0; i < basis.size(); i++){
      basis[i].setRadius(radii[i]);
	  if (verbose) std::cout << " Radius " << i << " " << radii[i] << std::endl;
    }		
	
	
	
  }
			

  bool uniformCenterDist(int seed, int triesMax){
			
    Ran randomnumber(seed);
	int size = basis.size();
    std::vector<Vector3> center(size);
	int tries = 0;	
	
    for (int i = size-1; i > -1; i--){
	  tries = 0;			
      do {					

		center[i].x = randomnumber.doub()*lat.getA();
		center[i].y = randomnumber.doub()*lat.getB();
		center[i].z = randomnumber.doub()*lat.getC();
		basis[i].setCenter(center[i]);
		if (verbose) std::cout << " Setting center " << i << " = " << center[i] << std::endl;
		if (tries > triesMax*(size-i)){return false;}
		tries++;
      } while( basisOverlap(basis[i] , i ) );
	
    }
	
	return true;
			
  }
  
	
  void correctEdges(){
			
    Basis temp;
    Vector3 ptemp;

    real a = lat.getA();
    real b = lat.getB();
    real c = lat.getC();
	int isize = basis.size();

    for (int ind = 0; ind < isize; ind++){
					
	  basis[ind].setRegion(ind+1);				
      if (objectOutside(basis[ind])){
	
if ( basisType == 1){//spheres	
	for (int i = -1; i <= 1; i++){
	  for (int j = -1; j <= 1; j++){
	    for (int k = -1; k <= 1; k++){
		  if (!(i == 0 && j == 0 && k ==0)){
	      temp = basis[ind]; 
	      ptemp = basis[ind].getCenter(); 
	      Vector3 cnew = Vector3( ptemp.x + a*i, ptemp.y + b*j, ptemp.z + c*k);
	      temp.setCenter(cnew);
	      basis.push_back(temp);
		  }
	    }}}
	}
	
	else{
	for (int i = -1; i <= 1; i++){
	  for (int j = -1; j <= 1; j++){
		  if (!(i == 0 && j == 0)){
	      temp = basis[ind]; 
	      ptemp = basis[ind].getCenter(); 
	      Vector3 cnew = Vector3( ptemp.x + a*i, ptemp.y + b*j, 0.0);
		  temp.setCenter(cnew);
		  if(circleIntersectsLattice(temp.getCenter().x, temp.getCenter().y, temp.getRadius())){
		  basis.push_back(temp);
		  }
		  }
	    }}
	
	}
	
      }
					
    }
	
	lat.setBasisSize(basis.size());
			
  }
			
 
 //only works for cylinders atm
 //to make it work for spheres we must add a j,k and l (triple for loop) to check for z overlaps
    bool basisOverlap(Basis & o, int oi){
      real a = lat.getA();
      real b = lat.getB();
      real c = lat.getC();			
	  for (int l = -1; l < 1; l++){
	  for (int m = -1; m < 1; m++){
      Vector3 center1 = o.getCenter() + Vector3(a*l,b*m,0.0);
      real radius1 = o.getRadius();

      for (int i = oi + 1; i < basis.size(); i++){
	
	  for (int j = -1; j < 1; j++){
	  for (int k = -1; k < 1; k++){
	  Vector3 center2 = basis[i].getCenter() + Vector3(a*j,b*k,0.0);
	  real radius2 = basis[i].getRadius();
	  real dist = (center1 - center2).magnitude();
      
	  if (dist < radius1 + radius2){
	
	   return true;
					
	  }
				}}
      }
	  }}
      return false;
     }
			

    bool objectOutside(Basis & o){
			
      Vector3 xhat(1.0,0.0,0.0);
      Vector3 yhat(0.0,1.0,0.0);
      Vector3 zhat(0.0,0.0,1.0);
      Vector3 center = o.getCenter();
      real radius = o.getRadius();
	  Vector3 r1 = 	center + xhat*radius;
	  Vector3 r2 = 	center - xhat*radius;
	  Vector3 r3 = 	center + yhat*radius;
	  Vector3 r4 = 	center - yhat*radius;
	  Vector3 r5 = 	center + zhat*radius;
	  Vector3 r6 = 	center - zhat*radius;
      if ( !lat.inLatticeCell(r1) ) {return true;}
      else if ( !lat.inLatticeCell(r2) ) {return true;}
      else if ( !lat.inLatticeCell(r3) ) {return true;}
      else if ( !lat.inLatticeCell(r4) ) {return true;}
      else if ( !lat.inLatticeCell(r5) && basisType == 1) {return true;}
      else if ( !lat.inLatticeCell(r6) && basisType == 1) {return true;}
      else {return false;}
			
    }
	
	bool circleIntersectsLattice(real x, real y, real r){
		
		real recthalfWidth = lat.getA()/2.0;
		real recthalfHeight = lat.getB()/2.0;
		real dist2rectcen_x = fabs(x - recthalfWidth);
		real dist2rectcen_y = fabs(y - recthalfHeight);

		if (dist2rectcen_x > (recthalfWidth + r)) { return false; }
		if (dist2rectcen_y > (recthalfHeight + r)) { return false; }

		if (dist2rectcen_x <= (recthalfWidth)){ return true; } 
		if (dist2rectcen_y <= (recthalfHeight)){ return true; }

		real cornerDistance_sq = (dist2rectcen_x - recthalfWidth)*(dist2rectcen_x - recthalfWidth) +
							 (dist2rectcen_y - recthalfHeight)*(dist2rectcen_y - recthalfHeight);

		return (cornerDistance_sq <= (r*r));
	}
	
	void setRegions(){
	
		for (int i = 0; i < basis.size(); i++){
			basis[i].setRegion(i);
		}
	
	}
	
	//only works for cylinders atm
	real getVI(){
	
	std::vector<real> radiiSoFar; 
	real tempRadius, VI = 0.;
	bool add;
	
		for (int i = 0; i < basis.size(); i++){
			tempRadius = basis[i].getRadius();
			add = true;
			
			for (int j = 0; j < radiiSoFar.size(); j++){
				if ( real_equal(tempRadius,radiiSoFar[j],EPSILON) ){add = false;}
			}
			
			if(add) radiiSoFar.push_back(tempRadius);
		}
		
			
		for (int j = 0; j < radiiSoFar.size(); j++){
			VI += PI*radiiSoFar[j]*radiiSoFar[j];
		}

	return VI/(lat.getA()*lat.getB());
	
	}
	
	
	
  };
