	// /************************************************************
	// *
	// *	RANDOMLY PACKED SQUARE LATTICE
	// *
	// *   This is only for objects with a radius.
	// *   ls = lattice_size, num = number of randomly placed objects
	// *
	// *
	// *************************************************************/
	
	class RPSLatticeInitializer {
	
		public:
	
			template <class OBJECT>
			RPSLattice(int num, Lattice<OBJECT> & lattice){
				lattice.basis.resize(num);
				lattice.setLatticeVectors(Vector3(1.0,0.0,0.0),Vector3(0.0,1.0,0.0),Vector3(0.0,0.0,1.0) );
			}
			
			~RPSLattice(){}
			
			template <class OBJECT>
			void uniformRadialDist(Lattice<OBJECT> & lattice, int seed, double max_radius, double min_radius){
			
					Ran randomnumber(seed);
					
					vector<double> radii(basis.size());
					
					double deviate;

					for (int i = 0; i < num_of_cyl_vox; i++) {
						deviate = (max_radius-min_radius)*randomnumber.doub() + min_radius;
						if (deviate > 0.0) {radii[i] = deviate;}
						else {i--;}
					}
					
					sort(radii.begin(), radii.end());
					
					for (int i = 0; i < basis.size(); i++){
						basis[i].setRadius(radii[i]);
					}					
			
			
			}
			
			template <class OBJECT>
			void gaussianRadialDist(Lattice<OBJECT> & lattice, int seed, double mean, double stddev){
			
					Normaldev nd(mean, stddev, seed);
					
					vector<double> radii(basis.size());
					
					double deviate;

					for (int i = 0; i < basis.size(); i++) {
						deviate = nd.dev();
						if (deviate > 0.0) {radii[i] = deviate;}
						else {i--;}
					}		
					
					sort(radii.begin(), radii.end());
					
					for (int i = 0; i < basis.size(); i++){
						basis[i].setRadius(radii[i]);
					}
			
			}
			
			template <class OBJECT>
			void gammaRadialDist(Lattice<OBJECT> & lattice, int seed, double alpha, double beta){
			
					Gammadev gd(alpha, beta, seed);
					
					vector<double> radii(basis.size());
					
					double deviate;

					for (int i = 0; i < num_of_cyl_vox; i++) {
						deviate = gd.dev();
						if (deviate > 0.0) {radii[i] = deviate;}
						else {i--;}
					}
					
					sort(radii.begin(), radii.end());
					
					for (int i = 0; i < basis.size(); i++){
						basis[i].setRadius(radii[i]);
					}			
			
			}
			
			template <class OBJECT>
			void uniformCenterDist(Lattice<OBJECT> & lattice, int seed){
			
				Ran randomnumber(seed);
				vector<Point> center(basis.size());
				
				for (int i = basis.size()-1; i > -1; i--){
					
					center[i].x = randomnumber.doub()*a;
					center[i].y = randomnumber.doub()*b;
					center[i].z = randomnumber.doub()*c;
					basis[i].setCenter(center[i]);
					
					while( basisOverlap(basis[i]) ){
						
						center[i].x = randomnumber.doub()*a;
						center[i].y = randomnumber.doub()*b;
						center[i].z = randomnumber.doub()*c;

						basis[i].setCenter(center[i]);

					}
					
				}
			
			}
			
			template <class OBJECT>
			void correctEdges(Lattice<OBJECT> & lattice){
			
				Object temp;
				Point ptemp;
				for (int ind = 0; ind < basis.size(); ind++){
					
					if (objectOutside(basis[ind])){
					
						for (int i = -1; i <= 1; i++){
						for (int j = -1; j <= 1; j++){
						for (int k = -1; k <= 1; k++){
							temp = basis[ind]; 
							ptemp = basis[ind].getCenter(); 
							Point cnew = Point( ptemp.x + a*i, ptemp.y + b*j, ptemp.z + c*k);
							temp.setCenter(cnew);
							basis.push_back(temp);
						}}}

					}
					
				}
			
			}
			
			template <class OBJECT>
			bool basisOverlap(OBJECT & o){
			
				for (int i = 0; i < basis.size(); i++){
				
					if (o.overlap(basis[i]){
					
						return true;
					
					}
				
				}
				return false;
			}
			
			template <class OBJECT>
			bool objectOutside(OBJECT & o){
			
				Vector3 xhat(1.0,0.0,0.0);
				Vector3 yhat(0.0,1.0,0.0);
				Vector3 zhat(0.0,0.0,1.0);
				Point p = o.getCenter();
				Vector3 center( p.x, p.y, p.z );
				
				if ( !inLattice(center + xhat*o.radius) ) {return false;}
				else if ( !inLattice(center - xhat*o.radius) ) {return false;}
				else if ( !inLattice(center + yhat*o.radius) ) {return false;}
				else if ( !inLattice(center - yhat*o.radius) ) {return false;}
				else if ( !inLattice(center + zhat*o.radius) ) {return false;}
				else if ( !inLattice(center - zhat*o.radius) ) {return false;}
				else {return true;}
			
			}
	
	};