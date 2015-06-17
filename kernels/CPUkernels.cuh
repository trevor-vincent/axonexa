template <class FieldType, class Basis>
void updateWalkersMagCPU(

				 const SimuParams *par, 
				 const Basis* basis, 
				 const FieldType * B, 

				 real* Mx, 
				 real* My, 
				 real* Mz
				 
				 )

{

for (int tid = 0; tid < par->number_of_particles; tid++){


  real phi, theta;
  real speed = sqrt(6.0*basis->getD()/par->timestep);
	
  Vector3 r = basis[0].unifRandCPU();
  
  for (int i = 0; i < par->measurements; i++){

    Mx[tid + i*par->number_of_particles] = par->mx_initial;
    My[tid + i*par->number_of_particles] = par->my_initial;
    Mz[tid + i*par->number_of_particles] = par->mz_initial;
    
	//printf ( " Kernel Minitial = %.6f %.6f %.6f \n",  Mx[tid + i*par->number_of_particles], My[tid + i*par->number_of_particles], Mz[tid + i*par->number_of_particles] );

    updateMag_noT2_rotate(  &Mx[tid + i*par->number_of_particles] ,  
			    &My[tid + i*par->number_of_particles] , 
			    &Mz[tid + i*par->number_of_particles] , 
			    B[i](r,(real)0.0), 
			    par->timestep );
 
  }
	
  for (int i = 1; i < par->steps; i++){
  
    Vector3 ri = r;
	
	
    phi = 2.0*PI*unifRandCPP();
    theta = acos(2.0*unifRandCPP() - 1);
	
	r += Vector3(speed*par->timestep*sin(theta)*cos(phi), speed*par->timestep*sin(theta)*sin(phi), speed*par->timestep*cos(theta));
   
#if defined SPECULAR_REFLECTION

    boundaryNormalCPU(ri,r, speed, basis, par->timestep);

#else

    if ( !basis[0].inside(r) )
      {
		r = ri;
      }    
    
#endif
    for (int j = 0; j < par->measurements; j++){
	  
      updateMag_noT2_rotate(  &Mx[tid + j*par->number_of_particles] ,  
							&My[tid + j*par->number_of_particles] , 
							&Mz[tid + j*par->number_of_particles] ,
							B[j](r,i*par->timestep), 
							par->timestep );
	//printf ( " Kernel Mafter = %.6f %.6f %.6f \n",  Mx[tid + j*par->number_of_particles], My[tid + j*par->number_of_particles], Mz[tid + j*par->number_of_particles] );

    } 
	
  }
  std::cout << " final r = " << r << std::endl; 
}

}

template <class FieldType, class Basis>
void updateWalkersPhaseCPU(

				 const SimuParamsPhase *par, 
				 const Basis* basis, 
				 const FieldType * G, 
				 real* phase
				 
				 )

{

for (int tid = 0; tid < par->number_of_particles; tid++){


  real phi, theta;
  real speed = sqrt(6.0*basis->getD()/par->timestep);
  real stepLength = speed*par->timestep;
  
  Vector3 r = basis[0].unifRandCPU();

  for (int i = 0; i < par->measurements; i++){
	phase[tid + i*par->number_of_particles] = GAMMA*(G[i](0.0)*r)*par->timestep;
  }
	
  for (int i = 1; i < par->steps; i++){
  
    Vector3 ri = r;
	
#if defined LATTICE_PICKING
	
	r += Vector3(stepLength*randomSign(),stepLength*randomSign(),stepLength*randomSign());
	
#else

    phi = 2.0*PI*unifRandCPP();
    theta = acos(2.0*unifRandCPP() - 1);
	
	r += Vector3(speed*par->timestep*sin(theta)*cos(phi), speed*par->timestep*sin(theta)*sin(phi), speed*par->timestep*cos(theta));
#endif 

#if defined SPECULAR_REFLECTION

    boundaryNormalCPU(ri,r, speed, basis, par->timestep);

#else

    if ( !basis[0].inside(r) )
      {
		r = ri;
      }    
    
#endif

    for (int j = 0; j < par->measurements; j++){
	  phase[tid + j*par->number_of_particles] += GAMMA*(G[j](par->timestep*i)*r)*par->timestep;
    } 
	
  }

}

}

template <class Basis>
__host__ void boundaryNormalCPU(Vector3 & ri, Vector3 & r, const real currentSpeed,  const Basis* basis, const real timestep){

	real v = 0.0;
	real accumtime = 0.0;
	
	while (basis->intersect(ri, r, v)){
			
				r = (r-ri)*v + ri;

				//determine incident vector
				const Vector3 I = (r-ri);
				const real Imag = I.magnitude();

				//calculate accumulated time 
				accumtime += Imag/currentSpeed;
    
				//determine normal for reflection calculation
				const Vector3 n = basis->getNormal(r);
    
				//set last position to boundary position before rflection
				ri = r;

				//reflect, reflection travel time is the remaining time (timestep - accumtime)
				r += ( I - n*2.0*(n*I) )*(currentSpeed/Imag)*( timestep - accumtime );	
		}

}
