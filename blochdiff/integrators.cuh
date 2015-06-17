template<class Basis, class Fieldtype>
class rotate_NoRelaxation{

	__device__ __host__ void integrate(
								real* Mx, 
								real* My, 
								real* Mz, 
								const Basis *basis,
								const FieldType & Bfield, 
								const Vector3 & r, 
								int i, 
								real h
							  )
				  
	{
		Vector3 B = Bfield(r,i*h);
		Vector3 M(*Mx, *My, *Mz);
		real Bmag = B.magnitude();
		Vector3 C = (B*(1.0/Bmag))*(M*(B*(1.0/Bmag)));
		Vector3 u = M-C;
		Vector3 finalM =  (C + u*cos(-GAMMA*Bmag*timestep) + ((B*(1.0/Bmag)) % u)*sin(-GAMMA*Bmag*timestep));
		*Mx = finalM.x;
		*My = finalM.y;
		*Mz = finalM.z;		
	}


};


class rotate_NoRelaxation_NoTimestep{

void integrate(


)


	void integrate(
					real* Mx, 
					real* My, 
					real* Mz, 
					const Vector3 & B, 
					real timestep
				  )
				  
	{
		Vector3 M(*Mx, *My, *Mz);
		real Bmag = B.magnitude();
		Vector3 C = (B*(1.0/Bmag))*(M*(B*(1.0/Bmag)));
		Vector3 u = M-C;
		Vector3 finalM =  (C + u*cos(-GAMMA*Bmag*timestep) + ((B*(1.0/Bmag)) % u)*sin(-GAMMA*Bmag*timestep));
		*Mx = finalM.x;
		*My = finalM.y;
		*Mz = finalM.z;		
	}


};