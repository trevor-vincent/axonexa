//#include "Objects.h"

#ifndef PARTICLE_H_
#define PARTICLE_H_

class Particle {

  public:
    // main particle info
    Vector3 position;
    double speed;
    double T2factor;

    // kept for records
    Vector3 unbounded_position;
    Vector3 initial_position;
    int region;

  public:
    Particle(double x, double y, double z) {

        position = Vector3(x, y, z);
        initial_position = Vector3(x, y, z);
        unbounded_position = Vector3(x, y, z);
        T2factor = 1.0;
    }

    friend std::ostream &operator<<(std::ostream &out, Particle &p) {

        out << "Current Position: " << p.position << std::endl
            << "Unbounded Position: " << p.unbounded_position << std::endl
            << "Speed: " << p.speed << std::endl;
        //<< "Region: " << p.region << std::endl;

        return out;
    }

    Particle() {}
    ~Particle() {}
};

#endif
