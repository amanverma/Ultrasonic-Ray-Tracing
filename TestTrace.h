#include <iosfwd>
#include <vector>
#include <fstream>
#include "Def.h"
#include "vector.h"

struct ray{vector3d startPoint; vector3d endPoint;};
struct transducer{vector3d position; double intensity[30] ;};
struct material {double diffuse; double specular; double power;double reflectance;};

