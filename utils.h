#ifndef __UTILS_H__
#define __UTILS_H__

// hk point east is lat
// hk point north is lon

#include <cmath>
#include <iostream>

namespace bus {
namespace updater{

struct Hk80Point
{
    double east_;
    double north_;
};

struct Wgs84Point
{
    Wgs84Point() {}
    Wgs84Point(double lat, double lon):lat_(lat), lon_(lon)  {}
    double lon_;
    double lat_;
};

void convertHK80GridToCartesian (const Hk80Point& hk_point, Wgs84Point& wgs_point);

}
}

#endif //__UTILS_H__
