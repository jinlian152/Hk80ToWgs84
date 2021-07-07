#include <cmath>
#include <iostream>
#include "utils.h"

const int g_a = 6378388;
const double g_f = 0.0033670033670034;
const double g_e2 = 0.0067226700223333;
const int g_m0 = 1;

const double g_lat0 = 0.38942018399288;
const double g_lng0 = 1.99279173737270;
const double g_n0 = 819069.80;
const double g_e0 = 836694.05;

double getMeridianDist(double lat) {
    double a0 = 1 - g_e2 / 4 - 3 * g_e2 * g_e2 / 64;
    double a2 = 3 * (g_e2 + g_e2 * g_e2 / 4) / 8;
    double a4 = 15 * g_e2 * g_e2 / 256;
    return g_a * (a0 * lat - a2 * sin(2 * lat) + a4 * sin(4 * lat));
}

double self_abs(double value) {
    return value >= 0 ? value : (0 - value);
}

double doBisectIter(double m, double x1, double x2, double epsilon = 0.00000001) {
    printf("m=%lf x1=%lf    x2=%lf\n", m, x1, x2);
    double y1 = m - getMeridianDist(x1);
    double y2 = m - getMeridianDist(x2);

    printf("y1=%lf    y2=%lf\n", y1, y2);
    if (y1 * y2 > 0) {
        return -__INT_MAX__ * 1.0;
    }

    if (y1 > 0) {
        return doBisectIter(m, x2, x1, epsilon);
    }

    double x = (x1 + x2) / 2;
    double y = m - getMeridianDist(x);

    if (self_abs(y) <= epsilon) {
        return x;
    }

    return doBisectIter(m, (y < 0) ? x : x1, (y > 0 )? x : x2, epsilon);
}

double rad2deg(double x) {
    return x * 180 / M_PI;
}

void convertHK80GridToCartesian (const Hk80Point& hk_point, Wgs84Point& wgs_point) {
    //printf("x1=%lf  x2=%lf\n", hk_point.north_, hk_point.east_);
    double north = hk_point.north_;
    double east = hk_point.east_;
    printf("north=%lf east=%lf\n", north, east);

    double dn = north - g_n0;
    double de = east - g_e0;

    double m0 = getMeridianDist(g_lat0);
    double m = (dn + m0) / g_m0;
    printf("m0=%lf dn=%lf m=%lf\n", m0, dn, m);

    double lat_p = doBisectIter(m, -10000, 10000);
    if ((-__INT_MAX__ * 1.0) == lat_p) {
        wgs_point.lat_ = __INT_MAX__;
        wgs_point.lon_ = __INT_MAX__;
        return;
    }

    double v_p = g_a / pow(1 - g_e2 * pow(sin(lat_p), 2), 0.5);
    double p_p = g_a * (1 - g_e2) / pow(1 - g_e2 * pow(sin(lat_p), 2), 1.5);
    double phi_p = v_p / p_p;

    double lat = lat_p - tan(lat_p) * pow(de / g_m0, 2) / (2 * p_p * v_p);
    double lng = g_lng0 + de / (g_m0 * v_p * cos(lat_p)) - (pow(de, 3) / (6 * pow(g_m0 * v_p, 3) * cos(lat_p))) * (phi_p + 2 * pow(tan(lat_p), 2));
    wgs_point.lat_ = rad2deg(lat) - 5.5 / 3600;
    wgs_point.lon_ = rad2deg(lng) + 8.8 / 3600;
}

// int main() {
//     bus::updater::Hk80Point hk_point;
//     hk_point.east_ = 844202.18750000;
//     hk_point.north_ = 819959.56250000;
//     bus::updater::Wgs84Point wgs_point;
//     convertHK80GridToCartesian(hk_point, wgs_point);
//     std::cout << "lat=" << wgs_point.lat_ << "   lng=" << wgs_point.lon_ << std::endl;
// }
