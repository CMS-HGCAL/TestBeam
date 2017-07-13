#ifndef HGCalTBWireChamberData_H
#define HGCalTBWireChamberData_H

struct WireChamberData {
  explicit WireChamberData(int _ID, double _x, double _y, double _z): ID(_ID), x(_x), y(_y), z(_z) {};
  WireChamberData(): ID(1), x(0), y(0), z(0) {};
  int ID;
  double x;
  double y;
  double z;
  bool goodMeasurement_X;
  bool goodMeasurement_Y;
  bool goodMeasurement;
};

typedef std::vector<WireChamberData> WireChambers;


#endif