#ifndef HGCalTBMultiWireChamberData_H
#define HGCalTBMultiWireChamberData_H

struct MultiWireChamberData {
  explicit MultiWireChamberData(int _ID, double _x, double _y, double _z): ID(_ID), x(_x), y(_y), z(_z) {};
  MultiWireChamberData(): ID(1), x(0), y(0), z(0) {};
  int ID;
  double x;
  double y;
  double z;
};

typedef std::vector<MultiWireChamberData> MultiWireChambers;


#endif