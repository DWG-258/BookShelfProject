#ifndef OBJECTS_H
#define OBJECTS_H
#include "common.hpp"

class Module;
class SiteRow;
class Pin;
class Net;
class PlaceData;

class Module
{
public:
    Module()
    {
        Init();
    }
    int idx;
    string name;
    float width;
    float height;
    float area;
    float orientation;
    bool isMacro;
    bool isFixed;
    bool isFiller;
    vector<shared_ptr<Pin>> modulePins;
    vector<shared_ptr<Net>> nets;
    void Init()
    {
        idx = -1;
        center.SetZero();
        width = 0;
        height = 0;
        area = 0;
        orientation = 0;
        isMacro = false;
        isFiller = false;
        isFixed = false;
     }
    POS_2D center; 
 
};
class SiteRow 
{
public:
    SiteRow()
    {
        bottom = 0;
        height = 0;
        step = 0;
        start.SetZero();
        end.SetZero();
    }

    double bottom;             
    double height;             
    double step;                
    POS_2D start;               
    POS_2D end;                 
};

class Pin
{
public:
    Pin()
    {
        init();
    }

    void init()
    {
        idx = -1;
        offset.SetZero();

    }
    int idx;
    weak_ptr<Module> module;
    weak_ptr<Net> net;
    POS_2D offset;
};

class Net {
public:
  Net(string &_name) : name(_name) { init(); }
  string name;
  vector<shared_ptr<Pin>> netPins;

  void init() { netPins.clear(); }
};

#endif