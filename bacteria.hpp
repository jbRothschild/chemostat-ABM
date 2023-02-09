#ifndef BACTERIA_H
#define BACTERIA_H
#include <stdio.h>      /* printf */
#include <math.h>       /* sqrt */

#define AREA_BACT 4.427928249198705
#define PI 3.14159265359
#define BACT1_M 4.56
#define BACT2_M 5.56
#define BACT1_G 0.01155  // 1h doubling time
#define BACT2_G 0.01155

double radiusConstantArea(double length)
{
  return (-2 * length + sqrt(4 * pow(length, 2.) + 4 * AREA_BACT * (PI - 4)))
          / ( 2 * ( PI - 4 ));
}

class Bacteria
{
  public:
    double radius;
    double max_length;
    double inertia;
    double growth_rate;
    double length;
    Bacteria(double r = 0.1, double m = 2.0, double i = 5.0, double g = 0.01,
             double l = 1.0){
             radius = r; 
             max_length = m;
             inertia = i;
             growth_rate = g;
             length = l;
    }
};

class EColiGFP: public Bacteria
{
  public:
    EColiGFP(double r = radiusConstantArea(BACT1_M), double m = BACT1_M,  
             double i = 5.0, double g = BACT1_G, double l = 2.0)
                       : Bacteria(r, m, i, g, l){ }
};

class EColiMCh1: public Bacteria
{
  public:
    EColiMCh1(double r = radiusConstantArea(BACT1_M), double m = BACT1_M,  
              double i = 5.0, double g = BACT1_G, double l = 2.0)
                       : Bacteria(r, m, i, g, l){ } // g = 0.009242??
};

class EColiMCh2: public Bacteria
{
  public:
    EColiMCh2(double r = radiusConstantArea(BACT1_M), double m = BACT1_M,  
              double i = 5.0, double g = BACT1_G, double l = 2.0)
                       : Bacteria(r, m, i, g, l){ }
};

class EColiA22: public Bacteria
{
  public:
    EColiA22(double r = radiusConstantArea(BACT2_M), double m = BACT2_M,  
             double i = 5.0, double g = BACT2_G, double l = 2.0)
                       : Bacteria(r, m, i, g, l){ }
};

class BSubt: public Bacteria
{
  public:
    BSubt(double r = 0.83 / 2.0 + 0.1, double m = 7.95,
          double i = 5.0, double g = 0.0039, double l = 3.00)
                       : Bacteria(r, m, i, g, l){ }
    double radius = 0.83 / 2.0 + 0.1;
    double max_length = 7.95;
    double inertia = 5.0;
    double growth_rate = 0.0039;
    double length = 3.00;
};

#endif
