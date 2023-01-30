#ifndef BACTERIA_H
#define BACTERIA_H

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
    EColiGFP(double r = 0.92 / 2.0 + 0.05, double m = 4.56,  // 0.009242
             double i = 5.0, double g = 0.01155, double l = 2.25)
                       : Bacteria(r, m, i, g, l){ }
};

class EColiMCh1: public Bacteria
{
  public:
    EColiMCh1(double r = 0.92 / 2.0 + 0.05, double m = 4.56,  // g = 0.008557
              double i = 5.0, double g = 0.01155, double l = 2.25)
                       : Bacteria(r, m, i, g, l){ }
};

class EColiMCh2: public Bacteria
{
  public:
    EColiMCh2(double r = 0.92 / 2.0 + 0.05, double m = 4.56,
              double i = 5.0, double g = 0.009242, double l = 2.25)
                       : Bacteria(r, m, i, g, l){ }
};

class EColiA22: public Bacteria
{
  public:
    EColiA22(double r = 1.20 / 2.0 + 0.05, double m = 4.10,
             double i = 5.0, double g = 0.0173, double l = 1.25)
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
