/*H**********************************************************************
* FILENAME :        bacteria.hpp
*
* DESCRIPTION :
*       File where different bacterial species/strains are defined. You can
*       define any new bacteria strains here and add them to different
*       initialize functions in des.cpp
*
* PUBLIC FUNCTIONS :
*       double    radiusConstantArea( length )
*       class     Bacteria( radius,
*                           max_length,
*                           inertia,
*                           growth_rate,
*                           length,
*                           damping )
*       class     EColiGFP( Bacteria )
*       class     EColiMCh1( Bacteria )
*       class     EColiMCh2( Bacteria )
*       class     EColiA22( Bacteria )
*       class     BSubt( Bacteria )
*
* NOTES :
*       This whole code works purely for bacteria that are spherocylinders.
*       For now, the only characteristics of a bacterial strain are radius,
*       max length, inertia, growth rate, starting length and damping. You
*       can define other characteristics, but recognize that they have no
*       meaning in the simulation code until you add it yourself 
*
*
* AUTHOR :    Jeremy Rothschild        START DATE :    06-2021
*
* CHANGES :
*
*
*H*/

#ifndef BACTERIA_H
#define BACTERIA_H
#include <stdio.h>      /* printf */
#include <math.h>       /* sqrt */

#define AREA_BACT 4.427928249198705 // area of a bacteria. for the data worked
                                    // with Tianyi, we assumed the bacteria had
                                    // all the same maximal area
#define PI 3.14159265359
#define BACT1_M 4.56 // bacteria maximal length 1
#define BACT2_M 4.16 // bacteria maximal length 2
#define BACT1_G 0.012836// bacteria growth rate 1: 54min double time, ln(2) / 54
#define BACT2_G 0.01155 // bacteria growth rate 2: 1h doubling time, ln(2) / 60

double radiusConstantArea(double length)
{
  // assuming the bacteria all have a similar maximum area (rougly speaking)
  // we've used this function to calculate the appropriate radius they should
  // have for different lengths
  return (-2 * length + sqrt(4 * pow(length, 2.) + 4 * AREA_BACT * (PI - 4)))
          / ( 2 * ( PI - 4 ));
}

class Bacteria
{
  public:
    double radius;      // the radius of the cell
    double max_length;  // the maximum length of the cell
    double inertia;     // the moment of inertia. for these bacteria (which 
                        // cylindrical rods in 2D), I've calculated it to be 5.
    double growth_rate; // average growth rate fo the strain
    double length;      // length at initilization of the cell
    double damping;     // damping coefficient in the des simulation.
    Bacteria(double r = 0.1, double m = 2.0, double i = 5.0, double g = 0.01,
             double l = 1.0, double d = 1.0){
             radius = r; 
             max_length = m;
             inertia = i;
             growth_rate = g;
             length = l;
             damping = d;
    }
};

class EColiGFP: public Bacteria
{
  // a very simple E. Coli with GFP, the 'null' strain of the experiments.
  // uses variable BACT1_M and BACT1_G as parameters.
  public:
    EColiGFP(double r = radiusConstantArea(BACT1_M), double m = BACT1_M,  
             double i = 5.0, double g = BACT1_G, double l = 2.0,
             double d = 1.0)
                       : Bacteria(r, m, i, g, l, d){ }
};

class EColiMCh1: public Bacteria
{
  // a very simple E. Coli with mCherry, the second strain of the experiments
  // which has same  maxlength (BACT1_M) as eGFP but different growth
  // rate (BACT2_G)
  public:
    EColiMCh1(double r = radiusConstantArea(BACT1_M), double m = BACT1_M,  
              double i = 5.0, double g = BACT2_G, double l = 2.0,
              double d = 1.0)
                       : Bacteria(r, m, i, g, l, d){ }
};

class EColiMCh2: public Bacteria
{
  // a very simple E. Coli with mCherry, the third strain of the experiments
  // which has same max length (BACT1_M) as eGFP and similar growth rate
  // (BACT1_G). essentially, when Tianyi found a good strain that was similar
  // to eGFP
  public:
    EColiMCh2(double r = radiusConstantArea(BACT1_M), double m = BACT1_M,  
              double i = 5.0, double g = BACT1_G, double l = 2.0,
              double d = 1.0)
                       : Bacteria(r, m, i, g, l, d){ }
};

class EColiA22: public Bacteria
{
  // variant of the E. Coli mCherry, that is affected by high levels of A22.
  // the strain has a different  max length (BACT2_M) and growth rate (BACT2_G)
  // as eGFP. this strain is thicker and grows more slowly
  public:
    EColiA22(double r = radiusConstantArea(BACT2_M), double m = BACT2_M,  
             double i = 5.0, double g = BACT2_G, double l = 2.0,
             double d = 2.0)
                       : Bacteria(r, m, i, g, l, d){ }
};

class BSubt: public Bacteria
{
  // a B. Subtilus variant that I never fully used in simulation, but that I
  // wanted to define here to show how you can define these bacteria.
  public:
    BSubt(double r = 0.83 / 2.0 + 0.1, double m = 7.95,
          double i = 5.0, double g = 0.0039, double l = 3.00,
          double d = 1.0)
                       : Bacteria(r, m, i, g, l, d){ }
};

#endif
