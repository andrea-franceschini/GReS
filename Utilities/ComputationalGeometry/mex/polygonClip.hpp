#ifndef POLYGON_CLIP_2D_HPP
#define POLYGON_CLIP_2D_HPP

#include <vector>
#include <array>
#include <cmath>

using Point2 = std::array<double,2>;
using Polygon = std::vector<Point2>;

// intersect two lines (x1,y1)-(x2,y2) with (x3,y3)-(x4,y4)
inline void intersectLines(
    double x1, double y1, double x2, double y2,
    double x3, double y3, double x4, double y4,
    double& xInt, double& yInt)
{
    double numX = (x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4 - y3*x4);
    double numY = (x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4 - y3*x4);
    double den  = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);

    if (std::fabs(den) < 1e-10) { xInt=x2; yInt=y2; return; }
    xInt = numX/den;
    yInt = numY/den;
}

// clip polygon against one edge
inline void clipAgainstEdge(Polygon& poly,
                            double xc1,double yc1,double xc2,double yc2)
{
    Polygon out;
    const size_t n = poly.size();
    if (n==0) return;

    for(size_t i=0;i<n;++i){
        size_t k=(i+1)%n;
        const auto &P=poly[i], &Q=poly[k];
        double i_pos=(xc2-xc1)*(P[1]-yc1)-(yc2-yc1)*(P[0]-xc1);
        double k_pos=(xc2-xc1)*(Q[1]-yc1)-(yc2-yc1)*(Q[0]-xc1);
        double xInt,yInt;
        if(i_pos>=0 && k_pos>=0){ out.push_back(Q); }
        else if(i_pos<0 && k_pos>=0){
            intersectLines(P[0],P[1],Q[0],Q[1],xc1,yc1,xc2,yc2,xInt,yInt);
            out.push_back({xInt,yInt}); out.push_back(Q);
        }
        else if(i_pos>=0 && k_pos<0){
            intersectLines(P[0],P[1],Q[0],Q[1],xc1,yc1,xc2,yc2,xInt,yInt);
            out.push_back({xInt,yInt});
        }
    }
    poly.swap(out);
}

// remove duplicates and check validity
inline bool validatePolygon(Polygon& poly){
    const double tol=1e-8;
    if(poly.size()<3) return false;
    Polygon filtered; filtered.reserve(poly.size());
    filtered.push_back(poly.front());
    for(size_t i=1;i<poly.size();++i){
        double dx=poly[i][0]-filtered.back()[0];
        double dy=poly[i][1]-filtered.back()[1];
        if(std::sqrt(dx*dx+dy*dy)>tol) filtered.push_back(poly[i]);
    }
    if(filtered.size()>1){
        double dx=filtered.front()[0]-filtered.back()[0];
        double dy=filtered.front()[1]-filtered.back()[1];
        if(std::sqrt(dx*dx+dy*dy)<tol) filtered.pop_back();
    }
    if(filtered.size()<3) return false;
    poly.swap(filtered);
    return true;
}

// clip polygon with convex clipper polygon
inline bool clipPolygon(Polygon& poly,const Polygon& clipper){
    for(size_t i=0;i<clipper.size();++i){
        size_t k=(i+1)%clipper.size();
        clipAgainstEdge(poly,clipper[i][0],clipper[i][1],clipper[k][0],clipper[k][1]);
        if(poly.empty()) return false;
    }
    return validatePolygon(poly);
}

#endif // POLYGON_CLIP_2D_HPP

