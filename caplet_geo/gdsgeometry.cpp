/*
CREATED : Jan 31, 2013
AUTHOR  : Yu-Chung Hsiao
EMAIL   : project.caplet@gmail.com

This file is part of CAPLET.

CAPLET is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CAPLET is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CAPLET.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "gdsgeometry.h"

#include "debug.h"

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;

//****
//*
//* Point
//*
//*
Point::Point()
    : x(0), y(0), z(0), dir(FLAT), len(0)
{ }

Point::Point(int xx, int yy)
    : x(xx), y(yy), z(0), dir(FLAT), len(0)
{ }

Point::Point(int xx, int yy, int zz)
    : x(xx), y(yy), z(zz), dir(FLAT), len(0)
{ }

Point Point::operator-(const Point& other)
{
    return Point( x - other.x, y - other.y, z - other.z );
}

bool Point::operator ==(const Point& other)
{
    if ( x==other.x && y==other.y && z==other.z ){
        return true;
    }else{
        return false;
    }
}

//****
//*
//* Polygon
//*
//*
Polygon::Polygon( const PointList::allocator_type &allo )
    : PointList(allo)
{ }

Polygon::Polygon(
        Polygon::size_type              n,
        const Point                     &value,
        const Polygon::allocator_type   &allo)
    : PointList(n, value, allo)
{ }

Polygon::Polygon(
        Polygon::iterator               first,
        Polygon::iterator               last,
        const Polygon::allocator_type   &allo)
    : PointList(first, last, allo)
{ }

Polygon::Polygon( const Polygon& poly)
    : PointList(poly)
{ }

//**
//* isManhattan
//* - return true if it is a Manhattan polygon
bool Polygon::isManhattan() const
{
    //* consider the case which is at least a rectangle and odd number of vertices
    if (this->size()<4){
        return false;
    }

    list<Point>::const_iterator eachVertexIt = ++this->begin();
    list<Point>::const_iterator prevVertexIt =   this->begin();
    bool dirFlag = ( (*prevVertexIt).x == this->back().x ) ? true : false;
    for ( ; eachVertexIt != this->end();
         ++eachVertexIt, ++prevVertexIt)
    {
        if (dirFlag==true && (*prevVertexIt).y==(*eachVertexIt).y ){
            dirFlag = false;
            continue;
        }else if (dirFlag==false && (*prevVertexIt).x==(*eachVertexIt).x ){
            dirFlag = true;
            continue;
        }else{
            return false;
        }
    }
    return true;
}



//****
//*
//* Rectangle
//*
//*

//**
//* Default constructor
Rectangle::Rectangle()
    : normal(FLAT), x1(0), x2(0), y1(0), y2(0), z1(0), z2(0)
{ }

Rectangle::Rectangle(Dir normal, int x1, int x2, int y1, int y2, int z1, int z2)
    : normal(normal), x1(x1), x2(x2), y1(y1), y2(y2), z1(z1), z2(z2)
{ }

//**
//* Copy constructor: default shallow copy
Rectangle::Rectangle( const Rectangle &rect )
{
    normal = rect.normal;
    x1     = rect.x1;
    x2     = rect.x2;
    y1     = rect.y1;
    y2     = rect.y2;
    z1     = rect.z1;
    z2     = rect.z2;
}

//**
//* Constructor
Rectangle::Rectangle( const Polygon &poly ) throw (ShapeTransformationError)
    : normal(ZP), z1(0), z2(0)
{
    if ( poly.size() > 5 ) {
        stringstream ss;
        ss << "Polygon contains more than five points, size = " << poly.size();
        throw ShapeTransformationError(ss.str());
    }
    list<Point>::const_iterator p1 = poly.begin();
    list<Point>::const_iterator p2 = ++(poly.begin());
    list<Point>::const_reverse_iterator p4 = poly.rbegin();

    if ( p1->x == p2->x ) {
        //* if p1 and p2 are of the same x coord
        //* then p4 will be of different x coord
        if ( p1->x < p4->x ){
            x1 = p1->x;
            x2 = p4->x;
        }else{
            x1 = p4->x;
            x2 = p1->x;
        }
        if ( p1->y < p2->y ){
            y1 = p1->y;
            y2 = p2->y;
        }else{
            y1 = p2->y;
            y2 = p1->y;
        }
    }else{
        //* or the other way around
        if ( p1->x < p2->x ){
            x1 = p1->x;
            x2 = p2->x;
        }else{
            x1 = p2->x;
            x2 = p1->x;
        }
        if ( p1->y < p4->y ){
            y1 = p1->y;
            y2 = p4->y;
        }else{
            y1 = p4->y;
            y2 = p1->y;
        }
    }
}

//**
//* area
//* - negative value means not a rectangle
double Rectangle::area() const{
    //* compute area
    //* convert to double to avoid possible overflow
    if (z1==z2){
        return (double)(x2-x1) * (double)(y2-y1);
    }
    if (x1==x2){
        return (double)(y2-y1) * (double)(z2-z1);
    }
    if (y1==y2){
        return (double)(x2-x1) * (double)(z2-z1);
    }
    return -1;
}

//**
//* geometry difference (only for those which result in rectangles)
//* - 2D in x-y plane so far (same z)
Rectangle Rectangle::operator- (const Rectangle &rect ) const{
    //* z-dir
    if ( (this->normal==ZP||this->normal==ZM)&&(rect.normal==ZP||rect.normal==ZM)
            && (this->z1 == rect.z1) && this->isOverlapping3d(rect)==true ){

        //* check: rect has no corner inside this
        if ( rect.hasCornerInside(*this)==true ){
            cerr << "WARNING: Rectangle::operator- encounters corner inside cases " << endl;
        }

        Rectangle diffRect(*this);
        if ( this->x1 < rect.x1 && rect.x1 < this->x2 ){
            diffRect.x2 = rect.x1;
        }
        else if ( this->x1 < rect.x2 && rect.x2 < this->x2 ){
            diffRect.x1 = rect.x2;
        }
        else if ( this->y1 < rect.y1 && rect.y1 < this->y2 ){
            diffRect.y2 = rect.y1;
        }
        else if ( this->y1 < rect.y2 && rect.y2 < this->y2 ){
            diffRect.y1 = rect.y2;
        }
        return diffRect;
    }
    return *this;
}


//**
//* isOverlapping
//* - consider x-y plane only
bool Rectangle::isOverlapping( const Rectangle &rect ) const{
    return x2 > rect.x1  &&  x1 < rect.x2  &&  y2 > rect.y1  &&  y1 < rect.y2;
}

//**
//* isOverlapping3d
//* - consider rect overlapping in 3d
bool Rectangle::isOverlapping3d( const Rectangle &rect ) const{
    //* z-dir
    if ( (this->z1==this->z2)&&(rect.z1==rect.z2)&&(this->z1==rect.z1) ){
        return x2 > rect.x1  &&  x1 < rect.x2  &&  y2 > rect.y1  &&  y1 < rect.y2;
    }
    //* x-dir
    if ( (this->x1==this->x2)&&(rect.x1==rect.x2)&&(this->x1==rect.x1) ){
        return y2 > rect.y1  &&  y1 < rect.y2  &&  z2 > rect.z1  &&  z1 < rect.z2;
    }
    //* y-dir
    if ( (this->y1==this->y2)&&(rect.y1==rect.y2)&&(this->y1==rect.y1) ){
        return x2 > rect.x1  &&  x1 < rect.x2  &&  z2 > rect.z1  &&  z1 < rect.z2;
    }
    return false;
}


//**
//* hasCornerInside
//* - return true if this has a corner in rect
bool Rectangle::hasCornerInside ( const Rectangle &rect ) const{

    //* Check if any corner of *thisRectIt is inside viaRect
    if ( (rect.x1 < this->x1 && this->x1 < rect.x2)||
         (rect.x1 < this->x2 && this->x2 < rect.x2)){
        if ( (rect.y1 < this->y1 && this->y1 < rect.y2)||
             (rect.y1 < this->y2 && this->y2 < rect.y2)){
            return true;
        }
    }
    return false;
}
//**
//* print
void Rectangle::print() const{
    cout << "Rectangle: " << endl;
    cout << "  normal = " << normal << endl;
    cout << "  area   = " << area() << endl;
    cout << "  ( " << x1 << ", " << x2 << ", " << y1 << ", " << y2 << ", " << z1 << ", " << z2 << " )" << endl;
}


//****
//*
//* RectangleList
//*
//*

RectangleList::RectangleList(const allocator_type &allo)
    :list<Rectangle>(allo)
{ }

RectangleList::RectangleList(size_type n, const Rectangle &value, const allocator_type &allo)
    :list<Rectangle>(n, value, allo)
{ }

RectangleList::RectangleList( iterator first, iterator last, const allocator_type &allo)
    :list<Rectangle>(first, last, allo)
{ }

RectangleList::RectangleList(const RectangleList &rectList)
    :list<Rectangle>()
{
    this->insert(this->begin(), rectList.begin(), rectList.end());
}

//**
//* RectangleList::merge
//* - CURRENTLY DOES NOT SUPPORT SUBLAYERS
//* - looks fine but seems not optimal (some may be combined but do not).
//* - works for all Manhattan directions
void RectangleList::merge()
{
    for ( RectangleList::iterator rectIit = this->begin();
          rectIit != --this->end(); ++rectIit )
    {
        RectangleList::iterator rectJit = rectIit;
        for ( ++rectJit;
              rectJit != this->end(); )
        {
            //* only combine rects with the same normal
            if ( rectIit->normal != rectJit->normal || rectIit==rectJit){
                ++rectJit;
                continue;
            }

            bool flagDeleteRectJ = false;

            //* same x range
            if ( rectIit->x1 == rectJit->x1 && rectIit->x2 == rectJit->x2 ){
                //* adjacent from topview above
                if ( rectIit->y2 == rectJit->y1 ){
                    rectIit->y2 = rectJit->y2;
                    flagDeleteRectJ = true;
                }
                //* adjacent from topview below
                else if ( rectIit->y1 == rectJit->y2 ){
                    rectIit->y1 = rectJit->y1;
                    flagDeleteRectJ = true;
                }
            }
            //* same y range
            else if ( rectIit->y1 == rectJit->y1 && rectIit->y2 == rectJit->y2 ){
                //* adjacent from topview right
                if ( rectIit->x2 == rectJit->x1 ){
                    rectIit->x2 = rectJit->x2;
                    flagDeleteRectJ = true;
                }
                //* adjacent from topview left
                else if( rectIit->x1 == rectJit->x2 ){
                    rectIit->x1 = rectJit->x1;
                    flagDeleteRectJ = true;
                }
            }

            if ( flagDeleteRectJ==true ){
                this->erase(rectJit);
                rectJit = this->begin();
            }else{
                ++rectJit;
            }
        }
    }
}

//**
//* RectangleList::decompose
//* - modify rectList of x-y plane rectangles
//*   such that the resulting rectList consists of disjoint rectangles
void decomposeXdir(Rectangle& rectI, Rectangle &rectJ, RectangleList &decomposedRectJ);
void decomposeYdir(Rectangle& rectI, Rectangle &rectJ, RectangleList &decomposedRectJ);
void RectangleList::decompose(){
    //* if fewer than two rects, no need to decompose
    if (this->size()<2){
        return;
    }

    RectangleMap rectMap;

    for ( RectangleList::iterator eachRectIt = this->begin();
          eachRectIt != this->end(); ++eachRectIt )
    {
        rectMap.insert(pair<double, Rectangle>( (*eachRectIt).area(), *eachRectIt ));
    }
//    printRectMap(rectMap);


    //* use rectI to decompose rectJ
    for ( RectangleMap::iterator itemIit = rectMap.begin();
          itemIit!=--rectMap.end(); ++itemIit )
    {
        Rectangle& rectI = (*itemIit).second;

        RectangleList decomposedRectJ;
        RectangleMap::iterator itemJit = itemIit;
        for ( ++itemJit; itemJit!=rectMap.end(); rectMap.erase(itemJit++) )
        {
            bool flagDecomposeXFirst = false;
            Rectangle& rectJ = (*itemJit).second;

            //* if rectI and rectJ overlap
            if ( rectI.isOverlapping(rectJ) )
            {
                int overlapLenX = min(rectI.x2, rectJ.x2) - max(rectI.x1, rectJ.x1);
                int overlapLenY = min(rectI.y2, rectJ.y2) - max(rectI.y1, rectJ.y1);

                int remainingLenX = rectJ.x2 - rectJ.x1 - overlapLenX;
                int remainingLenY = rectJ.y2 - rectJ.y1 - overlapLenY;

                if ( overlapLenX == rectI.x2 - rectI.x1 ){
                    flagDecomposeXFirst = true;
                }else if( overlapLenY == rectI.y2 - rectI.y1 ){
                    flagDecomposeXFirst = false;
                }else if( remainingLenY < remainingLenX ){
                    flagDecomposeXFirst = true;
                }

                //* perform decomposition
                if ( flagDecomposeXFirst == true ){
                    //* decompose x-dir first
                    decomposeXdir(rectI, rectJ, decomposedRectJ);
                    //* decompose y-dir later
                    decomposeYdir(rectI, rectJ, decomposedRectJ);
                }else{
                    //* decompose y-dir first
                    decomposeYdir(rectI, rectJ, decomposedRectJ);
                    //* decompose x-dir later
                    decomposeXdir(rectI, rectJ, decomposedRectJ);
                }

            }else{ //* no overlap
                decomposedRectJ.push_back(rectJ);
            }
        }

        //* insert decomposed results
        for ( RectangleList::iterator eachRectIt = decomposedRectJ.begin();
              eachRectIt != decomposedRectJ.end(); ++eachRectIt )
        {
            rectMap.insert( pair<double, Rectangle>((*eachRectIt).area(), *eachRectIt) );
        }
    }

    this->clear();
    for ( RectangleMap::iterator itemIt = rectMap.begin();
          itemIt!=rectMap.end(); ++itemIt )
    {
        this->push_back( itemIt->second );
    }
}
void decomposeXdir(Rectangle& rectI, Rectangle &rectJ, RectangleList &decomposedRectJ){
    if ( rectJ.x1 < rectI.x1 && rectI.x1 < rectJ.x2 ){
        //* if rectI.xmin is between rectJ.xmin and rectJ.xmax
        decomposedRectJ.push_back(rectJ);
        decomposedRectJ.back().x2 = rectI.x1;
        rectJ.x1 = rectI.x1;
    }
    if ( rectJ.x1 < rectI.x2 && rectI.x2 < rectJ.x2 ){
        //* if rectI.xmax is between rectJ.xmin and rectJ.xmax
        decomposedRectJ.push_back(rectJ);
        decomposedRectJ.back().x1 = rectI.x2;
        rectJ.x2 = rectI.x2;
    }
}
void decomposeYdir(Rectangle& rectI, Rectangle &rectJ, RectangleList &decomposedRectJ){
    if ( rectJ.y1 < rectI.y1 && rectI.y1 < rectJ.y2 ){
        //* if rectI.ymin is between rectJ.ymin and rectJ.ymax
        decomposedRectJ.push_back(rectJ);
        decomposedRectJ.back().y2 = rectI.y1;
        rectJ.y1 = rectI.y1;
    }
    if ( rectJ.y1 < rectI.y2 && rectI.y2 < rectJ.y2 ){
        //* if rectI.ymax is between rectJ.ymin and rectJ.ymax
        decomposedRectJ.push_back(rectJ);
        decomposedRectJ.back().y1 = rectI.y2;
        rectJ.y2 = rectI.y2;
    }
}


//****
//*
//* Conductor
//*
//*


Conductor::Conductor()
    : nMetal(0), nVia(0), nLayer(0), layer(0)
{ }

Conductor::Conductor(int nMetal1, int nVia1)
    : nMetal(nMetal1), nVia(nVia1), nLayer(nMetal1+nVia1), layer(nLayer, vector< RectangleList >(6))
{ }

//**
//* operator +=
Conductor& Conductor::operator+= ( const Conductor& rhs ) throw ( ConductorLayerNotCompatibleError )
{///
    if ( nMetal != rhs.nMetal || nVia != rhs.nVia ){
        throw ConductorLayerNotCompatibleError(nMetal, nVia, rhs.nMetal, rhs.nVia);
    }
    for ( int i=0; i<nLayer; ++i ){
        for ( int j=0; j<nDir; ++j ){
            layer[i][j].insert(layer[i][j].end(), rhs.layer[i][j].begin(), rhs.layer[i][j].end());
        }
    }
    return *this;
}

//**
//* isContaining
bool Conductor::isContaining(const Rectangle rect, const int layerIndex )
{
    //* check if rect is facing the +z dir
    if ( rect.normal != Z ){
        return false;
    }
    for ( RectangleList::iterator eachRectIt = layer[layerIndex][TOP].begin();
          eachRectIt != layer[layerIndex][TOP].end(); ++eachRectIt ){
        if ( eachRectIt->isOverlapping(rect) == true ){
            return true;
        }
    }
    return false;
}



//**
//* generateVia
//* - ASSUMPTION: no metal rect has a corner inside the via
//* - For side walls, simply generate wall rectangles
//* - For top or bottom, if isDecomposed is false, then generate rectangles
//* -                    if isDecomposed is true,
//* - For the top and the bottom metal layers,
//*   1. Search for the overlapping rectangle
//*   2. Generate four rectangles corresponding to the decomposed rectangle
//*   3. Erase the original rectangle (iterator advanced by 1)
//*      and insert the decomposed rectangles in front of the original rectangle
//*   4. Erase the decomposed rectangles with 0 area
//*   o. Modify via size if that exceeds the connected metal boundary (SHRINK_VIA)
//*      This is the only place that modifies rect.
Rectangle innerDecompose(RectangleList &rectList, RectangleList::iterator &thisRectIt, Rectangle viaRect);

void Conductor::generateVia(
        Rectangle       rect,
        const int       viaIndex,
        const int       *const *viaDef,
        const int       *const *viaConnect,
        const bool      isDecomposed )
{
    if ( viaIndex >= nVia ){
        //* viaIndex out of range
        return;
    }
    int viaTotalIndex       = viaIndex + nMetal;
    int bottomHeight        = viaDef[viaIndex][0];
    int topHeight           = viaDef[viaIndex][1];
    int bottomMetalIndex    = viaConnect[viaIndex][0];
    int topMetalIndex       = viaConnect[viaIndex][1];

    //* generate top and bottom rectangle without decomposition
    if (isDecomposed==false) {
        //* generate top rectangle
        layer[viaTotalIndex][TOP].push_back(rect);
        Rectangle& topRect    = layer[viaTotalIndex][TOP].back();
        topRect.z1            = topHeight;
        topRect.z2            = topHeight;
        topRect.normal        = ZP;

        //* generate bottom rectangle
        layer[viaTotalIndex][BOTTOM].push_back(rect);
        Rectangle& bottomRect = layer[viaTotalIndex][BOTTOM].back();
        bottomRect.z1         = bottomHeight;
        bottomRect.z2         = bottomHeight;
        bottomRect.normal     = ZM;
    }
    //* generate top and bottom rectangle with decomposition
    else{
        int heightArray[] = {topHeight, bottomHeight};
        int metalIndexArray[] = {topMetalIndex, bottomMetalIndex};
        Dir metalRectDirArray[] = {BOTTOM, TOP};
        Dir viaRectDirArray[] = {TOP, BOTTOM};
        ::Dir viaRectNormalArray[] = {ZP, ZM};
        int numberOfCases = 2;

        for (int i=0; i<numberOfCases; ++i){
            int         height = heightArray[i];
            int         metalIndex = metalIndexArray[i];
            Dir         metalRectDir = metalRectDirArray[i];
            Dir         viaRectDir = viaRectDirArray[i];
            ::Dir       viaRectNormal = viaRectNormalArray[i];

            Rectangle   tempRect(rect);

            RectangleList::iterator end = layer[metalIndex][metalRectDir].end();
            for ( RectangleList::iterator eachRectIt = layer[metalIndex][metalRectDir].begin();
                  eachRectIt != end; ){
                if ( eachRectIt->isOverlapping(tempRect) == false ){
                    ++eachRectIt;
                    continue;
                }
                //* rect overlaps with *eachRectIt
                //* Decompose *eachRectIt
                //* Note that *eachRectIt is erased which implies eachRectIt is advanced by 1
                tempRect = innerDecompose(layer[metalIndex][metalRectDir], eachRectIt, tempRect);
                if ( tempRect.area()==0 ){
                    break;
                }
            }
            //* Cover the via if it does not connect to a metal layer
            if ( tempRect.area()!=0 ){
                layer[viaTotalIndex][viaRectDir].push_back(tempRect);
                Rectangle& rectRef    = layer[viaTotalIndex][viaRectDir].back();
                rectRef.z1            = height;
                rectRef.z2            = height;
                rectRef.normal        = viaRectNormal;
            }
        }//END OF FOR
    }//END OF ELSE

    //* generate side walls
    //* generate left rectangle
    layer[viaTotalIndex][LEFT].push_back(Rectangle(rect));
    Rectangle& leftRect = layer[viaTotalIndex][LEFT].back();
    leftRect.x2 = rect.x1;
    leftRect.z1 = bottomHeight;
    leftRect.z2 = topHeight;
    leftRect.normal = XM;

    //* generate right rectangle
    layer[viaTotalIndex][RIGHT].push_back(Rectangle(rect));
    Rectangle& rightRect = layer[viaTotalIndex][RIGHT].back();
    rightRect.x1 = rect.x2;
    rightRect.z1 = bottomHeight;
    rightRect.z2 = topHeight;
    rightRect.normal = XP;

    //* generate back rectangle
    layer[viaTotalIndex][BACK].push_back(Rectangle(rect));
    Rectangle& backRect = layer[viaTotalIndex][BACK].back();
    backRect.y2 = rect.y1;
    backRect.z1 = bottomHeight;
    backRect.z2 = topHeight;
    backRect.normal = YM;

    //* generate front rectangle
    layer[viaTotalIndex][FRONT].push_back(Rectangle(rect));
    Rectangle& frontRect = layer[viaTotalIndex][FRONT].back();
    frontRect.y1 = rect.y2;
    frontRect.z1 = bottomHeight;
    frontRect.z2 = topHeight;
    frontRect.normal = YP;
}

//**
//* innerDecompose
//* - Aux function for Conductor::generateVia
//* - Decompose a rectangle by carving out via
//* - additional rectangles are inserted before thisRectIt
//* - thisRectIt iterator will be advanced by one
Rectangle innerDecompose(RectangleList &rectList, RectangleList::iterator &thisRectIt, Rectangle viaRect){
    int xlen = thisRectIt->x2 - thisRectIt->x1;
    int ylen = thisRectIt->y2 - thisRectIt->y1;

    //* Check if any corner of *thisRectIt is inside viaRect
    if ( thisRectIt->hasCornerInside(viaRect)==true ){
        cerr << "WARNING: generateVia needs to handle the corners inside via" << endl;
    }

    //* Check if any edge is inside viaRect
    Rectangle rectOutside(viaRect);

    if ( viaRect.x1 < thisRectIt->x1 && thisRectIt->x1 < viaRect.x2 ){
        rectOutside.x2 = thisRectIt->x1;
        viaRect.x1     = thisRectIt->x1;
    }
    else if ( viaRect.x1 < thisRectIt->x2 && thisRectIt->x2 < viaRect.x2 ){
        rectOutside.x1 = thisRectIt->x2;
        viaRect.x2     = thisRectIt->x2;
    }
    else if ( viaRect.y1 < thisRectIt->y1 && thisRectIt->y1 < viaRect.y2 ){
        rectOutside.y2 = thisRectIt->y1;
        viaRect.y1     = thisRectIt->y1;
    }
    else if ( viaRect.y1 < thisRectIt->y2 && thisRectIt->y2 < viaRect.y2 ){
        rectOutside.y1 = thisRectIt->y2;
        viaRect.y2     = thisRectIt->y2;
    }
    else{
        //* viaRect is totally contained in *thisRectIt
        rectOutside = Rectangle();
    }

    //* For viaRect which is completely contained in *thisRectIt
    //* left rect
    RectangleList::iterator leftRectIt = rectList.insert(thisRectIt, *thisRectIt);
    leftRectIt  ->x2 = viaRect.x1;
    //* right rect
    RectangleList::iterator rightRectIt = rectList.insert(thisRectIt, *thisRectIt);
    rightRectIt ->x1 = viaRect.x2;
    //* bottom rect
    RectangleList::iterator bottomRectIt = rectList.insert(thisRectIt, *thisRectIt);
    bottomRectIt->y2 = viaRect.y1;
    //* top rect
    RectangleList::iterator topRectIt = rectList.insert(thisRectIt, *thisRectIt);
    topRectIt   ->y1 = viaRect.y2;

    if ( xlen > ylen ){
        //* cut through y-dir
        //* bottom rect
        bottomRectIt->x1 = viaRect.x1;
        bottomRectIt->x2 = viaRect.x2;
        //* top rect
        topRectIt   ->x1 = viaRect.x1;
        topRectIt   ->x2 = viaRect.x2;
    }
    else{// xlen <= ylen
        //* cut through x-dir
        //* left rect
        leftRectIt  ->y1 = viaRect.y1;
        leftRectIt  ->y2 = viaRect.y2;
        //* right rect
        rightRectIt ->y1 = viaRect.y1;
        rightRectIt ->y2 = viaRect.y2;
    }
    //* remove *thisRectIt
    thisRectIt = rectList.erase(thisRectIt);

    //* remove rects with 0 area
    if ( leftRectIt->area()==0 ){
        rectList.erase(leftRectIt);
    }
    if ( rightRectIt->area()==0 ){
        rectList.erase(rightRectIt);
    }
    if ( bottomRectIt->area()==0 ){
        rectList.erase(bottomRectIt);
    }
    if ( topRectIt->area()==0 ){
        rectList.erase(topRectIt);
    }
    return rectOutside;
}



void Conductor::print() const {
    cout << "Conductor: " << endl;
    cout << "  nMetal = " << nMetal << endl;
    cout << "  nVia   = " << nVia   << endl;
    cout << "  nLayer = " << nLayer << endl;

    for ( int i=0; i<nLayer; ++i ){
        cout << "    Layer: " << i << endl;
        for ( int j=0; j<nDir; ++j ){
            cout << "      Dir : " << j << "  size: " << layer[i][j].size() << endl;
#ifdef CONDUCTOR_PRINT_ALL_RECTS
            for ( RectangleList::const_iterator eachRectIt = layer[i][j].begin();
                  eachRectIt != layer[i][j].end(); ++eachRectIt )
            {
                eachRectIt->print();
            }
#endif
        }
    }
}

//**
//* checkSelfOverlapping
bool Conductor::checkSelfOverlapping(const int *const *viaConnect)
{
    //* check metal/via overlap within the same layer
    const int dirArrayLength = 6;
    for ( int i=0; i<nLayer; ++i ){
        for ( int dirIndex=0; dirIndex<dirArrayLength; ++dirIndex ){
            RectangleList::const_iterator end1 = --layer[i][dirIndex].end();
            RectangleList::const_iterator end2 =   layer[i][dirIndex].end();
            for ( RectangleList::const_iterator each1 = layer[i][dirIndex].begin();
                  each1 != end1; ++each1 ){
                RectangleList::const_iterator each2 = each1;
                for (  ++each2; each2 != end2; ++each2){
                    if (each1->isOverlapping3d(*each2)==true) {
                        each1->print();
                        each2->print();
                        cout << "Overlap:" << endl;
                        cout << endl;
                        return true;
                    }
                }
            }
        }
    }
    //* check metal/via interfaces
    const Dir    dirViaArray  [] = {TOP, BOTTOM};
    const Dir    dirMetalArray[] = {BOTTOM, TOP};
    const int    dirInterfaceArrayLength = 2;
    const int    viaConnectIndexArray[] = {1, 0};
    for ( int totalIndex=nMetal; totalIndex<nLayer; ++totalIndex ){
        const int viaIndex = totalIndex-nMetal;
        for ( int dirIndex=0; dirIndex<dirInterfaceArrayLength; ++dirIndex ){
            Dir dirVia   = dirViaArray[dirIndex];
            Dir dirMetal = dirMetalArray[dirIndex];
            const int metalIndex = viaConnect[viaIndex][viaConnectIndexArray[dirIndex]];
            for ( RectangleList::const_iterator eachViaRectIt = layer[totalIndex][dirVia].begin();
                  eachViaRectIt != layer[totalIndex][dirVia].end(); ++eachViaRectIt ){
                for (   RectangleList::const_iterator eachMetalRectIt
                        = layer[metalIndex][dirMetal].begin();
                        eachMetalRectIt != layer[metalIndex][dirMetal].end();
                        ++eachMetalRectIt){
                    if (eachViaRectIt->isOverlapping3d(*eachMetalRectIt)==true) {
                        cout << "Overlap:" << endl;
                        eachViaRectIt->print();
                        eachMetalRectIt->print();
                        cout << endl;

                        layer[nLayer-1][TOP].push_back(*eachViaRectIt);
                        layer[nLayer-1][TOP].back().z1 = 1000;
                        layer[nLayer-1][TOP].back().z2 = 1000;
                        layer[nLayer-1][TOP].push_back(*eachMetalRectIt);
                        layer[nLayer-1][TOP].back().z1 = 1000;
                        layer[nLayer-1][TOP].back().z2 = 1000;
                        return true;
                    }
                }
            }
        }

    }
    return false;
}


//**
//* checkZeroAreaRectangle
//* - return true if this has any rectangle with zero area
bool Conductor::checkZeroAreaRectangle()
{
    for (int i=0; i<nLayer; ++i) {
        for (int j=0; j<nDir; ++j) {
            RectangleList::const_iterator end = layer[i][j].end();
            for (RectangleList::const_iterator each = layer[i][j].begin();
                    each != end; ++each){
                if (each->area()==0) {
                    each->print();
                    return true;
                }
            }
        }
    }
    return false;
}


//****
//*
//* RectangleGL
//*
//*

//* RectangleGL
RectangleGL::RectangleGL()
    : xn(0.0), yn(0.0), zn(0.0), x1(0.0), x2(0.0), y1(0.0), y2(0.0), z1(0.0), z2(0.0),
      shapeType(FLAT_TYPE), shapeDir(FLAT_SHAPE), shapeNormalDistance(0.0), shapeShift(0.0)
{ }

RectangleGL::RectangleGL(const Rectangle & rect, float unit)
    : xn(0.0), yn(0.0), zn(0.0),
      x1(rect.x1*unit), x2(rect.x2*unit),
      y1(rect.y1*unit), y2(rect.y2*unit),
      z1(rect.z1*unit), z2(rect.z2*unit),
      shapeType(FLAT_TYPE), shapeDir(FLAT_SHAPE), shapeNormalDistance(0.0), shapeShift(0.0)
{
    switch(rect.normal){
    case XP:
        xn= 1.0;
        break;
    case XM:
        xn=-1.0;
        break;
    case YP:
        yn= 1.0;
        break;
    case YM:
        yn=-1.0;
        break;
    case ZP:
        zn= 1.0;
        break;
    case ZM:
        zn=-1.0;
        break;
    case FLAT:
        zn= 1.0;
    default:
        cerr << "ERROR: impossible direction" << endl;
    }
}

bool RectangleGL::isOverlappingProjection(const RectangleGL &rect) const
{
    if (this->zn!=0 && rect.zn!=0){
        //* z-dir
        return x2 > rect.x1  &&  x1 < rect.x2  &&  y2 > rect.y1  &&  y1 < rect.y2;
    }
    if (this->xn!=0 && rect.xn!=0){
        //* x-dir
        return y2 > rect.y1  &&  y1 < rect.y2  &&  z2 > rect.z1  &&  z1 < rect.z2;
    }
    if (this->yn!=0 && rect.yn!=0){
        //* y-dir
        return z2 > rect.z1  &&  z1 < rect.z2  &&  x2 > rect.x1  &&  x1 < rect.x2;
    }

    //* debug info
    cerr << "ERROR: impossible direction in RectangleGL::isOverlappingProjection()" << endl;
    cerr << "(" << xn << "," << yn << "," << zn << ")" << endl;
    return false;
}

bool RectangleGL::isOverlapping(const RectangleGL &rect) const
{
    if (this->zn!=0 && rect.zn!=0 && z1==rect.z1){
        //* z-dir
        return x2 > rect.x1  &&  x1 < rect.x2  &&  y2 > rect.y1  &&  y1 < rect.y2;
    }
    if (this->xn!=0 && rect.xn!=0 && x1==rect.x1){
        //* x-dir
        return y2 > rect.y1  &&  y1 < rect.y2  &&  z2 > rect.z1  &&  z1 < rect.z2;
    }
    if (this->yn!=0 && rect.yn!=0 && y1==rect.y1){
        //* y-dir
        return z2 > rect.z1  &&  z1 < rect.z2  &&  x2 > rect.x1  &&  x1 < rect.x2;
    }

    return false;
}
bool RectangleGL::isOverlappingOrEdgeNeighboring(const RectangleGL &rect) const
{
    if (this->zn!=0 && z1==rect.z1){
        //* z-dir
        return (x2 >= rect.x1  &&  x1 <= rect.x2  &&  y2 > rect.y1  &&  y1 < rect.y2)||
               (x2 > rect.x1  &&  x1 < rect.x2  &&  y2 >= rect.y1  &&  y1 <= rect.y2)  ;
    }
    if (this->xn!=0 && x1==rect.x1){
        //* x-dir
        return (y2 >= rect.y1  &&  y1 <= rect.y2  &&  z2 > rect.z1  &&  z1 < rect.z2)||
               (y2 > rect.y1  &&  y1 < rect.y2  &&  z2 >= rect.z1  &&  z1 <= rect.z2)  ;
    }
    if (this->yn!=0 && y1==rect.y1){
        //* y-dir
        return (z2 >= rect.z1  &&  z1 <= rect.z2  &&  x2 > rect.x1  &&  x1 < rect.x2)||
               (z2 > rect.z1  &&  z1 < rect.z2  &&  x2 >= rect.x1  &&  x1 <= rect.x2)  ;
    }

    return false;
}

bool RectangleGL::isEmpty() const
{
    return ( xn==0 && yn==0 && zn==0 );
}

bool RectangleGL::operator ==(const RectangleGL &rect) const
{
    return (xn==rect.xn && yn==rect.yn && zn==rect.zn &&
            x1==rect.x1 && x2==rect.x2 && y1==rect.y1 && y2==rect.y2 && z1==rect.z1 && z2==rect.z2);
}

bool RectangleGL::isCoincidental(const RectangleGL &support, const float margin) const
{
    float xmargin = margin*(support.x2-support.x1);
    float ymargin = margin*(support.y2-support.y1);

    return (xn==support.xn && yn==support.yn && zn==support.zn &&

            support.x1 <= x1 && x1 <= support.x1+xmargin &&
            support.x2-xmargin <= x2 && x2 <= support.x2 &&

            support.y1 <= y1 && y1 <= support.y1+ymargin &&
            support.y2-ymargin <= y2 && y2 <= support.y2 &&

            z1==support.z1 && z2==support.z2 );
}

//**
//* isContaining
//* - 3D
//* - same direction and elevation
bool RectangleGL::isContaining(const RectangleGL &rect) const
{
    if ( (zn!=0 && zn==rect.zn && z1==rect.z1 && x1<=rect.x1 && rect.x2<=x2 && y1<=rect.y1 && rect.y2<=y2) ||
         (yn!=0 && yn==rect.yn && y1==rect.y1 && z1<=rect.z1 && rect.z2<=z2 && x1<=rect.x1 && rect.x2<=x2) ||
         (xn!=0 && xn==rect.xn && x1==rect.x1 && y1<=rect.y1 && rect.y2<=y2 && z1<=rect.z1 && rect.z2<=z2) ){
        return true;
    }
    return false;
}

RectangleGL RectangleGL::intersectProjection(const RectangleGL &rect) const
{
    RectangleGL result;
    result.xn = this->xn;
    result.yn = this->yn;
    result.zn = this->zn;

    if (result.xn!=0){
        //* normal x-dir
        result.x1 = this->x1;
        result.x2 = this->x2;
        result.y1 = max(this->y1, rect.y1);
        result.y2 = min(this->y2, rect.y2);
        result.z1 = max(this->z1, rect.z1);
        result.z2 = min(this->z2, rect.z2);
    }
    else if (result.yn!=0){
        //* normal y-dir
        result.x1 = max(this->x1, rect.x1);
        result.x2 = min(this->x2, rect.x2);
        result.y1 = this->y1;
        result.y2 = this->y2;
        result.z1 = max(this->z1, rect.z1);
        result.z2 = min(this->z2, rect.z2);
    }
    else if (result.zn!=0){
        //* normal z-dir
        result.x1 = max(this->x1, rect.x1);
        result.x2 = min(this->x2, rect.x2);
        result.y1 = max(this->y1, rect.y1);
        result.y2 = min(this->y2, rect.y2);
        result.z1 = this->z1;
        result.z2 = this->z2;
    }
//    if (*this == result){ // marked as empty for RectangleGL.isEmpty()
//        result.xn=0;
//        result.yn=0;
//        result.zn=0;
//    }
    result.shapeShift = 1; // marked as inserted rect
    return result;
}


RectangleGL RectangleGL::intersectArchOnFlat(const RectangleGL &flat) const
{
    RectangleGL result(*this);

    //* only consider the intersection with flat
    if (flat.shapeType != FLAT_TYPE){
        result.xn = 0;
        result.yn = 0;
        result.zn = 0;
        return result;
    }

    //* only consider overlapping cases
    if (this->isOverlapping(flat)==false){
        result.xn = 0;
        result.yn = 0;
        result.zn = 0;
        return result;
    }

    //* only consider the edge of decaying head is contained in flat
    if ( (shapeDir == X_DECAY && shapeNormalDistance > 0 && !( flat.x1 <= x1 && x1 < flat.x2 )) || //* decay in +x
         (shapeDir == X_DECAY && shapeNormalDistance < 0 && !( flat.x1 < x2 && x2 <= flat.x2 )) || //* decay in -x
         (shapeDir == Y_DECAY && shapeNormalDistance > 0 && !( flat.y1 <= y1 && y1 < flat.y2 )) || //* decay in +y
         (shapeDir == Y_DECAY && shapeNormalDistance < 0 && !( flat.y1 < y2 && y2 <= flat.y2 )) || //* decay in -y
         (shapeDir == Z_DECAY && shapeNormalDistance > 0 && !( flat.z1 <= z1 && z1 < flat.z2 )) || //* decay in +z
         (shapeDir == Z_DECAY && shapeNormalDistance < 0 && !( flat.z1 < z2 && z2 <= flat.z2 )) ){ //* decay in -z
        result.xn = 0;
        result.yn = 0;
        result.zn = 0;
        return result;
    }

    if (result.xn!=0 ){
        //* normal x-dir
        result.y1 = max(this->y1, flat.y1);
        result.y2 = min(this->y2, flat.y2);
        result.z1 = max(this->z1, flat.z1);
        result.z2 = min(this->z2, flat.z2);
    }
    else if (result.yn!=0){
        //* normal y-dir
        result.x1 = max(this->x1, flat.x1);
        result.x2 = min(this->x2, flat.x2);
        result.z1 = max(this->z1, flat.z1);
        result.z2 = min(this->z2, flat.z2);
    }
    else if (result.zn!=0){
        //* normal z-dir
        result.x1 = max(this->x1, flat.x1);
        result.x2 = min(this->x2, flat.x2);
        result.y1 = max(this->y1, flat.y1);
        result.y2 = min(this->y2, flat.y2);
    }
    return result;
}

void RectangleGL::print() const
{
    cout << "normal = (" << xn << "," << yn << "," << zn << ")" << endl;
    cout << "coords = (" << x1 << "," << x2 << ", "
                         << y1 << "," << y2 << ", "
                         << z1 << "," << z2 << ")" << endl;
    cout << "shape  = (" << shapeType << ", " << shapeDir << ", "
                         << shapeNormalDistance << ", " << shapeShift << endl;
}

void RectangleGL::printCapletLine(ostream &out) const
{
    #ifdef COMBINED_SHAPE
    switch(this->shapeType){
    case FLAT_TYPE:
        out << "F 1 ";
        break;
    case ARCH_TYPE:
        out << "A 0 ";
        break;
    case SIDE_TYPE:
        out << "S 0 ";
        break;
    }
    #endif
    #ifndef COMBINED_SHAPE
    switch(this->shapeType){
    case FLAT_TYPE:
        out << "F 1 ";
        break;
    case ARCH_TYPE:
        out << "A 1 ";
        break;
    case SIDE_TYPE:
        out << "S 1 ";
        break;
    }
    #endif
    out << setw(14) << setprecision(6) << this->x1
        << setw(14) << setprecision(6) << this->x2
        << setw(14) << setprecision(6) << this->y1
        << setw(14) << setprecision(6) << this->y2
        << setw(14) << setprecision(6) << this->z1
        << setw(14) << setprecision(6) << this->z2;

    if (xn!=0){
        out << " 0 ";
    }
    else if (yn!=0){
        out << " 1 ";
    }
    else if (zn!=0){
        out << " 2 ";
    }

    out << this->shapeDir;
    out << setw(14) << setprecision(6) << this->shapeNormalDistance
        << setw(14) << setprecision(6) << this->shapeShift;
    out << endl;
}

void RectangleGL::printCapletLineFlat(ostream &out) const
{
    out << "F 1 ";
    out << setw(14) << setprecision(6) << this->x1
        << setw(14) << setprecision(6) << this->x2
        << setw(14) << setprecision(6) << this->y1
        << setw(14) << setprecision(6) << this->y2
        << setw(14) << setprecision(6) << this->z1
        << setw(14) << setprecision(6) << this->z2;

    if (xn!=0){
        out << " 0 ";
    }
    else if (yn!=0){
        out << " 1 ";
    }
    else if (zn!=0){
        out << " 2 ";
    }

    out << FLAT_SHAPE;
    out << setw(14) << setprecision(6) << this->shapeNormalDistance
        << setw(14) << setprecision(6) << this->shapeShift;
    out << endl;
}

//****
//*
//* ConductorFloat
//*
//*
ConductorFP::ConductorFP()
    : nMetal(0), nVia(0), nLayer(nMetal+nVia), layer(0, DirRectangleGLList(nDir, RectangleGLList()))
{ }

ConductorFP::ConductorFP(const int nMetal, const int nVia)
    : nMetal(nMetal), nVia(nVia), nLayer(nMetal+nVia), layer(nLayer, DirRectangleGLList(nDir, RectangleGLList()))
{
    cout << "ConductorFP::ConductorFP(const int nLayer)" << endl;
}

ConductorFP::ConductorFP(const Conductor &cond, const float unit)
    : nMetal(cond.nMetal), nVia(cond.nVia), nLayer(nMetal+nVia), layer(nLayer, DirRectangleGLList(nDir, RectangleGLList()))
{
    for ( int layerIndex = 0; layerIndex < nLayer; ++layerIndex ){
        for ( unsigned dirIndex = 0; dirIndex < nDir; ++dirIndex){
            for ( RectangleList::const_iterator eachRect = cond.layer[layerIndex][dirIndex].begin();
                    eachRect != cond.layer[layerIndex][dirIndex].end(); ++eachRect ){
                layer[layerIndex][dirIndex].push_back(RectangleGL(*eachRect, unit));
            }
        }
    }
}

//**
//* ConductorFP::size()
//* - return total number of RectangleGL
unsigned ConductorFP::size() const
{
    unsigned count = 0;
    for (unsigned layerIndex=0; layerIndex<layer.size(); ++layerIndex){
        for (unsigned dirIndex=0; dirIndex<ConductorFP::nDir; ++dirIndex){
            count += layer[layerIndex][dirIndex].size();
        }
    }
    return count;
}

ConductorFPList::ConductorFPList(const allocator_type &allo)
    : list<ConductorFP>(allo)
{ }

ConductorFPList::ConductorFPList(
        size_type n,
        const ConductorFP &value,
        const allocator_type &allo)
    : list<ConductorFP>(n, value, allo)
{ }

ConductorFPList::ConductorFPList(
        iterator first,
        iterator last,
        const allocator_type &allo)
    : list<ConductorFP>(first, last, allo)
{ }

ConductorFPList::ConductorFPList(const ConductorFPList &condFGList)
    : list<ConductorFP>()
{
    this->insert(this->begin(), condFGList.begin(), condFGList.end());
}

void ConductorFPList::constructFrom(const ConductorList &condList, const float unit)
{
    this->clear();
    for ( ConductorList::const_iterator eachCond = condList.begin();
            eachCond != condList.end(); ++eachCond){
        this->push_back(ConductorFP(*eachCond, unit));
    }
}

ConductorFPList::ConductorFPList(const ConductorList &condList, const float unit)
    : list<ConductorFP>()
{
    this->constructFrom(condList, unit);
}


//****
//*
//* RectangleGLList
//*
//*

RectangleGLList::RectangleGLList(const allocator_type &allo)
    : list<RectangleGL> (allo)
{ }

RectangleGLList::RectangleGLList(size_type n, const RectangleGL &value, const allocator_type &allo)
    : list<RectangleGL> (n, value, allo)
{ }

RectangleGLList::RectangleGLList(iterator first, iterator last, const allocator_type &allo)
    : list<RectangleGL> (first, last, allo)
{ }

RectangleGLList::RectangleGLList(const RectangleGLList &rectList)
    : list<RectangleGL> (rectList)
{ }

//**
//* mergeProjection Ver1.0
//* - Only used when 'this' RectangleGLList contains only one signed direction
//*   (being as part of a condcutorFPList)
//* - If each contains what is after, then erase the latter one.
//* - If only overlapping with the same boundaries, then merge (update *each and erase *after)
//* - Inner loop returns to the first one if any modification happens to the list
//* - does not support sublayer
void RectangleGLList::mergeProjection()
{
    iterator first = this->begin();

    //* Find the first projection
    for ( iterator each=this->begin(); each!=end(); ++each){
        if (each->shapeShift!=0){
            first = each;
            break;
        }
    }

    //* From the first projection, find
    for ( iterator each = first; each != end(); ++each ){
        for ( iterator after = first; after != end(); ++after){

            //* Inner loop starts from the next of each
            if (after==each){
                continue;
            }

            //* Debug: must be in the same direction (regardless of sign)
            if (each->xn*after->xn==0 && each->yn*after->yn==0 && each->zn*after->zn==0){
                cerr << "ERROR: should be in the same direction (regardless of sign)" << endl;
                continue;
            }

            //* Check if *each and *after are on the same surface
            //  and either overlapping or edge neighboring
            if (each->isOverlappingOrEdgeNeighboring(*after)==false){
                continue;
            }

            //* If containing, erase later
            bool flagEraseAfter = false;
            if (each->isContaining(*after)==true){
                flagEraseAfter = true;
            }

            //* If overlapping, not containing, same sign of direction
            else if (each->zn!=0 && each->z1==after->z1){
                //* z-dir
                if ( each->x1==after->x1 && each->x2==after->x2 ){
                    //* same xrange
                    each->y1 = min(each->y1, after->y1);
                    each->y2 = max(each->y2, after->y2);
                    flagEraseAfter = true;
                }
                else if ( each->y1==after->y1 && each->y2==after->y2 ){
                    //* same yrange
                    each->x1 = min(each->x1, after->x1);
                    each->x2 = max(each->x2, after->x2);
                    flagEraseAfter = true;
                }
            }
            else if (each->xn!=0 && each->x1==after->x1){
                //* x-dir
                each->y1 = min(each->y1, after->y1);
                each->y2 = max(each->y2, after->y2);
                flagEraseAfter = true;
            }
            else if (each->yn!=0 && each->y1==after->y1){
                //* y-dir
                each->x1 = min(each->x1, after->x1);
                each->x2 = max(each->x2, after->x2);
                flagEraseAfter = true;
            }

            if (flagEraseAfter==true){
                if (after==first){
                    after = erase(after);
                    first = after;
                }
                else{
                    after = erase(after);
                    after = first;
                }
            }
        }
    }
}


//**
//* mergeProjection Ver1.1
//* - Only used when 'this' RectangleGLList contains only one signed direction
//*   (being as part of a condcutorFPList)
//* - Honor the shapeNormalDistance info
//* - Assume the grid size is 1e-9 (zero)
//* - If each contains what is after, then erase the latter one.
//* - If only overlapping with the same boundaries, then merge (update *each and erase *after)
//* - does not support sublayer
void RectangleGLList::mergeProjection1_1(const float projectionMergeDistance)
{
    const float zero = projectionMergeDistance;
    iterator first = this->begin();

    //* Find the first projection
    for ( iterator each=this->begin(); each!=end(); ++each){
        if (each->shapeShift!=0){
            first = each;
            break;
        }
    }

    //* From the first projection, find
    for ( iterator each = first; each != end(); ++each ){
        if (each->shapeShift==0){
            continue;
        }
        for ( iterator after = first; after != end(); ++after){
            //* Inner loop starts from the next of each
            if (after==each){
                continue;
            }

            if (after->shapeShift==0){
                continue;
            }

            //* Check if *each and *after are on the same surface
            //  and either overlapping or edge neighboring
            if (each->isOverlappingOrEdgeNeighboring(*after)==false){
                continue;
            }

            //* If containing and *after is farther than *each, erase *after later
            bool flagEraseAfter = false;
            if (each->isContaining(*after)==true &&
                each->shapeNormalDistance <= after->shapeNormalDistance ){

                flagEraseAfter = true;
            }

            //* If overlapping but not containing with the same sign of direction
            //* and projected from the same distance
            else if (each->zn * after->zn == 1){
                //* z-dir
                if ( each->x1==after->x1 && each->x2==after->x2 &&
                     abs(each->shapeNormalDistance-after->shapeNormalDistance)<zero ){
                    //* same xrange
                    each->y1 = min(each->y1, after->y1);
                    each->y2 = max(each->y2, after->y2);
                    flagEraseAfter = true;
                }
                else if ( each->y1==after->y1 && each->y2==after->y2 &&
                          abs(each->shapeNormalDistance-after->shapeNormalDistance)<zero ){
                    //* same yrange
                    each->x1 = min(each->x1, after->x1);
                    each->x2 = max(each->x2, after->x2);
                    flagEraseAfter = true;
                }
            }
            else if (each->xn * after->xn == 1){
                //* x-dir
                if ( abs(each->shapeNormalDistance-after->shapeNormalDistance)<zero ){
                    each->y1 = min(each->y1, after->y1);
                    each->y2 = max(each->y2, after->y2);
                    flagEraseAfter = true;
                }
            }
            else if (each->yn * after->yn == 1){
                //* y-dir
                if ( abs(each->shapeNormalDistance-after->shapeNormalDistance)<zero ){
                    each->x1 = min(each->x1, after->x1);
                    each->x2 = max(each->x2, after->x2);
                    flagEraseAfter = true;
                }
            }

            if (flagEraseAfter==true){
                if (after==first){
                    after = erase(after);
                    first = after;
                }
                else{
                    after = erase(after);
                    after = first;
                }
            }
        }
    }
}


//**
//* insertProjectedOverlappingRectangleGL
RectangleGLList::IteratorList RectangleGLList::insertProjectedOverlappingRectangleGL(const RectangleGL &rect, const float distance)
{
    RectangleGLList::IteratorList itList;

    RectangleGLList::iterator dummyEnd = this->insert(this->end(), RectangleGL()); // anchor for inserting projected rect
    for ( RectangleGLList::iterator eachRectIt = this->begin();
          eachRectIt != dummyEnd && eachRectIt->shapeShift==0; ++eachRectIt){

        RectangleGL &eachRect = *eachRectIt;

        if ( !(//* x-dir normal
             ( rect.xn > 0 && eachRect.xn < 0 && rect.x1 < eachRect.x1 && eachRect.x1 < rect.x1+distance ) ||
             ( rect.xn < 0 && eachRect.xn > 0 && eachRect.x1 < rect.x1 && rect.x1 < eachRect.x1+distance ) ||
             //* y-dir normal
             ( rect.yn > 0 && eachRect.yn < 0 && rect.y1 < eachRect.y1 && eachRect.y1 < rect.y1+distance ) ||
             ( rect.yn < 0 && eachRect.yn > 0 && eachRect.y1 < rect.y1 && rect.y1 < eachRect.y1+distance ) ||
             //* z-dir normal
             ( rect.zn > 0 && eachRect.zn < 0 && rect.z1 < eachRect.z1 && eachRect.z1 < rect.z1+distance ) ||
             ( rect.zn < 0 && eachRect.zn > 0 && eachRect.z1 < rect.z1 && rect.z1 < eachRect.z1+distance )
             ) ) {
            continue;
        }

        //* find and insert the projection of rect onto this list
        if (eachRect.isOverlappingProjection(rect)==true){
            RectangleGL projection = eachRect.intersectProjection(rect);
            if (projection.xn!=0){
                projection.shapeNormalDistance = abs(rect.x1-eachRect.x1);
            }
            else if (projection.yn!=0){
                projection.shapeNormalDistance = abs(rect.y1-eachRect.y1);
            }
            else if (projection.zn!=0){
                projection.shapeNormalDistance = abs(rect.z1-eachRect.z1);
            }

            RectangleGLList::iterator it = this->insert( this->end(), projection );

            //* Ver1.0 obsolete
            /*
            if (it->isEmpty()==true){
                //* if the inserted projection is *eachRect itself, delete it
                this->erase(it);
            }
            else{
                //* keep the iterator
                itList.push_back(it);
            }
            */

            //* Ver1.1
            itList.push_back(it);
        }
    }
    this->erase(dummyEnd);
    return itList;
}

void RectangleGLList::absorbCommonSupport()
{
    for ( iterator each = begin(); each!=end(); ++each ){
        iterator after = each;
        for ( ++after; after!=end(); ){
            if (after==each){
                ++after;
                continue;
            }
            if (*after == *each){
                after = this->erase(after);
            }
            else{
                ++after;
            }
        }
    }
}

void RectangleGLList::absorbCommonSupport(RectangleGLList::IteratorList &itList)
{
    if (itList.empty()){
        return;
    }

    for ( IteratorList::iterator eachIt = itList.begin();
          eachIt != --itList.end(); ++eachIt){
        IteratorList::iterator afterIt = eachIt;
        for ( ++afterIt; afterIt!=itList.end(); ){
            RectangleGL &eachRect = **eachIt;
            RectangleGL &afterRect = **afterIt;
            if (eachRect==afterRect){
                *afterIt = this->erase(*afterIt);
                afterIt = itList.erase(afterIt);
            }
            else{
                ++afterIt;
            }
        }
    }
}

//* - Assume single signed direction
//* - Combine all projections.
//* - If the combined projection coincides the underlying rectangle,
//*   remove the farthest projection and combine again.
//* - Repeat until no bad projection.
bool removeOneFarthestBadProjection(RectangleGLList &rectList,
                                    const RectangleGLList::iterator lastSupportIt, const float margin);
void RectangleGLList::removeBadProjection(float margin) {

    //* Find the first projection
    iterator lastSupportIt = this->begin();
    for ( ; lastSupportIt!=end() && lastSupportIt->shapeShift==0; ++lastSupportIt){
    }
    if (lastSupportIt==end()){
        return;
    }
    else{
        --lastSupportIt;
    }

    while( removeOneFarthestBadProjection(*this, lastSupportIt, margin) ){
        ;
    }
}

bool removeOneFarthestBadProjection(RectangleGLList &rectList,
                                    const RectangleGLList::iterator lastSupportIt, const float margin) {

    RectangleGLList::iterator firstProjectionIt = lastSupportIt;
    ++firstProjectionIt;

    //* Make a copy of projections
    RectangleGLList projectionList(firstProjectionIt, rectList.end());

    //* Merge projections
    projectionList.mergeProjection();

    //* Check if any merged projection coincides underlying rectangles
    for (RectangleGLList::iterator projectionIt=projectionList.begin();
         projectionIt != projectionList.end(); ++projectionIt){
        //* each copied projection
        for (RectangleGLList::iterator supportIt=rectList.begin(); supportIt!=firstProjectionIt; ++supportIt) {
            //* each underlying rectangle
            if (projectionIt->isCoincidental(*supportIt, margin)) {
                //* Smell bad projections
                //* Search for the farthest projection component
                float farthestDistance = 0;
                RectangleGLList::iterator farthestProjectionIt;
                for (RectangleGLList::iterator compIt=firstProjectionIt; compIt!=rectList.end(); ++compIt){
                    if (projectionIt->isContaining(*compIt) &&
                        compIt->shapeNormalDistance > farthestDistance ) {

                        farthestDistance = compIt->shapeNormalDistance;
                        farthestProjectionIt = compIt;
                    }
                }
                //* Erase the farthest projection component
                rectList.erase(farthestProjectionIt);
                return true;
            }
        }
    }
    return false;
}


//* Ver1.0 obsolete
//void RectangleGLList::markCommonSupport()
//{
//    for ( iterator each = begin(); each!=--end(); ++each ){
//        if (each->isEmpty()){
//            continue;
//        }
//        iterator after = each;
//        for ( ++after; after!=end(); ++after){
//            if (after->isEmpty()){
//                continue;
//            }
//            if (*after == *each){
//                after->xn = 0;
//                after->yn = 0;
//                after->zn = 0;
//            }
//        }
//    }
//}


//**
//* hasOverlappingRectangle
bool RectangleGLList::hasCommonSupport() const
{
    for ( RectangleGLList::const_iterator each = this->begin();
          each != --this->end(); ++each){
        RectangleGLList::const_iterator after = each;
        for ( ++after; after != this->end(); ++after){
            if (*each == *after){
                return true;
            }
        }
    }
    return false;
}





