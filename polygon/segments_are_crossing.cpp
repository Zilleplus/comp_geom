/*
 * simple algorithms to find if 2 segments are crossing.
 * A segment is defines as a straight line with a start and end position
 */

#include<iostream>
#include<tuple>
#include<string>
#include<optional>
#include<limits>
#include<assert.h>

struct Point{
    double x;
    double y;
    Point(double x, double y):x(x),y(y){}
};
struct Segment{
    Point start;
    Point end;
    Segment(Point start, Point end):start(start),end(end){}
};
std::string toStringPoint(Point p)
{
    return "["+std::to_string(p.x)+";"+std::to_string(p.y)+"]";
}

std::string toStringSegment(Segment s)
{
    return "(start="+toStringPoint(s.start) + "; end=" + toStringPoint(s.end) + ")";
}

bool isCrossingInOnePoint(
        Segment first,
        Segment second)
{
    // Check if the y axis is in range:
    if(std::min(first.start.y, first.end.y) <= std::max(second.start.y, second.end.y)
       && std::max(first.start.y, first.end.y) >= std::min(second.start.y, second.end.y))
    {
        // Check if the y axis is in range:
        if(std::min(first.start.x, first.end.x) <= std::max(second.start.x, second.end.x)
           && std::max(first.start.x, first.end.x) >= std::min(second.start.x, second.end.x))
        {
            return true;
        }
    }
    
    return false;
}

std::optional<Point> crossingPoint(Segment first, Segment second)
{
    // makes sure the ranges are fine, and make sure they are not parallel.
    if(!isCrossingInOnePoint(first,second))
    {
        return {};
    }
    // If we know that there is one unique point we can define the problem as:
    // -> forall u v in range[0;1]
    // (first.end-first.start)*u + (first.start) = (second.end-second.start)*v + (second.start)
    // put it to zero and take the gradient:
    // -> set: fdiff = first.end-first.start
    // -> set: sdiff = second.end-second.start
    // [fdiff,-sdiff]*[u;v] = second.start - first.start
    
    // solve the equation
    double fdiffX = first.end.x - first.start.x;
    double fdiffY = first.end.y - first.start.y;
    double sdiffX = second.end.x - second.start.x;
    double sdiffY = second.end.y - second.start.y;
    
    // if the dim(x)=2 
    // if A = [[a;c],[b;d]]
    // det(A)= a*b - b*c
    double a = fdiffX;
    double c = fdiffY;
    double b = -sdiffX;
    double d = -sdiffY;
    // invers matrix  = 1/abs(det(A)) * [[d;-c],[-b;a]]
    std::cout << "[" << std::to_string(a) << " " << std::to_string(b) << "]" << std::endl;
    std::cout << "[" << std::to_string(c) << " " << std::to_string(d) << "]" << std::endl;
    std::cout << std::to_string(a*d) << "..." << std::to_string(-b*c) << std::endl;
    double detA = a*d-b*c;
    if(std::abs(detA) < std::numeric_limits<double>::epsilon())
    {
        std::cout 
            << "error det=0 should not be possible due to the check at the start of this method" 
            << std::endl;
        return {};
    }
    double detA_1 = 1.0/std::abs(detA);

    double rhsX = -(second.start.x - first.start.x);
    double rhsY = second.start.y - first.start.y;

    double xCross = detA_1*(a*rhsX + -b*rhsY);
    double yCross = detA_1*(-c*rhsX + d*rhsY);

    return Point(xCross,yCross);
}

void demoSegmentsCrossing()
{
    std::cout << "Segments that do cross:" << std::endl;
    Segment a1(Point(0,0),Point(1,1));
    Segment a2(Point(0,1),Point(1,0));
    std::cout << "a1=" + toStringSegment(a1) << std::endl;
    std::cout << "a2=" + toStringSegment(a2) << std::endl;
    std::cout 
        <<" -> should cross on [0.5;0.5] isCrossingInOnePoint=" 
        << std::to_string(isCrossingInOnePoint(a1,a2)) 
        << std::endl;
    std::cout << "crossPoint=" << toStringPoint((*crossingPoint(a1, a2))) << std::endl;
}

// This is a version using determinant, can also be done using the cross product.
// But then we need a sine, and a 3th dimension.
bool segmentsAreParallel(Segment s1, Segment s2)
{
    // If the determinant of the matrix with the two direction vectors is zero, 
    // then lines are parrallel.
    // As this means that they are not dependent.
    // They might also be collinear.
    // Should we normalize them first? Could be more stable.
    double directionS1x = s1.end.x -s1.start.x;
    double directionS1y = s1.end.y -s1.start.y;

    double directionS2x = s2.end.x -s2.start.x;
    double directionS2y = s2.end.y -s2.start.y;

    // Might be interesting to put this inside a seperate method next time.
    double a = directionS1x;
    double b = directionS2x;
    double c = directionS1y;
    double d = directionS2y;
    double det = a*d-b*c;

    // if det < 0 -> they are parallel.
    return std::abs(det) < std::numeric_limits<double>::epsilon();
}

double norm2(Point p)
{
    return p.x*p.x + p.y*p.y;
}

Point scale(Point p, double scale)
{
    return Point(p.x * scale, p.y * scale);
}

Point normalize(Point p)
{
    auto vecSize = norm2(p);
    if(vecSize < std::numeric_limits<double>::epsilon())
    {
        return p;
    }
    auto res = scale(p, 1.0/vecSize);
    return res;
}

bool segmentsAreCollinear(Segment s1, Segment s2)
{
    if(!segmentsAreParallel(s1,s2)) // direction is parallel
    {return false;}

    auto offset1 = normalize(s1.start);
    auto offset2 = normalize(s2.start);

    return std::abs(offset1.x - offset2.x) < std::numeric_limits<double>::epsilon() 
        && std::abs(offset1.y - offset2.y) < std::numeric_limits<double>::epsilon();
}

void demoAreCollinear()
{
    Segment a1(Point(0,1),Point(1,1));
    Segment a2(Point(0,0),Point(1,0)); // a1 and a2 are parallel but not colliner
    Segment a3(Point(1,0),Point(2,0)); // a2 and a3 are not parallel

    std::cout 
        << "segment:" << toStringSegment(a1) << std::endl
        << "and segment:" << toStringSegment(a2) << std::endl
        << "should not be collinear=> " 
        << std::to_string(segmentsAreCollinear(a1,a2))
        << std::endl;
    std::cout 
        << "segment:" << toStringSegment(a2) << std::endl
        << "and segment:" << toStringSegment(a3) << std::endl
        << "should be collinear= " 
        << std::to_string(segmentsAreCollinear(a2,a3))
            << std::endl;
}

void demoAreParallel()
{
    Segment a1(Point(0,1),Point(1,1));
    Segment a2(Point(0,0),Point(1,0)); // a1 and a2 are parallel
    Segment a3(Point(0,1),Point(1,0)); // a2 and a3 are not parallel

    std::cout 
        << "should be parallel= " 
        << std::to_string(segmentsAreParallel(a1,a2))
        << std::endl;
    std::cout 
        << "should not be parallel= " 
        << std::to_string(segmentsAreParallel(a2,a3))
            << std::endl;
}

void printWall(std::string name = "")
{
    for(int i=0;i<20;i++)
    {
        std::cout << "#";
    }

    std::cout << name ;

    for(int i=0;i<20;i++)
    {
        std::cout << "#";
    }

    std::cout << std::endl;
}

int main()
{
    printWall(" demo Collinear ");
    demoAreCollinear();
    printWall(" demo parallel ");
    demoAreParallel();
    printWall(" demo crossing ");
    demoSegmentsCrossing();
}
