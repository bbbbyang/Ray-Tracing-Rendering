#ifndef OBJECT_H
#define OBJECT_H

#include <math.h>
#include <iomanip>
using namespace std;

#define ROWS 512
#define COLS 512
unsigned char image[ROWS][COLS];
float xmin = 0.0175;
float ymin = -0.0175;
float xmax = -0.0175;
float ymax = 0.0175;
float focal = 0.05;
float Ip = 200.0;

class Matrix;				//Forward declaration
class Ray;					//Forward declaration

class Point{
	public:
		float point[4];
	public:
		//constructor
		Point(){
			for(int i = 0; i < 3; i++)
				point[i] = 0;
			point[3] = 1;
		}
		//constructor
		Point(float a, float b, float c){
			point[0] = a;
			point[1] = b;
			point[2] = c;
			point[3] = 1;
		}
		//destructor
		~Point(){
		}
		//set the value of a point
		void SetPoint(float a, float b, float c){
			point[0] = a;
			point[1] = b;
			point[2] = c;
		}
		//print the properties of the point
		void Print(){
			cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<endl;
			cout<<endl;
		}
		//get the coordinate value of a point in corresponding coordinate
		void TranslatePoint(Matrix&);
};

Point LRP(5, 9, 8);					//the position of light

class Vector{
	public:
		float vector[4];
	public:
		//constructor
		Vector(){
			for(int i = 0; i < 4; i++)
				vector[i] = 0;
		}
		//destructor
		~Vector(){
		}
		//set the value of a vector
		void SetVector(float a, float b, float c)
		{
			vector[0] = a;
			vector[1] = b;
			vector[2] = c;
		}
		//unitizale the vector
		void UnitizationVector(){
			float squaremode = vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2];
			squaremode = sqrt(squaremode);
			vector[0] = vector[0]/squaremode;
			vector[1] = vector[1]/squaremode;
			vector[2] = vector[2]/squaremode;
		}
		//cross puduct of two vectors
		void CrossProduct(Vector V1, Vector V2){
			float x = V1.vector[1]*V2.vector[2] - V1.vector[2]*V2.vector[1];
			float y = V1.vector[2]*V2.vector[0] - V1.vector[0]*V2.vector[2];
			float z = V1.vector[0]*V2.vector[1] - V1.vector[1]*V2.vector[0];
			this->SetVector(x, y, z);
		}
		//product with other vector
		float ProductVector(Vector V1){
			float ref = 0;
			for(int i = 0; i < 3; i++)
				ref += vector[i] * V1.vector[i];
			return ref;
		}
		//mutiple the vector
		Vector Large(float ref){
			Vector newvector;
			for(int i = 0; i < 3; i++)
				newvector.vector[i] = ref * vector[i];
			return newvector;
		}
		//vector of two point subtraction
		void MinPoint(Point P1, Point P2){
			for(int i = 0; i < 3; i++)
				vector[i] = P1.point[i] - P2.point[i];
		}
		//vector of two vectors subtraction
		void MinVector(Vector P1, Vector P2){
			for(int i = 0; i < 3; i++)
				vector[i] = P1.vector[i] - P2.vector[i];
		}
		//print the properties of the vector
		void Print(){
			cout<<vector[0]<<" "<<vector[1]<<" "<<vector[2]<<endl;
			cout<<endl;
		}
};

class Matrix{
	public:
		float mar[4][4];                 //(N+1)th-order Matrix
	public:
		//constructor, default value of element in Matrix
		Matrix(){
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++)
				{
					mar[i][j] = 0;
					if( i == j)
						mar[i][j] = 1;
				}
		}
		//destructor
		~Matrix(){
		}
		//set the value of matrix
		void SetMatrix(Vector V1, Vector V2, Vector V3, Point P){
			
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++)
					switch (i){
						case 0: mar[0][j] = V1.vector[j];
						case 1: mar[1][j] = V2.vector[j];
						case 2: mar[2][j] = V3.vector[j];
						case 3: mar[3][j] = P.point[j];
					}
		}
		//set the kth colunm of the Matrix
		void SetMatrix(Point P, int k){
			Point V;
			V.SetPoint(-P.point[0], -P.point[1], -P.point[2]);
			for(int i = 0; i < 4; i++)
				mar[i][k] = V.point[i];
		}
		//inverse the Translation Matrix
		void InverseTransforMatrix(Matrix m){
			for(int i = 0; i < 3; i++)
				mar[i][3] = -m.mar[i][3];
		}
		//inverse the Rotation Matrix
		void InverseRotationMatrix(Matrix m){
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++)
					mar[i][j] = m.mar[j][i];
		}
		//mutiply of two matrixs
		void MultiplProduct(Matrix m, Matrix n){
			float value;
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++){
					value = 0;
					for(int k = 0; k < 4; k++){
						value += m.mar[i][k] * n.mar[k][j];
					}
					mar[i][j] = value;
				}
		}
		//print the Matrix
		void Print(){
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++)
				{
					cout<<setw(10);         //set the print width of the Matrix element
					cout<<setiosflags(ios::fixed)<<setprecision(5)<<mar[i][j]<<" ";
					if( j == 3)
						cout<<endl<<endl;
				}
		}
};

//get the coordinate value of a point in corresponding coordinate, product with a translate matrix
void Point::TranslatePoint(Matrix& M){
			float value;
			Point V;
			for(int i = 0; i < 4; i++){
					value = 0;
					for(int k = 0; k < 4; k++)
						value += M.mar[i][k] * point[k];
					V.point[i] = value;
				}
				for(int i = 0; i < 3; i++)
					point[i] = V.point[i];
			}


class NewCoordinate{
	protected:
		Vector n;                   //n unit vector
		Vector u;                   //u unit vector
		Vector v;                   //v unit vector
	public:
		//constructor
		NewCoordinate(){
		}
		//destructor
		~NewCoordinate(){
		}
		//get new coordinate with two vectors
		void GetNewCoordinate(Vector V1, Vector V2){
			//get the vector n
			n = V2;
			n.UnitizationVector();
	
			//corss product of V1 and V2
			//get the vector u and unitlize it
			u.CrossProduct(V1, V2);
			u.UnitizationVector();
			
			//corss product of n and u	
			//get the unit vector of v
			v.CrossProduct(n, u);
		}
		//return the unit vector along the axis of the new coordinate
		Vector NReturn(){
			return n;
		}
		Vector UReturn(){
			return u;
		}
		Vector VReturn(){
			return v;
		}
};

class Sphere{
	public:
		Point center;		//the center of the sphere
		float radius;		//the radius of the sphere
		float kd;
		float w1;			//percentage of contribution of shading value
		float w2;			//percentage of contribution of reflection value
		float w3;			//percentage of contribution of refraction value
	public:
		//constructor
		Sphere(){
			w1 = 1.0;
			w2 = 0.0;
			w3 = 0.0;
		}
		//set sphere value
		void SetSphere(float a, float b, float c, float d, float e){
			center.SetPoint(a, b, c);
			radius = d;
			kd = e;
		}
		//destructor
		~Sphere(){
		}
		//intersection with a ray
		float RaySphereIntersection(Ray* const);
		//print the sphere
		void Print(){
			center.Print();
			cout<<radius<<" "<<kd<<endl;
		}
};

class Poly4{
	public:
		Point vertex[4];							//four vertex points
		Vector normal;								//normal vector of the polygon
		float kd;
		float w1;									//percentage of contribution of shading value
		float w2;									//percentage of contribution of reflection value
		float w3;									//percentage of contribution of refraction value
	public:
		//constructor
		Poly4(){
			w1 = 1.0;
			w2 = 0.0;
			w3 = 0.0;
		}
		//set polygon value
		void SetPoly4(float A[4][3], Vector V, float k){
			for(int i = 0; i < 4; i++)
				vertex[i].SetPoint(A[i][0], A[i][1], A[i][2]);
			//initialize the normal vector
			normal = V;
			kd = k;
			//normal[3] equals the the D in the plane function (Ax + By + Cz + D = 0)
			normal.vector[3] = -normal.vector[0] * vertex[0].point[0] - normal.vector[1] * vertex[0].point[1] - normal.vector[2] * vertex[0].point[2];
		}
		//destructor
		~Poly4(){
		}
		//project the polygon to the plane, return the lagest number of element in normal vector
		int Projection()
		{
			int large = 0;
			int count = 0;
			//get the largest value of the normal vector
			large = max(abs(normal.vector[0]), max(abs(normal.vector[1]), abs(normal.vector[2])));
			//get the number of largest value
			for(int i = 0; i < 3; i++)
				if(abs(normal.vector[i]) == large)
					count = i;
			return count;
		}
		//intersection with a Ray
		float RayPolyIntersection(Ray* const);
		//print the polygon
		void Print(){
			vertex[0].Print();
			vertex[1].Print();
			vertex[2].Print();
			vertex[3].Print();
			normal.Print();
			cout<<kd<<endl;;
		}
};

//object contains sphere and polygon
class Obj:public Sphere, public Poly4{
	public:
		int flag;					//distinct which is sphere or polygon
	public:
		void SetObj(float a, float b, float c, float d, float e){
			SetSphere(a, b, c, d, e);
			flag = 1; 				//sphere
		}
		void SetObj(float A[4][3], Vector V, float k){
			SetPoly4(A, V, k);
			flag = 2;				//polygon
		}
		float RayObjIntersection(Ray* const);
};
	
class Ray{
	public:
		Point P;							//VRP
		Vector V;							//the unit vector from VRP to the point in screen
	public:
		//constructor
		Ray(Point O){
			P = O;							//initialize the VRP
		}
		//destructor
		~Ray(){
		}
		//Ray construction
		void ConstructionOfRay(int i, int j, Matrix M){
			//transform the point [i, j] in the image plane to the [m, n] in the screen of camera
			float m;
			float n;
			m = xmin;
			n = ymax;
			n -= i*0.035/511;
			m -= j*0.035/511;
			Point P1;						//the point in the screen
			P1.SetPoint(m, n, focal);
			P1.TranslatePoint(M);			//transform the P1 from camrea coordianate to the world coordinate
			V.MinPoint(P1, P);				//get the vector from two point, P1 and VRP
			V.UnitizationVector();			//unitize the vector
		}
		//ray intersect with the object
		float ObjectIntersection(Obj obj){
			return obj.RayObjIntersection(this);
		}
		//get the shading value of the intersection point
		float Shading(Point t, Vector Nor, float kd){
			Vector L;
			float ref;
			//the direction of the L is from intersection point in sphere to the LRP 
			L.SetVector(LRP.point[0] - t.point[0], LRP.point[1] - t.point[1], LRP.point[2] - t.point[2]);
			L.UnitizationVector();			//unitilze the vector
			ref = L.ProductVector(Nor);		//product with normail vector
			ref = ref * Ip * kd + 50;				//get the shading value
			return ref;
		}
		//compute the shadow
		float Shadow(Point t, int pointer, Obj obj[]){
			Ray ShadowRay(LRP);							//shadowray emit from LRP to the intersection point
			ShadowRay.V.MinPoint(t, LRP);
			ShadowRay.V.UnitizationVector();
			int counter = 0;
			float shadowref;							//shadowray intersect with other objects
			float lengthref;							//the distance between LRP and intersection point
			
			lengthref = ShadowRay.ObjectIntersection(obj[pointer]);
			
			for(int i = 0; i < 7; i++){
				if(i != pointer){
					shadowref = ShadowRay.ObjectIntersection(obj[i]);
					if((shadowref < lengthref)&&(shadowref > 0))	//shadowray has intersection point with other objects
						counter++;									//and the distance is smaller than lengthref
				}		
			}
			if(counter != 0)
				return 50;
				else
					return 0;
		}
/*		Ray RayRefraction(Point t, Vector nor, int level){
			Ray rayrefraction(t);
			Vector L;
			Vector normal;
			float angle;
			
			L = V.Large(-1);
			angle = nor.ProductVector(L);
			angle = sqrt(1 - pow(angle, 2));
			normal.CrossProduct(L, nor);
			normal.vector[3] = -normal.vector[0] * t.point[0] - normal.vector[1] * t.point[1] - normal.vector[2] * t.point[2];
		}*/		//under construction
		Ray RayReflection(Point t, Vector normalvector){
			Ray rayreflection(t);			//reflection ray start from intersection point t
			Vector L;						//the opposite direction vector of incident ray
			Vector ref;
			float angle;					//cosine value of the angle between L and N
			
			L = V.Large(-1);
			angle = normalvector.ProductVector(L) * 2;		//2 * cos(N * L)
			ref = normalvector.Large(angle);				//2 * cos(N * L) * N
			
			rayreflection.V.MinVector(ref, L);				//2 * cos(N * L) * N - l
			rayreflection.V.UnitizationVector();
			
			return rayreflection;
		}
		//Ray tracing
		float RayTracing(Obj obj[], int level){
			float re[7];
			float kd;
			int pointer;					//record the number of the shading object
			float min;
			float x;						//intersection point x value
			float y;						//intersection point y value
			float z;						//intersection point z value
			float w1;
			float w2;
			float w3;
			Vector nor;
			Point t;						//intersection point with object
			if(level < 1)
				return 0;
				else{
					for(int i = 0; i < 7; i++)
						re[i] = ObjectIntersection(obj[i]);
					//get the distance from the point to intersection point
					min = 0;
					pointer = 0;
					
					//give min a value of any non-zone distance from re[]
					for(int i = 0; i < 7; i++){
						if(re[i] > 0)
						min = re[i];
					}
					//get the smallest non-zone distance from re[]
					for(int i = 0; i < 7; i++){
						if((min >= re[i])&&(re[i] > 0)){
							min = re[i];
							pointer = i;
						}
					}
					if(0 == min)
						return 0;
							else{
								x = P.point[0] + re[pointer] * V.vector[0];
								y = P.point[1] + re[pointer] * V.vector[1];
								z = P.point[2] + re[pointer] * V.vector[2];
								t.SetPoint(x, y ,z);							//intersection point t
								if(obj[pointer].flag == 1){						//flag == 1 the objec is sphere
									nor.MinPoint(t, obj[pointer].center);		//get normal vector
									nor.UnitizationVector();
									//get the parameter of the object
									kd = obj[pointer].Sphere::kd;
									w1 = obj[pointer].Sphere::w1;
									w2 = obj[pointer].Sphere::w2;
									w3 = obj[pointer].Sphere::w3;
								}
									else{										//the object is polygon
										nor = obj[pointer].Poly4::normal;
										kd = obj[pointer].Poly4::kd;
										w1 = obj[pointer].Poly4::w1;
										w2 = obj[pointer].Poly4::w2;
										w3 = obj[pointer].Poly4::w3;
									}
								float shadowvalue = Shadow(t, pointer, obj);	//shaow function
								
								if(w1 > 0)										//shading condition
									x = Shading(t, nor, kd);
								if((w2 > 0)&&(level > 1)){						//reflection condition
									Ray rayreflection = RayReflection(t, nor);
									y = rayreflection.RayTracing(obj, level - 1);
								}
								
				/*				if((w3 > 0 )&&(level > 1)){						//refraction condition
									Ray rayrefraction = RayRefraction(t, nor , level);
									z = rayrefraction.RayTracing(obj, level - 1);
								}*/
								return ((x *w1 + y * w2 + z * w3)- shadowvalue);
							}
				}
		}
};

float Sphere::RaySphereIntersection(Ray*const ray){
	//intersect with sphere Ax^2 + By^2 +Cz^2 = R^2;  A = 1
	float B = 0;
	float C = 0;
	float dis;
	for(int i = 0; i < 3; i++){
		B += ray->V.vector[i] * (ray->P.point[i] - center.point[i]);
		C += pow((ray->P.point[i] - center.point[i]), 2);
	}
	B = B * 2;											//get value of B
	C = C - pow(radius, 2);								//get value of C
	
	dis = pow(B, 2) - 4 * C;							//delta of the intesction equation
	if(dis < 0)
		return 0;										//delta < 0 no intersction
		else{
			float t1;
			t1 = (-B - pow(dis, 0.5))/2;				//the smaller root t1
			if(t1 >= 0)
				return t1;
				else
					return t1 = (- B + pow(dis, 0.5))/2;//if samller root is negetive, then calculate the bigger one
			}
}

float Poly4::RayPolyIntersection(Ray* const ray){
	int dis = Projection();								//projection plane
	float projectionvertex[5][2];						//projection vertex, 5th point equal the first point
	int count = 0;										//count the positive number of root
	float root = 0;
	float M = 0;
	float N = 0;
	float t2;
	float x;
	float y;
	float z;
	float k;
	Point test;											//intersection point
	//intersect with plane Ax + By + Cz + D = 0
	for(int i = 0; i <3; i++){
		M -= normal.vector[i] * ray->P.point[i];
		N += normal.vector[i] * ray->V.vector[i];
	}
	M -= normal.vector[3];
	if(0 == N)
		return 0;
		else{
			t2 = M / N;									//get the intersection distance
			x = ray->P.point[0] + ray->V.vector[0] * t2;
			y = ray->P.point[1] + ray->V.vector[1] * t2;
			z = ray->P.point[2] + ray->V.vector[2] * t2;
			test.SetPoint(x, y, z);						//intersection point
			//projection with disfferent cases to (U, V)
			//move the intersection point to the origin
			for(int i = 0; i < 4; i++)
				switch(dis){
					case 0:
						projectionvertex[i][0] = vertex[i].point[1] - test.point[1];
						projectionvertex[i][1] = vertex[i].point[2] - test.point[2];
						break;
					case 1:
						projectionvertex[i][0] = vertex[i].point[0] - test.point[0];
						projectionvertex[i][1] = vertex[i].point[2] - test.point[2];
						break;
					case 2:
						projectionvertex[i][0] = vertex[i].point[0] - test.point[0];
						projectionvertex[i][1] = vertex[i].point[1] - test.point[1];
						break;
				}
			//5th point equal the first point
			projectionvertex[4][0] = projectionvertex[0][0];
			projectionvertex[4][1] = projectionvertex[0][1];
			
			//the ith and (i+1)th point construct a line
			for(int i = 0; i < 4; i++){
				if(0 == (projectionvertex[i][0] - projectionvertex[i+1][0]))		//the line parallel to V, with no root in V axis
					root = 0;
					else{
						k = (projectionvertex[i][1] - projectionvertex[i+1][1]) / (projectionvertex[i][0] - projectionvertex[i+1][0]);
						//the two vertex of the line is on the two sides of the V axis get root
						if((projectionvertex[i][0] * projectionvertex[i+1][0]) < 0)
							root = projectionvertex[i][1] - k * projectionvertex[i][0];
							else
								root = 0;
					}
				//if root is positive, count++
				if(root > 0)
					count++;
			}
			//if count is even, the intersection point is outside
			if(0 == count%2)
				return 0;
				else					//else the intersection point is inside
					return t2;
		}
}

float Obj::RayObjIntersection(Ray* const ray)
{
	if(1 == flag)
		return Sphere::RaySphereIntersection(ray);
		else
			return Poly4::RayPolyIntersection(ray);
}
#endif
