#include <iostream>
#include "object.h"

/* run this program using the console pauser or add your own getch, system("pause") or input loop */

int main(int argc, char** argv) {
	
	Point Default;
	Point VRP;
//	Point LRP;
	Vector VPN;
	Vector VUP;
	
	float c;
	
	VRP.SetPoint(5, 5, 30);
	VPN.SetVector(0, 0, -1);
	VUP.SetVector(0, 1, 0);
//	LRP.SetPoint(-10, 10 ,2);
	
	NewCoordinate corView;						// get the new coordinate of camera
	corView.GetNewCoordinate(VUP, VPN);
	
	Matrix RWC;									// Rotation Matrix of World to camera cordinate
	Matrix TWC;									// Translation Matrix of World to camera cordinate
	Matrix MWC;									// Transformation Matrix of World to camera cordiante
	Matrix RCW;									// Roatation Matrix of camera to World cordinate
	Matrix TCW;									// Translation Matrix of camera to World cordinate
	Matrix MCW;									// Transformation Matrix of camera to World cordinate
	
	RWC.SetMatrix(corView.UReturn(), corView.VReturn(), corView.NReturn(), Default);				// set Rotation Matrix
	TWC.SetMatrix(VRP,3);						// set Translation Matrix
	MWC.MultiplProduct(RWC, TWC);
//	MWC.Print();
	
	//Matrix of camera condinate transform to World condinate	
	RCW.InverseRotationMatrix(RWC);				// inverse matrix of RWC
	TCW.InverseTransforMatrix(TWC);				// inverse matrix of TWC
	MCW.MultiplProduct(TCW, RCW);
//	MCW.Print();

	Vector n;
	Obj obj[7];
	
	float ver[4][3] = {	0.0, 0.0, 0.0,		/* v0 */
		 				0.0, 0.0, 20.0,		/* v1 */
						10.0, 0.0, 20.0,	/* v2 */
						10.0, 0.0, 0.0,		/* v3 */
						};
	n.SetVector(0, 1, 0);
	obj[0].SetObj(ver, n, 0.8);						// initialize the Polygon
	
	float ver1[4][3] = { 	0.0, 0.0, 0.0,		/* v0 */
					    	0.0, 10.0, 0.0,		/* v1 */
							0.0, 10.0, 20.0,	/* v2 */
							0.0, 0.0, 20.0,		/* v3 */
							};
	n.SetVector(1, 0, 0);
	obj[1].SetObj(ver1, n, 0.8);					// initialize the Polygon
	
	float ver2[4][3] = {	0.0, 10.0, 0.0,		/* v0 */
		 					0.0, 10.0, 20.0,	/* v1 */
							10.0, 10.0, 20.0,	/* v2 */
							10.0, 10.0, 0.0,	/* v3 */
							};
	n.SetVector(0, -1, 0);
	obj[2].SetObj(ver2, n, 0.8); 					// initialize the Polygon
	
	float ver3[4][3] = {	10.0, 0.0, 0.0,		/* v0 */
		 					10.0, 10.0, 0.0,	/* v1 */
							10.0, 10.0, 20.0,	/* v2 */
							10.0, 0.0, 20.0,	/* v3 */
							};
	n.SetVector(-1, 0, 0);
	obj[3].SetObj(ver3, n, 0.8); 					// initialize the Polygon
	
	float ver4[4][3] = {	0.0, 0.0, 0.0,		/* v0 */
		 					10.0, 0.0, 0.0,		/* v1 */
							10.0, 10.0, 0.0,	/* v2 */
							0.0, 10.0, 0.0,		/* v3 */
							}; 
	n.SetVector(0, 0, 1);
	obj[4].SetObj(ver4, n, 0.8); 					// initialize the Polygon
	obj[4].Poly4::w1 = 0.1;
	obj[4].Poly4::w2 = 0.9;

	obj[5].SetObj(2.0, 2.0, 8.0, 2.0, 0.75); 		// initialize the Sphere
	obj[6].SetObj(8.0, 2.0, 12.0, 2.0, 0.75);
//	obj[5].Sphere::w1 = 0.4;
//	obj[5].Sphere::w2 = 0.6;
//	obj[6].Sphere::w3 = 1.0;

//	float ver5[4][3] = {	0.0, 0.0, 20.0,		/* v0 */
//		 					10.0, 0.0, 20.0,	/* v1 */
//							10.0, 10.0, 20.0,	/* v2 */
//							0.0, 10.0, 20.0,	/* v3 */
//							};
//	n.SetVector(0, 1, 0);
//	obj[7].SetObj(ver, n, 0.8);						// initialize the Polygon

	Ray ray(VRP);								// initialize the Ray
	
	for(int i = 0; i < ROWS; i++)
		for(int j = 0; j < COLS; j++){
			ray.ConstructionOfRay(i, j, MCW);	// Ray construction
			c = ray.RayTracing(obj, 2);			// intersect with two objection
			if(c > 0){
				image[i][j] = ((unsigned char)(int)c);
			}
		}
	

	FILE *outFile = fopen("image.raw", "wb");		// witre the image information into file of raw
    if(0 == outFile)
    {
        perror("ss.raw");
        return -1;
    }
/*
    for(int row=0; row<ROWS; ++row)
        for(int col=0; col<COLS; ++col)
            fwrite(image[row]+col, sizeof(unsigned char), 1, outFile);
*/	
	fwrite(image, COLS, ROWS, outFile);
    fclose(outFile);

	return 0;
}