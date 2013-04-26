// TestTrace.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include "Medium.h"
#include "PluginFactories.h"
#include "vtkCubeSource.h"
#include "vtkCylinderSource.h"
#include "vtkPolyData.h"
#include <vtkPointData.h>
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkImageWriter.h"
#include "Def.h"
#include "TestTrace.h"
#include "vector.h"
#include <vtkNew.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkPolyDataWriter.h>
#include "vtkSphereSource.h"
#include <fstream>
#include <math.h>
#include <vtkSmartPointer.h>
#include <deque>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkSmartPointer.h>
#include <vtkSTLWriter.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangleFilter.h>

using namespace std;
using namespace RayTracer;


vtkPolyData* CreateCube(char *inp2)
{
	std::string inputFilename = inp2;
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(inputFilename.c_str());
	reader->Update();
	vtkPolyData * pd = reader->GetOutput();

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInput(pd);
	
	
	vtkSmartPointer<vtkPolyDataNormals> pdn = vtkSmartPointer<vtkPolyDataNormals>::New();
	pdn->SetInput(triangleFilter->GetOutput());
	pd = pdn->GetOutput();
	pdn->Update();
	pd->Register(0);
	return pd;

	/*vtkCubeSource* cube = vtkCubeSource::New();
	cube->SetCenter(0,0,0);
	cube->SetBounds(-25,25,-25,25,-25,25);
	cube->Update();*/
	/*vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->RotateY(45);
	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformFilter->SetTransform(transform);
	transformFilter->SetInputConnection(cube->GetOutputPort());
	transformFilter->Update();
	cube->Update();*/
	/*vtkPolyDataWriter* wr = vtkPolyDataWriter::New();
	wr->SetInput(cube->GetOutput());	
	wr->SetFileName("E:\\cube.vtk");
	wr->Update();
	wr->Write();*/
	//return cube->GetOutput();
	
}

double Directivity(vector3d incomingRaydir) 
{
	return 1;
	vector3d mainLobeDirection = math::vector3d(-0.65,1,1);
	double v = (dot(unit(mainLobeDirection),unit(incomingRaydir)));
	double angle = acos(v);
	double pi = 3.14;
	double p= angle*180/pi;
	cout<<p<<"  ";
	if(p >71)
		return 0;
	return angle;; 
}

vector3d GetReflectedRayDirection(vector3d incident, vector3d normal)
{
	//IMPLEMENTATION OF PHONG SHADING MODEL//			
	double reflet = 2.0 * dot(unit(incident),unit(normal));
	vector3d reflected = unit(incident - reflet*normal);
	//reflected ray vector r = i - 2(i*n)n//
	return reflected;
	//"E:\Aman\TestTrace\TestTrace\centered cube 1 edge-center 50mm side, 40x10mm triangular notch.stl"
}

void ExecuteRayTracing(char *inp, char* inp2)
{
	//Create a PolyData
	vtkPolyData* cube = CreateCube(inp2);	
	vtkPolyDataWriter *wr = vtkPolyDataWriter::New();
	wr->SetInput(cube);
	wr->SetFileName("cube.vtk");
	wr->Update();


	//A Scan
	const int ASCAN_LENGTH = 1000; //us
	double sum_array[ASCAN_LENGTH];
	double aScanArray[ASCAN_LENGTH];	

	//Source
	transducer transmitter;
	const int PULSE_WIDTH = 5; //us
	const vector3d POSITION = math::vector3d(0,0,24);
	transmitter.position = POSITION;//Show the cases for differen transmitter pos
	double spherical_pulse[PULSE_WIDTH];
	double a = 1;
	double pi = 3.14;
	std::fill(spherical_pulse,spherical_pulse+PULSE_WIDTH,0);
	double velocity = 5; // mm/us

	//MEDIUM
	const char * library = "RayTracer";
	const char * type = "PolyDataMedium";
	boost::shared_ptr<RayTracer::Medium> medium;
	medium.reset(RayTracer::CreateMedium("BulletRayTracer", "BulletMedium", 5, cube));

	//New VTK Image 
	vtkImageData *image = vtkImageData::New();
	const int IMAGE_WIDTH = 20;
	const int IMAGE_HEIGHT= 20;
	image->SetDimensions(IMAGE_WIDTH,IMAGE_HEIGHT,ASCAN_LENGTH);
	image->SetSpacing(1.0,1.0,1.0);
	image->SetExtent(-10,10,-10,10,-10,990);
	image->SetOrigin(0,0,0);
	image->SetScalarTypeToDouble();
	image->SetNumberOfScalarComponents(1);
	image->AllocateScalars();
	image->Update();
	image->GetPointData()->GetScalars()->FillComponent(0,0);

	double op[PULSE_WIDTH];
	vector3d end = math::vector3d(0,0,-5000);
	const double CAMERA_POSITION = 24;
	vector3d intersection;
	vector3d normal;
	double dist2;
	material mesh = {mesh.diffuse = 0.7, mesh.power = 60, mesh.reflectance = 0.99, mesh.specular =0.6};
	double count = 0;
	ofstream ouput;
	ouput.open("E:\\ouput.txt");
	
	//Scanning
	for(int y = 0; y < IMAGE_HEIGHT ;y++)
	{
		for (int x = 0; x < IMAGE_WIDTH; x++)
		{	
			std::fill(sum_array, sum_array + ASCAN_LENGTH, 0);
			ouput<<"NEXT RAY BEGINS HERE  "<<x<<" "<<y<<endl;
			double coef = 1.0;	int level = 0; double pathLength = 0; double time = 0;
			/*We shoot a ray(startingpoint,direction) in the z direction
			& then we look for it's intersection with the medium*/
			vector3d startPos = math::vector3d(double(x),double(y),CAMERA_POSITION); //Camera Plane
			ray viewRay = {startPos,end};			
			vector3d viewRayDir = unit(viewRay.endPoint - viewRay.startPoint);
			do
			{
				std::fill(aScanArray, aScanArray + ASCAN_LENGTH, 0);
				std::fill(op, op + PULSE_WIDTH, 0);
				std::fill(transmitter.intensity, transmitter.intensity + PULSE_WIDTH, 0);
				//Look for closest Intersection
				if(!medium->Intersect(viewRay.startPoint,viewRay.endPoint,intersection,normal,dist2))
					break;
				//Find the intersection point & normal at the point
				vector3d newStart = intersection;				
				unit(normal);
				//Calculate Illumination
				//Shoot from intersection to source
				vector3d newViewRayDir = unit(transmitter.position - newStart);
				double t = (dot(newViewRayDir,newViewRayDir));
				if(t==0) //discard rays that hit directly from camera to transmitter
					break;
				if(fabs(dot(newViewRayDir,normal)) >= 0.0001 && t >= 0.0001) 
				{
					ray newViewRay;
					newViewRay.startPoint = newStart;				
					newViewRay.endPoint = newStart +  newViewRayDir*1000;
					//computation of shadows
					bool inShadow = false; 
					vector3d points,nml;
					double dis;
					if(medium->Intersect(newViewRay.startPoint,transmitter.position,points,nml,dis))
						inShadow = true; //this will not contribute to the illumination
					if (!inShadow) 
					{
						double directivityCoeff = Directivity(newViewRayDir);
						double len1 = sqrt(dist2);
						double len2 = sqrt(dis);
						pathLength += len1+len2;
						time = pathLength/velocity;
						//double numberOfElements = 10;
						double step_length = 360/PULSE_WIDTH;//numberOfElements;
						double drop = std::exp(-pathLength/100);
						//ouput<<"drop  "<<drop;
						for(int k = 0; k < PULSE_WIDTH; k++)
						{
							spherical_pulse[k] = a*sin(step_length*k*pi/180);
							transmitter.intensity[k] = spherical_pulse[k]*drop; //transmitter int amplitude is gauss pulse of 50 mic-sec
							//ouput<<"angle  "<<step_length*k<<"  "<<"ampl  "<<transmitter.intensity[k]<<endl;
						}
						count++;
						double lambert = fabs(dot(newViewRayDir,normal) * (coef));
						// lambert for diffuse component & reflection					
						//IMPLEMENTATION OF PHONG SHADING MODEL//			
						vector3d phongdir = GetReflectedRayDirection(newViewRayDir, normal);//phong Direction
						ouput<<phongdir<<endl;
						double dotl = (dot(phongdir,viewRayDir));
						double phong = std::max(dotl,0.0); //phong dependency on viewer direction
						phong =  double(std::powf(phong,mesh.power));
   						/*specular term affect is measured by the alpha power of cosine of angle 
						between viewer direction & reflection vector */
						double amplitude = directivityCoeff*lambert * mesh.diffuse;
						//ouput<<"phong "<<phong<<endl;
						double amplitude2 = directivityCoeff*phong*mesh.specular;
						std::fill(aScanArray, aScanArray + ASCAN_LENGTH, 0);
						for(int j = 0 ; j < PULSE_WIDTH; j ++)
						{
							op[j]+= transmitter.intensity[j] * amplitude;
							op[j]+= transmitter.intensity[j] *amplitude2;
							aScanArray[int(pathLength +j)] = op[j];
							op[j] = 0;
						}
						for(int k = 0; k < ASCAN_LENGTH; k++)
						{
							sum_array[k] += aScanArray[k];
						}
						ouput<<"time  "<<time<<"pathLength   "<<pathLength<<"level  "<<level<<" coef  "<<coef<<endl;
					}
				} 
				coef *= mesh.reflectance;
				//Run the doWhile loop for subsequent reflections...
				viewRay.startPoint = newStart; //intersection point
				viewRayDir = GetReflectedRayDirection(viewRayDir, normal);
				viewRay.endPoint = viewRay.startPoint + 2000*viewRayDir;
				level++;
			}
			while ((coef>0.0) &&(level<10));
			//Creating Image
			for(int j = 0; j < ASCAN_LENGTH; j++)
			{
				image->SetScalarComponentFromDouble(x-10,y-10,j-10,0,sum_array[j]);
				image->Update();
			}			
			for(int k = 0;k<ASCAN_LENGTH;k++)
			{
				ouput<<"scan  "<<sum_array[k]<<endl;
			}
		}
	}
	cout<<"count  "<<count;
	//Writing to Image 
	vtkNew<vtkStructuredPointsWriter> writer;
	writer->SetFileName(inp);
	writer->SetInput(image);
	writer->Update();
	writer->Write();
}


int _cdecl main(int argc, char* argv[])
{
	//Set Parameters from here itself 
	

	ExecuteRayTracing(argv[1],argv[2]);
	return 1;
}










