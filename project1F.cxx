#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <algorithm>
#include <cmath>
#define NORMALS

using std::cerr;
using std::endl;
int triangleID = -1;
double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}

vtkImageData * NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void WriteImage(vtkImageData *img, const char *filename)
{
	std::string full_filename = filename;
	full_filename += ".png";
	vtkPNGWriter *writer = vtkPNGWriter::New();
	writer->SetInputData(img);
	writer->SetFileName(full_filename.c_str());
	writer->Write();
	writer->Delete();
}

class Triangle
{
	public:
		double	X[3];
		double	Y[3];
		double	Z[3];	
		double  colors[3][3];
		double 	normals[3][3];
		double 	shading[3];
	public:
		void Print()
		{
			printf("TriangleID: %d, X0: %f, Y0: %f, Z0: %f, X1: %f, Y1: %f, Z1: %f, X2: %f, Y2: %f, Z2: %f\n ",
					triangleID,X[0],Y[0],Z[0],X[1],Y[1],Z[1],X[2],Y[2],Z[2]);
		}
	// would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  	public:
		unsigned char *buffer;
		double Zbuffer[1000*1000];
		int width, height;
		void SetPixels(double r, double c, double *color, double z, double shade)
		{
			//cerr << "Triangle " << triangleID << " is writing to row " << r << ", column " << c << endl;
			if (r < 0 || r >= height || c < 0 || c >= width)
				return;
			int index = r*width + c;
			if (z > Zbuffer[index])
			{
				buffer[3*index+0] = ceil_441(255*fmin(1,color[0]*shade));
				buffer[3*index+1] = ceil_441(255*fmin(1,color[1]*shade));
				buffer[3*index+2] = ceil_441(255*fmin(1,color[2]*shade));
				Zbuffer[index] = z;
			}
		}
};

std::vector<Triangle> GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

class Vertices
{
	public:
		double X;
		double Y;
		double Z;
		double color[3];
		double shade;
};

double lerp(double A, double B, double FA, double FB, double X)
{
	double FX;
	if (B == A)
		return FB;	//FX = FB or FA
	if (X == A)
		return FA;	//FX = FA;
	double t = (X - A) / (B - A);
	FX = FA + t*(FB - FA);
	return FX;
}

void RasterizeDownTriangle(Triangle tri, Screen *screen)
{
	Vertices Bottomvertex,Leftvertex,Rightvertex;
    double YMin,YMax;
    YMin = fmin(fmin(tri.Y[0],tri.Y[1]),tri.Y[2]);
    YMax = fmax(fmax(tri.Y[0],tri.Y[1]),tri.Y[2]);

    int count; // identifying bottom
    if (YMin == tri.Y[0]) {
    	Bottomvertex.Y = tri.Y[0];
    	Bottomvertex.X = tri.X[0];
		Bottomvertex.Z = tri.Z[0];
		Bottomvertex.shade = tri.shading[0];
		for (int i = 0; i < 3; i++)
			Bottomvertex.color[i] = tri.colors[0][i];
    	count = 0;
    } else {
    	if (YMin == tri.Y[1]) {
    		Bottomvertex.Y = tri.Y[1];
    		Bottomvertex.X = tri.X[1];
			Bottomvertex.Z = tri.Z[1];
			Bottomvertex.shade = tri.shading[1];
			for (int i = 0; i < 3; i++)
				Bottomvertex.color[i] = tri.colors[1][i];
    		count = 1;
    	} else { 
    		Bottomvertex.Y = tri.Y[2];
    		Bottomvertex.X = tri.X[2];
			Bottomvertex.Z = tri.Z[2];
			Bottomvertex.shade = tri.shading[2];
			for (int i = 0; i < 3; i++)
				Bottomvertex.color[i] = tri.colors[2][i];
    		count = 2;
    	}
    }

    // assigning the remaining vertices
    if (count == 0) {
    	if (tri.X[1] < tri.X[2]) {
    		Leftvertex.X = tri.X[1];
    		Leftvertex.Y = tri.Y[1];
			Leftvertex.Z = tri.Z[1];
			Leftvertex.shade = tri.shading[1];
    		Rightvertex.X = tri.X[2];
    		Rightvertex.Y = tri.Y[2];
			Rightvertex.Z = tri.Z[2];
			Rightvertex.shade = tri.shading[2];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[1][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[2][i];
			}
    	} else {
    		Leftvertex.X = tri.X[2];
    		Leftvertex.Y = tri.Y[2];
			Leftvertex.Z = tri.Z[2];
			Leftvertex.shade = tri.shading[2];
    		Rightvertex.X = tri.X[1];
    		Rightvertex.Y = tri.Y[1];
			Rightvertex.Z = tri.Z[1];
			Rightvertex.shade = tri.shading[1];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[2][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[1][i];
			}
    	}
    }
    if (count == 1) {
    	if (tri.X[0] < tri.X[2]) {
    		Leftvertex.X = tri.X[0];
    		Leftvertex.Y = tri.Y[0];
			Leftvertex.Z = tri.Z[0];
			Leftvertex.shade = tri.shading[0];
    		Rightvertex.X = tri.X[2];
    		Rightvertex.Y = tri.Y[2];
			Rightvertex.Z = tri.Z[2];
			Rightvertex.shade = tri.shading[2];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[0][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[2][i];
			}
    	} else {
    		Leftvertex.X = tri.X[2];
    		Leftvertex.Y = tri.Y[2];
			Leftvertex.Z = tri.Z[2];
			Leftvertex.shade = tri.shading[2];
    		Rightvertex.X = tri.X[0];
    		Rightvertex.Y = tri.Y[0];
			Rightvertex.Z = tri.Z[0];
			Rightvertex.shade = tri.shading[0];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[2][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[0][i];
			}
    	}
    }  
    if (count == 2) {
    	if (tri.X[0] < tri.X[1]) {
    		Leftvertex.X = tri.X[0];
    		Leftvertex.Y = tri.Y[0];
			Leftvertex.Z = tri.Z[0];
			Leftvertex.shade = tri.shading[0];
    		Rightvertex.X = tri.X[1];
    		Rightvertex.Y = tri.Y[1];
			Rightvertex.Z = tri.Z[1];
			Rightvertex.shade = tri.shading[1];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[0][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[1][i];
			}
    	} else {
    		Leftvertex.X = tri.X[1];
    		Leftvertex.Y = tri.Y[1];
			Leftvertex.Z = tri.Z[1];
			Leftvertex.shade = tri.shading[1];
    		Rightvertex.X = tri.X[0];
    		Rightvertex.Y = tri.Y[0];
			Rightvertex.Z = tri.Z[0];
			Rightvertex.shade = tri.shading[0];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[1][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[0][i];
			}
    	}
    }

	// scanlining
    double rowMin = ceil_441(YMin);
    double rowMax = floor_441(YMax);
    for (double r = rowMin; r <= rowMax; r++)
    {
      	//calculate the intersections
      	double leftEnd,rightEnd;
		leftEnd = lerp(Leftvertex.Y,Bottomvertex.Y,Leftvertex.X,Bottomvertex.X,r);
		rightEnd = lerp(Rightvertex.Y,Bottomvertex.Y,Rightvertex.X,Bottomvertex.X,r);
		double leftColor[3],rightColor[3];
		for (int j = 0; j < 3; j++) { // calculate multiple color channels for leftEnd and rightEnd
			leftColor[j] = lerp(Leftvertex.Y,Bottomvertex.Y,Leftvertex.color[j],Bottomvertex.color[j],r); 
			rightColor[j] = lerp(Rightvertex.Y,Bottomvertex.Y,Rightvertex.color[j],Bottomvertex.color[j],r);
		}
		//interpolating Z(leftEnd),Z(rightEnd)
		double ZLeft = lerp(Leftvertex.Y,Bottomvertex.Y,Leftvertex.Z,Bottomvertex.Z,r);
		double ZRight = lerp(Rightvertex.Y,Bottomvertex.Y,Rightvertex.Z,Bottomvertex.Z,r);
		double leftShade = lerp(Leftvertex.Y,Bottomvertex.Y,Leftvertex.shade,Bottomvertex.shade,r);
		double rightShade = lerp(Rightvertex.Y,Bottomvertex.Y,Rightvertex.shade,Bottomvertex.shade,r);
		//printf("Rasterizing along row %f with left end = %f (Z: %f, RGB = %f/%f/%f) and right end = %f (Z: %f, RGB = %f/%f/%f)\n",r,leftEnd,ZLeft,leftColor[0],leftColor[1],leftColor[2],rightEnd,ZRight,rightColor[0],rightColor[1],rightColor[2]);
      	for (double c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++)
      	{
			//interpolating Z for each column of each scanline
			double Zcol = lerp(leftEnd,rightEnd,ZLeft,ZRight,c);
			double columnColor[3];
			double columnShade = lerp(leftEnd,rightEnd,leftShade,rightShade,c);
			for (int j = 0; j < 3; j++)
				columnColor[j] = lerp(leftEnd,rightEnd,leftColor[j],rightColor[j],c); //calculating colors for each column of each scanline			
			//printf("\tGot fragment r = %f, c = %f, z = %f, color = %f/%f/%f\n",r,c,Zcol,columnColor[0],columnColor[1],columnColor[2]);
        	screen->SetPixels(r,c,columnColor,Zcol,columnShade);
      	}
    }
}

void RasterizeUpTriangle(Triangle tri, Screen *screen)
{
    Vertices Topvertex,Leftvertex,Rightvertex;
    double YMin,YMax;
    YMin = fmin(fmin(tri.Y[0],tri.Y[1]),tri.Y[2]);
    YMax = fmax(fmax(tri.Y[0],tri.Y[1]),tri.Y[2]);

    int count; // identifying top
    if (YMax == tri.Y[0]) {
    	Topvertex.Y = tri.Y[0];
    	Topvertex.X = tri.X[0];
		Topvertex.Z = tri.Z[0];
		Topvertex.shade = tri.shading[0];
		for (int i = 0; i < 3; i++)
			Topvertex.color[i] = tri.colors[0][i];
    	count = 0;
    } else {
    	if (YMax == tri.Y[1]) {
    		Topvertex.Y = tri.Y[1];
    		Topvertex.X = tri.X[1];
			Topvertex.Z = tri.Z[1];
			Topvertex.shade = tri.shading[1];
			for (int i = 0; i < 3; i++)
				Topvertex.color[i] = tri.colors[1][i];
    		count = 1;
    	} else { 
    		Topvertex.Y = tri.Y[2];
    		Topvertex.X = tri.X[2];
			Topvertex.Z = tri.Z[2];
			Topvertex.shade = tri.shading[2];
			for (int i = 0; i < 3; i++)
				Topvertex.color[i] = tri.colors[2][i];
    		count = 2;
    	}
    }
    // assigning the remaining vertices
    if (count == 0) {
    	if (tri.X[1] < tri.X[2]) {
    		Leftvertex.X = tri.X[1];
    		Leftvertex.Y = tri.Y[1];
			Leftvertex.Z = tri.Z[1];
			Leftvertex.shade = tri.shading[1];
    		Rightvertex.X = tri.X[2];
    		Rightvertex.Y = tri.Y[2];
			Rightvertex.Z = tri.Z[2];
			Rightvertex.shade = tri.shading[2];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[1][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[2][i];
			}
    	} else {
    		Leftvertex.X = tri.X[2];
    		Leftvertex.Y = tri.Y[2];
			Leftvertex.Z = tri.Z[2];
			Leftvertex.shade = tri.shading[2];
    		Rightvertex.X = tri.X[1];
    		Rightvertex.Y = tri.Y[1];
			Rightvertex.Z = tri.Z[1];
			Rightvertex.shade = tri.shading[1];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[2][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[1][i];
			}
    	}
    }
    if (count == 1) {
    	if (tri.X[0] < tri.X[2]) {
    		Leftvertex.X = tri.X[0];
    		Leftvertex.Y = tri.Y[0];
			Leftvertex.Z = tri.Z[0];
			Leftvertex.shade = tri.shading[0];
    		Rightvertex.X = tri.X[2];
    		Rightvertex.Y = tri.Y[2];
			Rightvertex.Z = tri.Z[2];
			Rightvertex.shade = tri.shading[2];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[0][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[2][i];
			}
    	} else {
    		Leftvertex.X = tri.X[2];
    		Leftvertex.Y = tri.Y[2];
			Leftvertex.Z = tri.Z[2];
			Leftvertex.shade = tri.shading[2];
    		Rightvertex.X = tri.X[0];
    		Rightvertex.Y = tri.Y[0];
			Rightvertex.Z = tri.Z[0];
			Rightvertex.shade = tri.shading[0];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[2][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[0][i];
			}
    	}
    } 
    if (count == 2) {
    	if (tri.X[0] < tri.X[1]) {
    		Leftvertex.X = tri.X[0];
    		Leftvertex.Y = tri.Y[0];
			Leftvertex.Z = tri.Z[0];
			Leftvertex.shade = tri.shading[0];
    		Rightvertex.X = tri.X[1];
    		Rightvertex.Y = tri.Y[1];
			Rightvertex.Z = tri.Z[1];
			Rightvertex.shade = tri.shading[1];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[0][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[1][i];
			}
    	} else {
    		Leftvertex.X = tri.X[1];
    		Leftvertex.Y = tri.Y[1];
			Leftvertex.Z = tri.Z[1];
			Leftvertex.shade = tri.shading[1];
    		Rightvertex.X = tri.X[0];
    		Rightvertex.Y = tri.Y[0];
			Rightvertex.Z = tri.Z[0];
			Rightvertex.shade = tri.shading[0];
			for (int i = 0; i < 3; i++) {
				Leftvertex.color[i] = tri.colors[1][i]; //identifying the color channels for LEFT/RIGHT vertices
				Rightvertex.color[i] = tri.colors[0][i];
			}
    	}
    }
	// scanlining
    double rowMin = ceil_441(YMin);
    double rowMax = floor_441(YMax);
    for (double r = rowMin; r <= rowMax; r++)
    {
	    //calculate the intersections
	    double leftEnd,rightEnd;
		leftEnd = lerp(Leftvertex.Y,Topvertex.Y,Leftvertex.X,Topvertex.X,r);
		rightEnd = lerp(Rightvertex.Y,Topvertex.Y,Rightvertex.X,Topvertex.X,r);
		double leftColor[3],rightColor[3];
		for (int j = 0; j < 3; j++) { // calculate multiple color channels for leftEnd and rightEnd
			leftColor[j] = lerp(Leftvertex.Y,Topvertex.Y,Leftvertex.color[j],Topvertex.color[j],r); 
			rightColor[j] = lerp(Rightvertex.Y,Topvertex.Y,Rightvertex.color[j],Topvertex.color[j],r);
		}
		//interpolating Z(leftEnd),Z(rightEnd)
		double ZLeft = lerp(Leftvertex.Y,Topvertex.Y,Leftvertex.Z,Topvertex.Z,r);
		double ZRight = lerp(Rightvertex.Y,Topvertex.Y,Rightvertex.Z,Topvertex.Z,r);
		//interpolating shading values
		double leftShade = lerp(Leftvertex.Y,Topvertex.Y,Leftvertex.shade,Topvertex.shade,r);
		double rightShade = lerp(Rightvertex.Y,Topvertex.Y,Rightvertex.shade,Topvertex.shade,r);
		//printf("Rasterizing along row %f with left end = %f (Z: %f, RGB = %f/%f/%f) and right end = %f (Z: %f, RGB = %f/%f/%f)\n",r,leftEnd,ZLeft,leftColor[0],leftColor[1],leftColor[2],rightEnd,ZRight,rightColor[0],rightColor[1],rightColor[2]);
     	for (double c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++)
      	{	
			//interpolating Z for each column of each scanline
			double Zcol = lerp(leftEnd,rightEnd,ZLeft,ZRight,c);
			double columnColor[3];
			//interpolating shading for each column
			double columnShade = lerp(leftEnd,rightEnd,leftShade,rightShade,c);
			//calculating colors for each column of each scanline
			for (int j = 0; j < 3; j++)
				columnColor[j] = lerp(leftEnd,rightEnd,leftColor[j],rightColor[j],c);		
			//printf("\tGot fragment r = %f, c = %f, z = %f, color = %f/%f/%f\n",r,c,Zcol,columnColor[0],columnColor[1],columnColor[2]);
        	screen->SetPixels(r,c,columnColor,Zcol,columnShade);
      	}
    }
}

void RasterizeArbitraryTriangle(Triangle tri, Screen *screen)
{
    Vertices topVertex,bottomVertex,middleVertex;
    double YMin,YMax,middleY;
    Triangle UpTri,DownTri;

    // figuring out Y-cordinates
    YMin = fmin(fmin(tri.Y[0],tri.Y[1]),tri.Y[2]);
    YMax = fmax(fmax(tri.Y[0],tri.Y[1]),tri.Y[2]);
	int indexYMin,indexYMax;
    for (int y = 0; y < 3; y++)	{  // finding index for YMin and YMax
		if (tri.Y[y] == YMin)
			indexYMin = y;
		if (tri.Y[y] == YMax)
			indexYMax = y;
    }
	for (int mid = 0; mid < 3; mid++) {
		if ((mid != indexYMin) && (mid != indexYMax))
			middleY = tri.Y[mid]; //assigning middleY
	}
	// Determine if this triangle is arbitrary or not
	if (middleY == YMax)
		RasterizeDownTriangle(tri,screen);
	else if (middleY == YMin)
		RasterizeUpTriangle(tri,screen);
	else { //if Arbitrary
		for (int j = 0; j < 3; j++)
		{
			if (tri.Y[j] == YMax)
			{		
				topVertex.Y = tri.Y[j];
				topVertex.X = tri.X[j];
				topVertex.Z = tri.Z[j];
				topVertex.shade = tri.shading[j];
				for (int c = 0; c < 3; c++)
					topVertex.color[c] = tri.colors[j][c]; //identifying which color channels(RGB) belong to which vertex
			}
			if (tri.Y[j] == YMin)
			{
				bottomVertex.Y = tri.Y[j];
				bottomVertex.X = tri.X[j];
				bottomVertex.Z = tri.Z[j];
				bottomVertex.shade = tri.shading[j];
				for (int c = 0; c < 3; c++)
					bottomVertex.color[c] = tri.colors[j][c]; //identifying which color channels belong to which vertex
			}
			if (tri.Y[j] == middleY)
			{
				middleVertex.Y = tri.Y[j];
				middleVertex.X = tri.X[j];
				middleVertex.Z = tri.Z[j];
				middleVertex.shade = tri.shading[j];
				for (int c = 0; c < 3; c++)
					middleVertex.color[c] = tri.colors[j][c]; //identifying which color channels belong to which vertex
			}
		}

		// calculate all the essentials for the split point
		double XSlope = lerp(bottomVertex.Y,topVertex.Y,bottomVertex.X,topVertex.X,middleY);
		double ZSlope = lerp(bottomVertex.Y,topVertex.Y,bottomVertex.Z,topVertex.Z,middleY); // Z interpolation
		double ShadeSlope = lerp(bottomVertex.Y,topVertex.Y,bottomVertex.shade,topVertex.shade,middleY);
		//splitting Up Triangle, check for vertical side
		UpTri.X[0] = XSlope;
		UpTri.X[1] = middleVertex.X;
		UpTri.X[2] = topVertex.X;
		UpTri.Y[0] = middleVertex.Y; // split point y
		UpTri.Y[1] = middleVertex.Y;
		UpTri.Y[2] = topVertex.Y;
		UpTri.Z[0] = ZSlope; //split point z
		UpTri.Z[1] = middleVertex.Z;
		UpTri.Z[2] = topVertex.Z;
		UpTri.shading[0] = ShadeSlope; //shade value at split point
		UpTri.shading[1] = middleVertex.shade;
		UpTri.shading[2] = topVertex.shade;
		for (int c = 0; c < 3; c++) 
		{	//color interpolation for split point, and remaining points
			UpTri.colors[0][c] = lerp(bottomVertex.Y,topVertex.Y,bottomVertex.color[c],topVertex.color[c],middleY);
			UpTri.colors[1][c] = middleVertex.color[c];
			UpTri.colors[2][c] = topVertex.color[c];
		}
		RasterizeUpTriangle(UpTri,screen);
		//splitting Down Triangle
		DownTri.X[0] = XSlope; // split point x
		DownTri.X[1] = middleVertex.X;
		DownTri.X[2] = bottomVertex.X;
		DownTri.Y[0] = middleVertex.Y; // split point y
		DownTri.Y[1] = middleVertex.Y;
		DownTri.Y[2] = bottomVertex.Y;
		DownTri.Z[0] = ZSlope; // split point z
		DownTri.Z[1] = middleVertex.Z;
		DownTri.Z[2] = bottomVertex.Z;
		DownTri.shading[0] = ShadeSlope; //shade value at split point
		DownTri.shading[1] = middleVertex.shade;
		DownTri.shading[2] = bottomVertex.shade;
		for (int c = 0; c < 3; c++)
		{	
			DownTri.colors[0][c] = lerp(bottomVertex.Y,topVertex.Y,bottomVertex.color[c],topVertex.color[c],middleY);
			DownTri.colors[1][c] = middleVertex.color[c];
			DownTri.colors[2][c] = bottomVertex.color[c];
		}
		RasterizeDownTriangle(DownTri,screen);
	}
}

double Dot(double *A, double *B)
{
	double result = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
	return result;
}

void Cross(double *A, double *B, double *result)
{
	result[0] = A[1]*B[2] - A[2]*B[1]; //result.X = A.Y*B.Z - A.Z*B.Y;
	result[1] = B[0]*A[2] - A[0]*B[2]; //result.Y = B.X*A.Z - A.X*B.Z;
	result[2] = A[0]*B[1] - A[1]*B[0]; //result.Z = A.X*B.Y - A.Y*B.X;
}

void Normalize(double *A)
{
	double norm = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
	for (int i = 0; i < 3; i++)
		A[i] = A[i] / norm;
}

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0] + ptIn[1]*A[1][0] + ptIn[2]*A[2][0] + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1] + ptIn[1]*A[1][1] + ptIn[2]*A[2][1] + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2] + ptIn[1]*A[1][2] + ptIn[2]*A[2][2] + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3] + ptIn[1]*A[1][3] + ptIn[2]*A[2][3] + ptIn[3]*A[3][3];
}

class Camera
{
 	public:
		double          near, far;
		double          angle;
		double          position[3];
		double          focus[3];
		double          up[3];

		Matrix          ViewTransform(void);
		Matrix          CameraTransform(void);
		Matrix          DeviceTransform(void);
};

Matrix Camera::DeviceTransform(void)
{
	int n = 1000;
	int m = 1000;
	Matrix result;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			result.A[i][j] = 0; //initializing the matrix, all 0's
	result.A[0][0] = n/2;
	result.A[1][1] = m/2;
	result.A[2][2] = 1;
	result.A[3][3] = 1;
	result.A[3][0] = n/2;
	result.A[3][1] = m/2;
	return result;	
}

Matrix Camera::ViewTransform(void)
{
	double alpha = angle;
	double n = near;
	double f = far;
	Matrix result;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			result.A[i][j] = 0; //initializing the matrix, all 0's
	result.A[0][0] = 1/tan(alpha/2);
	result.A[1][1] = 1/tan(alpha/2);
	result.A[2][2] = (f+n)/(f-n);
	result.A[2][3] = -1;
	result.A[3][2] = 2*f*n/(f-n);
	return result;
}

Matrix Camera::CameraTransform(void)
{
	//Calculate Camera Frame (O,v,u,w)
	double O[3],w[3],u[3],v[3],O_focus[3],Up[3],t[3];
	for (int i = 0; i < 3; i++)
	{
		Up[i] = up[i];
		O[i] = position[i];
		O_focus[i] = O[i] - focus[i];
		w[i] = O_focus[i];
		t[i] = 0 - O[i];
	}
	Cross(Up,O_focus,u);
	Cross(O_focus,u,v);
	//normalize v,u,w
	Normalize(v);
	Normalize(u);
	Normalize(w);
	//printf("Camera Frame: U = %f, %f, %f\n",u[0],u[1],u[2]);
	//printf("Camera Frame: V = %f, %f, %f\n",v[0],v[1],v[2]);
	//printf("Camera Frame: W = %f, %f, %f\n",w[0],w[1],w[2]);
	//printf("Camera Frame: O = %f, %f, %f\n",O[0],O[1],O[2]);
	//calculate Camera Transform Matrix
	Matrix result;
	for (int j = 0; j < 4; j++) {
		result.A[j][0] = u[j]; //or *(u+j)
		result.A[j][1] = v[j]; //or *(v+j)
		result.A[j][2] = w[j]; //or *(w+j)
		result.A[j][3] = 0;
		if (j == 3) {
			result.A[j][0] = Dot(u,t);
			result.A[j][1] = Dot(v,t);
			result.A[j][2] = Dot(w,t);
			result.A[j][3] = 1;
		}
	}
	return result;
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

Triangle TransformTrianglesToDeviceSpace(Triangle T,Matrix M)
{
	// apply M to each vertex of T
	double vIn[3][4]; // to store vertices
	for (int v = 0; v < 3; v++) 
	{	//pull out vertices to be transformed
		vIn[v][0] = T.X[v];
		vIn[v][1] = T.Y[v];
		vIn[v][2] = T.Z[v];
		vIn[v][3] = 1;
	}
	double vOut[3][4];
	M.TransformPoint(vIn[0],vOut[0]);
	M.TransformPoint(vIn[1],vOut[1]);
	M.TransformPoint(vIn[2],vOut[2]);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 4; j++)
			vOut[i][j] = vOut[i][j]/vOut[i][3]; // divide back by w
	//copying newly-calculated vertices to new triangle 
	for (int t = 0; t < 3; t++) 
	{	
		T.X[t] = vOut[t][0];
		T.Y[t] = vOut[t][1];
		T.Z[t] = vOut[t][2];
	}
	//printf("Transformed V0 from (%f, %f, %f) to (%f, %f, %f)\n",vIn[0][0],vIn[0][1],vIn[0][2],T.X[0],T.Y[0],T.Z[0]);
	//printf("Transformed V1 from (%f, %f, %f) to (%f, %f, %f)\n",vIn[1][0],vIn[1][1],vIn[1][2],T.X[1],T.Y[1],T.Z[1]);
	//printf("Transformed V2 from (%f, %f, %f) to (%f, %f, %f)\n",vIn[2][0],vIn[2][1],vIn[2][2],T.X[2],T.Y[2],T.Z[2]);
	return T;
}

void SaveImage(vtkImageData *img, int frame)
{
	char filename[32];
	sprintf(filename, "frame%03d", frame);
	WriteImage(img, filename);
}

void InitializeScreen(Screen *screen, vtkImageData *image)
{
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
	int npixels = 1000*1000;
	for (int i = 0 ; i < npixels*3 ; i++)
	{
		buffer[i] = 0;
	}
	screen->buffer = buffer;
	for (int j = 0 ; j < npixels; j++) {
		screen->Zbuffer[j] = -1;
	}
}

struct LightingParameters
{
    LightingParameters(void)
    {
		lightDir[0] = -0.6;
		lightDir[1] = 0;
		lightDir[2] = -0.8;
		Ka = 0.3;
		Kd = 0.7;
		Ks = 2.3;
		alpha = 2.5;
    };

    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};
LightingParameters lp;

double CalculatePhongShading(LightingParameters &lp, double *viewDirection, double *normal)
{
	double R[3];
	Normalize(lp.lightDir);
	Normalize(normal);
	for (int i = 0; i < 3; i++)
		R[i] = 2*Dot(lp.lightDir,normal)*normal[i] - lp.lightDir[i];
	Normalize(viewDirection);
	Normalize(R);
	double VdotR = Dot(viewDirection,R);
	double shading_amount = lp.Ka + lp.Kd*std::abs(Dot(lp.lightDir,normal)) + lp.Ks*fmax(0,pow(VdotR,lp.alpha));
	return shading_amount;
}

int main()
{
	vtkImageData *image = NewImage(1000, 1000);

	std::vector<Triangle> triangles = GetTriangles();
	Screen screen; 	//allocate the screen
	screen.width = 1000;
	screen.height = 1000;

	for (int j = 0; j < 1; j++) {
		InitializeScreen(&screen,image);
		Camera c = GetCamera(j,1000);
		Matrix cameraMatrix = c.CameraTransform(); //calculate camera-transform matrix
		Matrix viewMatrix = c.ViewTransform(); 	//Calculate View-Transform matrix
		Matrix deviceMatrix = c.DeviceTransform(); 		//Calculate Device-Transform matrix
		Matrix tempM = cameraMatrix.ComposeMatrices(cameraMatrix,viewMatrix);  		//Compose those 3 matrix into M
		Matrix M = tempM.ComposeMatrices(tempM,deviceMatrix);
		for (int i = 0; i < triangles.size(); i++) {
			triangleID = i;
			//printf("\nWorking on triangle %d\n",triangleID);
			Triangle T = triangles[i];
			//Calculate shading value AND view direction for each vertex
			for (int t = 0; t < 3; t++)
			{	
				double viewDir[3];
				viewDir[0] = c.position[0] - T.X[t];
				viewDir[1] = c.position[1] - T.Y[t];
				viewDir[2] = c.position[2] - T.Z[t];
				T.shading[t] = CalculatePhongShading(lp,viewDir,T.normals[t]);
			}
			Triangle DeviceTri = TransformTrianglesToDeviceSpace(T,M);
			RasterizeArbitraryTriangle(DeviceTri,&screen);
		}
		SaveImage(image,j);
	}
}