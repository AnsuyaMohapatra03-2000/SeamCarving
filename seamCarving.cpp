#include "opencv2/opencv.hpp"
#include <stdio.h>
#include <stdint.h>
#include <string.h>
using namespace cv;
using namespace std;

class SeamCarver {
		 Mat image;
	  public:
   SeamCarver(string name){
	       image = imread(name,-1);
   }               
   void showpic(){
		if (!image.data) {
			printf("No image data \n");
		}
       imshow("Gradient", image);
   }
     int width()  {
			return image.rows;
	 }
	 int height()  {
			return image.cols;
	 }
	 int calxdirection(int x,int y){
								Vec3b left,right;
								if(y==0){
									left = image.at<Vec3b>(x,image.cols-1);
								}else{
									left = image.at<Vec3b>(x,y-1);
								}
								if(y==image.cols-1){
									right = image.at<Vec3b>(x,0);
								}else{
									right = image.at<Vec3b>(x,y+1);
								}
								int b =  abs(static_cast<int>(left[0])-static_cast<int>(right[0])); 
								int g =  abs(static_cast<int>(left[1])-static_cast<int>(right[1])); 
								int r =  abs(static_cast<int>(left[2])-static_cast<int>(right[2])); 

								return pow(b,2)+pow(g,2)+pow(r,2);
	 }
	  int calydirection(int x,int y){
							Vec3b top,down;
							if(x==0){
								top = image.at<Vec3b>(image.rows-1,y);
							}else{
								top = image.at<Vec3b>(x-1,y);
							}
							if(x==image.rows-1){
								down = image.at<Vec3b>(image.rows-1,y);
							}else{
								down = image.at<Vec3b>(x+1,y);
							}
							int b =  abs(static_cast<int>(top[0])-static_cast<int>(down[0])); 
							int g =  abs(static_cast<int>(top[1])-static_cast<int>(down[1])); 
							int r =  abs(static_cast<int>(top[2])-static_cast<int>(down[2])); 

							return pow(b,2)+pow(g,2)+pow(r,2);
	 }
    double energy(int x, int y){
		   int xdir = calxdirection(x,y);
		   int ydir = calydirection(x,y);
		    return sqrt(xdir+ydir);
	}
int* findHorizontalSeam(double** cost, int r, int c) {
 int *path = new int[r];

    int minindex = 0;
    for (int i = 1; i < r; ++i) {
        if (cost[i][c - 1] < cost[minindex][c - 1]) {
            minindex = i;
        }
    }

    int x = minindex;
    int y = c - 1;
    path[y] = x;

    while (y != 0) {  
        if (x == 0) {
            if (cost[x][y - 1] < cost[x + 1][y - 1]) {
                y = y - 1;
            } else {
                x = x + 1;
                y = y - 1;
            }
        } else if (x == r - 1) {
            if (cost[x - 1][y - 1] < cost[x][y - 1]) {
                x = x - 1;
                y = y - 1;
            } else {
                y = y - 1;
            }
        } else {
            double a = cost[x - 1][y - 1];
            double b = cost[x][y - 1];
            double c = cost[x + 1][y - 1];

            if (a <= b && a <= c) {
                x = x - 1;
            } else if (c <= a && c <= b) {
                x = x + 1;
            }
            y = y - 1;
        }
        path[y] = x;
    }

    return path;
   }
    
   int* findVerticalSeam(double** cost,int r,int c) {
	int *path = new int[r];

    int minindex = 0;
    for (int i = 0; i < c; ++i) {
        if (cost[r - 1][i] < cost[r - 1][minindex]) {
            minindex = i;
        }
    }

    int x = r - 1; 
    int y = minindex;
    path[x] = y;

    while (x != 0) { 
        if (y == 0) {
            if (cost[x - 1][y] < cost[x - 1][y + 1]) {
                x = x - 1;
            } else {
                x = x - 1;
                y = y + 1;
            }
        } else if (y == c - 1) {
            if (cost[x - 1][y - 1] < cost[x - 1][y]) {
                x = x - 1;
                y = y - 1;
            } else {
                x = x - 1;
            }
        } else {
            double a = cost[x - 1][y - 1];
            double b = cost[x - 1][y];
            double c = cost[x - 1][y + 1];

            if (a <= b && a <= c) {
                x = x - 1;
                y = y - 1;
            } else if (b <= a && b <= c) {
                x = x - 1;
            } else if (c <= a && c <= b) {
                x = x - 1;
                y = y + 1;
            }
        }

        path[x] = y;
    }

    return path;
   }     

   double** findcostMatrix(double** energyval,int r,int c) {
			double** g = new double*[r];  

			for (int i = 0; i < r; ++i) {
				g[i] = new double[c];
				for (int j = 0; j < c; ++j) {
					g[i][j] = 0; 
				}
			}
			for (int j = 0; j < c; ++j) {
					g[0][j] = energyval[0][j];  
			}
		  for (int i = 1; i < r; ++i) {
				for (int j = 0; j < c; ++j) {
					if(j==0)
					     g[i][j] = energyval[i][j] + min(g[i-1][j],g[i-1][j+1]);
					else if(j==c-1)
					    g[i][j] = energyval[i][j] + min(g[i-1][j-1],g[i-1][j]); 
					else
					     g[i][j] = energyval[i][j] + min(min(g[i-1][j],g[i-1][j+1]),g[i-1][j-1]);
				}
			}
        return g;
   }           
double ** removepixelfromcost(int *ptr,double ** g,int row,int col){
	     double **ans;
		 double matrix[row][col-1];
		 int i=0;
		 int flag=0;
		 while(i<row){
			   int y=ptr[i];
				for(int r=0;r<row;r++){
					   flag=0;
						for(int c=0,q=0;c<col;c++){
							if(r==i && c==y){  flag=1; }
							if(c!=col-1 && flag==1){
								matrix[r][c] = g[r][c+1];
							}else{
								matrix[r][c] = g[r][c];
							}
						}
				}
				i=i+2;
		 }
		for (int i = 0; i < row; ++i) {
			ans[i] = new double[col];  
			for (int j = 0; j < col; ++j) {
				ans[i][j] = matrix[i][j];
			}
        }
		return ans;
}
void removeVerticalSeam(int* seam,int size) {
	     
			int rows = image.rows;
			int cols = image.cols;


			cv::Mat mod(rows, cols - 1, CV_8UC3);

			for (int row = 0; row < rows; ++row) {
				int seamIndex = seam[row]; 
				int modCol = 0; 

				for (int col = 0; col < cols; ++col) {
					if (col != seamIndex) {
						mod.at<cv::Vec3b>(row, modCol++) = image.at<cv::Vec3b>(row, col);
					}
				} 
              }

			image = mod;
   } 

void removeHorizontalSeam(int *seam, int size) {
    int rows = image.rows;
    int cols = image.cols;

    cv::Mat mod(rows - 1, cols, image.type());

    for (int c = 0; c < cols; ++c) {
        int seamIndex = seam[c]; 
        int modRow = 0; 

        for (int r = 0; r < rows; ++r) {
            if (r != seamIndex) {
                mod.at<cv::Vec3b>(modRow++, c) = image.at<cv::Vec3b>(r, c);
            }
        }
    }
    image = mod;
}

};
int main(int argc, char** argv)
{
	if(argc!=2){
		cout<<"Enter the path of the image file !\n";
		return 0;
	}
	string name=string(argv[1]);
	SeamCarver s(name);
	//cout<<"height :: "<<s.width()<<" width :: "<<s.height()<<"\n";

	int l,b;
	//cout<<"Enter the height of the image that you want ? ";
	cin>>l;
	//cout<<"Enter the width of the image that you want ? ";
	cin>>b;

 while(s.height()>b){
	             
	                     	int row = s.width();
							int col = s.height();
	
						double energyval[row][col];
						for(int r=0;r<row;r++)
						{
								for(int c=0;c<col;c++){
									energyval[r][c] =  s.energy(r,c);
								}
						}
						double** g = new double*[row]; 
                   
						for (int i = 0; i < row; ++i) {
								g[i] = new double[col];  
								for (int j = 0; j < col; ++j) {
									g[i][j] = energyval[i][j]; 
								}
						}
					
						double **p = s.findcostMatrix(g,row,col);
					
						int *ptr = s.findVerticalSeam(p,row,col);
				
						s.removeVerticalSeam(ptr,row); 
						
						for (int i = 0; i < row; ++i) {
							delete[] g[i];
						}
						for (int i = 0; i < row; ++i) {
							delete[] p[i];
						}
						delete[] p;
						delete[] g;
						delete[] ptr;
						
}
 while(s.width()>l){
	               
	                     	int row = s.width();
							int col = s.height();
					
						double energyval[row][col];
						for(int r=0;r<row;r++)
						{
								for(int c=0;c<col;c++){
									energyval[r][c] =  s.energy(r,c);
								}
						}
						double** g = new double*[row]; 
                        
						for (int i = 0; i < row; ++i) {
								g[i] = new double[col];  
								for (int j = 0; j < col; ++j) {
									g[i][j] = energyval[i][j];  
								}
						}
						double **p = s.findcostMatrix(g,row,col);
						int *ptr = s.findHorizontalSeam(p,row,col);
						s.removeHorizontalSeam(ptr,row); 
						for (int i = 0; i < row; ++i) {
							delete[] g[i];
							delete[] p[i];
						}
						delete[] p;
						delete[] g;
						delete[] ptr;
}
	s.showpic();

	   Mat original = imread(name,-1);
		if (! original.data) {
			printf("No image data \n");
		}
       imshow("original",  original);
	waitKey();
	return 0;
}

