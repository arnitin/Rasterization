/*  
    smooth.c
    Nate Robins, 1998

    Model viewer program.  Exercises the glm library.
*/

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <GL/glut.h>
#include "gltb.h"
#include "glm.h"
#include "dirent32.h"

#pragma comment( linker, "/entry:\"mainCRTStartup\"" )  // set the entry point to be main()



#define ILLUMINATION 0.689

#define DATA_DIR "data/"

int        rasterFlag = 0;
int        flatonly = 0; 
int        pointInfo[512*512];
char*      model_file = NULL;		/* name of the obect file */
GLint      entries = 0;			/* entries in model menu */
float      pixels[512*512*3];
double     linesMat[512][512];
GLuint     model_list = 0;		/* display list for object */
GLuint     material_mode = 0;		/* 0=none, 1=color, 2=material */
GLfloat    scale;		        /* original scale factor */
GLfloat    smoothing_angle = 90.0;	/* smoothing angle */
GLfloat    weld_distance = 0.00001;	/* epsilon for welding vertices */
GLdouble   pan_x = 0.0;
GLdouble   pan_y = 0.0;
GLdouble   pan_z = 0.0;
GLMmodel*  model;		        /* glm model data structure */
GLboolean  facet_normal = GL_FALSE;	/* draw with facet normal? */
GLboolean  bounding_box = GL_FALSE;	/* bounding box on? */
GLboolean  performance = GL_FALSE;	/* performance counter on? */
GLboolean  stats = GL_FALSE;		/* statistics on? */



#define T(x) (model->triangles[(x)])

#include <time.h>
#if defined(_WIN32)
#include <sys/timeb.h>
#define CLK_TCK 1000
#else
#include <limits.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#endif


double DIFF_COEFF = 0;
double SPEC_COEFF = 0;
double  AMB_COEFF = 0;

typedef struct PointStr{
  double x;
  double y;
  double z;
} Point;

Point  lightVec,dirVec;

void Minus(Point v1,Point v2, Point * result) {
  result->x = fabs(v1.x - v2.x);
  result->y = fabs(v1.y - v2.y);
  result->z = fabs(v1.z - v2.z);
}

void crossProduct(Point v1,Point v2, Point *result){
  
  result-> x = v1.y*v2.z - v2.y*v1.z;
  result-> y = (v1.x*v2.z - v2.x*v1.z) * (-1);
  result-> z = v1.x*v2.y - v2.x*v1.y;
}


float dotProduct(Point v1, Point v2) {
  return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

void calcBisector(Point * In,Point * a, Point * b) {
  float norm = sqrt((a->x + b->x)* (a->x+b->x) +  (a->y+b->y)*(a->y+b->y) + (a->z+b->z)*(a->z+b->z));
  In->x =  (a->x+b->x)/norm;
  In->y =  (a->y+b->y)/norm;
  In->z =  (a->z+b->z)/norm;
}

void getNormal(Point v1, Point v2,Point v3, Point *result){
  Point temp1,temp2,cp;
  float denom = 0;
  int flag = 0;
  if(result->x == -999){
    flag = 1;
  } 
  Minus(v1,v2,&temp1);
  Minus(v3,v2,&temp2); 
  cp.x = result-> x = temp1.y*temp2.z - temp2.y*temp1.z;
  cp.y = result-> y = (temp1.x*temp2.z - temp2.x*temp1.z) * (-1);
  cp.z = result-> z = temp1.x*temp2.y - temp2.x*temp1.y;
  
  denom = sqrt(result-> x * result->x + result->y * result->y + result->z * result->z);    
  result->x /= denom; 
  result->y /= denom;
  result->z /= denom; 
}

double area(Point one,Point two,Point three){
  return(fabs((one.x*two.y + two.x*three.y + three.x*one.y - one.x*three.y - two.x*one.y - three.x*two.y)/2));
}

//D = x1y2 + x2y3 + x3y1 - x1y3 - x2y1 - x3y2
 
void draw_line(double x0,double y0,double x1,double y1){
 
  x1 = (int)x1;
  x0 = (int)x0;
  y1 = (int)y1;
  y0 = (int)y0;

  int dx = abs(x1-x0), sx = x0<x1 ? 1 : -1;
  int dy = abs(y1-y0), sy = y0<y1 ? 1 : -1; 
  int err = (dx>dy ? dx : -dy)/2, e2;
  int y = 1;
  for(;;){
    y = 0;
    
    if(x0>511 || y0 > 511 || x0<0 || y0<0 ){
      if(x0>511){
	if(y0>511)
	  linesMat[511][511] = 1;
	else
	  linesMat[511][(int)y0] = 1;
      } else if (x0 < 0){
	if(y0>511)
	  linesMat[0][511] = 1;
	else
	  linesMat[0][(int)y0] = 1;
      } else if( y0 > 511){
      	if(x0>512)
	  linesMat[511][511] = 1;
	else
	  linesMat[(int)x0][511] = 1;
      } else {
      	if(x0>511)
	  linesMat[511][0] = 1;
	else
	  linesMat[(int)x0][0] = 1;
      }
      y =1;
    }
    
    if(y==0){
      linesMat[(int)x0][(int)y0] = 1;
    }
    
    if (x0==x1 && y0==y1) {
      break;
    }
    
    e2 = err;
    if (e2 >-dx) { 
      err -= dy; 
      x0 += sx; 
    }
    if (e2 < dy) {
      err += dx;
      y0 += sy; 
    }
  }
}

int checkifIntersects(Point A,Point B,Point  C,double i,double j){
  
  if(linesMat[(int)i][(int)j] == 1){
    linesMat[(int)i][(int)j] = 0;
    return 1;    
  } else {
    return 0;
  }
 
  double maxX,minX,maxY,minY;
  if ((i == abs(A.x) && j == abs(A.y)) || 
      (i == abs(B.x) && j == abs(B.y)) ||
      (i == abs(C.x) && j == abs(C.y))){
    return 0;
  }

}
  
void triangleCoord(int v1,int v2,int v3,double * A,double *B, double * C) {
  int vIndices[3];
  vIndices[0] = v1;
  vIndices[1] = v2;
  vIndices[2] = v3;

  GLdouble modelProj[4*4];
  GLdouble proj[4*4];
  GLint    view[4];
  double   tx[3], ty[3], tz[3];
  glGetDoublev(GL_MODELVIEW_MATRIX, modelProj);
  glGetDoublev(GL_PROJECTION_MATRIX, proj);
  glGetIntegerv(GL_VIEWPORT, view);
 
  gluProject(model->vertices[3*vIndices[0] + 0],
	     model->vertices[3*vIndices[0] + 1],
	     model->vertices[3*vIndices[0] + 2],
	     modelProj,
	     proj,
	     view,
	     &A[0],
	     &A[1],
	     &A[2]);
  gluProject(model->vertices[3*vIndices[1] + 0],
	     model->vertices[3*vIndices[1] + 1],
	     model->vertices[3*vIndices[1] + 2],
	     modelProj,
	     proj,
	     view,
	     &B[0],
	     &B[1],
	     &B[2]);
  gluProject(model->vertices[3*vIndices[2] + 0],
	     model->vertices[3*vIndices[2] + 1],
	     model->vertices[3*vIndices[2] + 2],
	     modelProj,
	     proj,
	     view,
	     &C[0],
	     &C[1],
	     &C[2]);  
}

void fillCoords(Point * In,double * Arr){
  In->x = Arr[0];
  In->y = Arr[1];
  In->z = Arr[2];
}

void fillCoordsStr(Point * In,Point  From){
  In->x = (From.x);
  In->y = (From.y);
  In->z = From.z;
}

float calculateLighting(Point bisector,Point normalIn) {

  float Ls = 0.0,Ld = 0,La = 0.0; 
  float dotProd = normalIn.x * bisector.x + normalIn.y * bisector.y + normalIn.z * bisector.z;
  
  if (dotProd > 0) {
    dotProd = pow(dotProd,100);
    Ld = SPEC_COEFF *  ILLUMINATION * dotProd;;
    if (Ld > 0) {
    }
  }
  
  dotProd = normalIn.x * (lightVec.x) + normalIn.y * lightVec.y + normalIn.z * lightVec.z;
  
  if (dotProd > 0) {    
    Ls = DIFF_COEFF *  ILLUMINATION * dotProd;
  }
  
  La = AMB_COEFF * ILLUMINATION;
  return (Ld+Ls+La);
}

void rasterization() {
 
  lightVec.x  =  -1;
  lightVec.y  =  0.5;
  lightVec.z  = 0;
  dirVec.x    = 0;
  dirVec.y    = 0;
  dirVec.z    = 1; 

  int    fill = 0;
  int    res = 0,prev = 0;;
  int    i,j,v1Index,v2Index,v3Index;
  int    k = 0,counter = 0;
  int    minx,miny,maxx,maxy; 
  int    left       = 0;
  int    right      = 511;
  int    leftpixel  = -1;
  int    rightpixel = -1;
  float  alpha,beta,gamma;
  Point  A,B,C,p,norA,norB,norC;
  double m   = 0, n = 0, areatot;
  double red = 1, green= 1, blue=1;
  double lightingA,lightingB,lightingC;
  static GLMgroup * group;
  static GLMtriangle * triangle;
  double a[3],b[3],c[3];
  double zBuffer[512][512];
  double c1,c2,c3;

  group = model->groups;
  for(n=0;n<512;n++) {
    for(m=0;m<512;m++) {
       counter = 512 * n + m;
	 k = 3 * counter;
	 pixels[k]   = 1;//red;
	 pixels[k+1] = 1;//green; 
	 pixels[k+2] = 1;//blue;
      //linesMat[(int)m][(int)n]=0;
      zBuffer[(int)m][(int)n]=10;
    }
  }

  memset(linesMat, 0, sizeof(linesMat[0][0]) * 512 * 512);
  // memset(zBuffer, 3E8, sizeof(zBuffer[0][0]) * 512 * 512);
  // memset(pixels, 0, sizeof(float) * 512 *512 *3);
  Point normal,bisector;

  for (i = 0; i <model->numtriangles; i++) {
    v1Index = (float)T(i).vindices[0];
    v2Index = (float)T(i).vindices[1];
    v3Index = (float)T(i).vindices[2];
    triangleCoord((int)v1Index,(int)v2Index,(int)v3Index,a,b,c);
    GLfloat* nor = &model->facetnorms[3 * (i+1)];
    normal.x = nor[0];
    normal.y = nor[1];
    normal.z = nor[2];


    norA.x = model->normals[3*v1Index + 0];
    norA.y =model->normals[3*v1Index + 1];
    norA.z =model->normals[3*v1Index + 2];
    
    norB.x = model->normals[3*v2Index + 0];
    norB.y =model->normals[3*v2Index + 1];
    norB.z =model->normals[3*v2Index + 2];
    
    norC.x = model->normals[3*v3Index + 0];
    norC.y =model->normals[3*v3Index + 1];
    norC.z =model->normals[3*v3Index + 2];
    
    fillCoords(&A,a);
    fillCoords(&B,b);
    fillCoords(&C,c);

    lightingA = calculateLighting(bisector,norA);
    lightingB = calculateLighting(bisector,norB);
    lightingC = calculateLighting(bisector,norC);

    double lighting = 0;
  
    calcBisector(&bisector ,&dirVec,&lightVec);
    if (flatonly) {
      lighting = calculateLighting(bisector,normal);
    }
    memset(linesMat, 0, sizeof(linesMat[0][0]) * 512 * 512);
    draw_line(A.x,A.y,C.x,C.y);
    draw_line(A.x,A.y,B.x,B.y);
    draw_line(B.x,B.y,C.x,C.y);

    minx = A.x<=B.x?A.x:B.x;
    minx = minx<=C.x?minx:C.x;

    maxx = A.x>=B.x?A.x:B.x;
    maxx = maxx>=C.x?maxx:C.x;

    maxy = A.y>=B.y?A.y:B.y;
    maxy = maxy>=C.y?maxy:C.y;

    miny = A.y<=B.y?A.y:B.y;
    miny = miny<=C.y?miny:C.y;

    k       = 0;
    counter = 0;

    if(minx<=0){
      minx = 0;
    }
    
    if(maxx>=511){
      maxx = 511;
    }

    if(miny<=0){
      miny = 0;
    }
    
    if(maxy>=511){
      maxy = 511;
    }


    for(n=miny;n<maxy;n++) {
      fill           = 0;
      res            = 0;
      left           = minx;;
      right          = maxx;
      leftpixel      = -1;
      rightpixel     = -1;
      
      while(left<=right && (leftpixel==-1 || rightpixel==-1)){
	if(leftpixel == -1){
	  if(linesMat[left][(int)n] == 1){
	    leftpixel = left;
	    linesMat[left][(int)n] = 0;	  
	  } else {
	    left++;
	  }
	   
	}
	if(rightpixel == -1){
	  if(linesMat[right][(int)n] == 1){
	    rightpixel = right;
	    linesMat[right][(int)n] = 0;
	  } else{
	    right--;
	  }
	}
      }

      if(leftpixel != -1 && rightpixel == -1){
	int fdsfs=0;
      }
      double Za,Zb,Zp;
      if(left<=right){
	int fillvals =0;
       	for(fillvals=leftpixel;fillvals<=rightpixel;fillvals++){	

	  p.x = fillvals;
	  p.y = n;
	  p.z = 0;
	  
	  c1 = area(p,A,B);
	  c2 = area(p,B,C);
	  c3 = area(p,C,A);
	  areatot = area(A,B,C);
	  
	  double alpha = c2/areatot;
	  double beta  = c3/areatot;
	  double gamma = c1/areatot;

	  Zp = A.z*alpha + B.z*beta + C.z*gamma;

	  counter = 512 * n + fillvals;
	  k = 3 * counter;	

	  if(Zp < zBuffer[fillvals][(int)n] && (Zp>0) && Zp!=0){
	    if (!flatonly) {
	      lighting = lightingA*alpha + lightingB*beta + lightingC*gamma;
	    }
	    zBuffer[fillvals][(int)n] = Zp;
	    pixels[k]   = 0.3 + lighting;
	    pixels[k+1] = 0.0;
	    pixels[k+2] = 0.0;	
	  } 
	}      
      }  
    }
  }
}	


float elapsed(void)
{
  static long begin = 0;
  static long finish, difference;
    
#if defined(_WIN32)
  static struct timeb tb;
  ftime(&tb);
  finish = tb.time*1000+tb.millitm;
#else
  static struct tms tb;
  finish = times(&tb);
#endif
    
  difference = finish - begin;
  begin = finish;
    
  return (float)difference/(float)CLOCKS_PER_SEC;
}

void
shadowtext(int x, int y, char* s) 
{
  int lines;
  char* p;
    
  glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, glutGet(GLUT_WINDOW_WIDTH), 
	  0, glutGet(GLUT_WINDOW_HEIGHT), -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glColor3ub(0, 0, 0);
  glRasterPos2i(x+1, y-1);
  for(p = s, lines = 0; *p; p++) {
    if (*p == '\n') {
      lines++;
      glRasterPos2i(x+1, y-1-(lines*18));
    }
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *p);
  }
  glColor3ub(0, 128, 255);
  glRasterPos2i(x, y);
  for(p = s, lines = 0; *p; p++) {
    if (*p == '\n') {
      lines++;
      glRasterPos2i(x, y-(lines*18));
    }
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *p);
  }
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glEnable(GL_DEPTH_TEST);
}

void lists(void)
{
  GLfloat ambient[] = { 0.2, 0.2, 0.2, 1.0 };
  GLfloat diffuse[] = { 0.8, 0.8, 0.8, 1.0 };
  GLfloat specular[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat shininess = 65.0;
  DIFF_COEFF = 0.8;
  SPEC_COEFF = 0.0;
  AMB_COEFF = 0.2;
  glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
  glMaterialf(GL_FRONT, GL_SHININESS, shininess);
    
  if (model_list)
    glDeleteLists(model_list, 1);
    
  /* generate a list */
  if (material_mode == 0) { 
    if (facet_normal)
      model_list = glmList(model, GLM_FLAT);
    else
      model_list = glmList(model, GLM_SMOOTH);
  } else if (material_mode == 1) {
    if (facet_normal)
      model_list = glmList(model, GLM_FLAT | GLM_COLOR);
    else
      model_list = glmList(model, GLM_SMOOTH | GLM_COLOR);
  } else if (material_mode == 2) {
    if (facet_normal)
      model_list = glmList(model, GLM_FLAT | GLM_MATERIAL);
    else
      model_list = glmList(model, GLM_SMOOTH | GLM_MATERIAL);
  }
}

void init(void)
{
  gltbInit(GLUT_LEFT_BUTTON);
    
  /* read in the model */
  model = glmReadOBJ(model_file);
  scale = glmUnitize(model);
  glmFacetNormals(model);
  glmVertexNormals(model, smoothing_angle);
    
  if (model->nummaterials > 0)
    material_mode = 2;
    
  /* create new display lists */
  lists();
    
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    
  glEnable(GL_DEPTH_TEST);
    
  glEnable(GL_CULL_FACE);
 
}

void reshape(int width, int height)
{
  gltbReshape(width, height);
    
  glViewport(0, 0, width, height);
    
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, (GLfloat)height / (GLfloat)width, 1.0, 128.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(0.0, 0.0, -3.0);
}

#define NUM_FRAMES 5
void display(void)
{
  static char s[256], t[32];
  static char* p;
  static int frames = 0;
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	    
  glPushMatrix();	    
  glTranslatef(pan_x, pan_y, 0.0);	    
  gltbMatrix();
  if(rasterFlag) {  
    rasterization();   
    glDrawPixels(512,512,GL_RGB,GL_FLOAT,pixels);    
  
  } else {
#if 0   /* glmDraw() performance test */
    if (material_mode == 0) { 
      if (facet_normal)
	glmDraw(model, GLM_FLAT);
      else
	glmDraw(model, GLM_SMOOTH);
    } else if (material_mode == 1) {
      if (facet_normal)
	glmDraw(model, GLM_FLAT | GLM_COLOR);
      else
	glmDraw(model, GLM_SMOOTH | GLM_COLOR);
    } else if (material_mode == 2) {
      if (facet_normal)
	glmDraw(model, GLM_FLAT | GLM_MATERIAL);
      else
	glmDraw(model, GLM_SMOOTH | GLM_MATERIAL);
    }
#else
    glCallList(model_list);
#endif
      
      
    glDisable(GL_LIGHTING);
    if (bounding_box) {
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_BLEND);
      glEnable(GL_CULL_FACE);
      glColor4f(1.0, 0.0, 0.0, 0.25);
      glutSolidCube(2.0);
      glDisable(GL_BLEND);
    }
  }	    
  glPopMatrix();
	    
  if (stats) {
    /* XXX - this could be done a _whole lot_ faster... */
    int height = glutGet(GLUT_WINDOW_HEIGHT);
    glColor3ub(0, 0, 0);
    sprintf(s, "%s\n%d vertices\n%d triangles\n%d normals\n"
	    "%d texcoords\n%d groups\n%d materials",
	    model->pathname, model->numvertices, model->numtriangles, 
	    model->numnormals, model->numtexcoords, model->numgroups,
	    model->nummaterials);
    shadowtext(5, height-(5+18*1), s);
  }
	    
  /* spit out frame rate. */
  frames++;
  if (frames > NUM_FRAMES) {
    sprintf(t, "%g fps", frames/elapsed());
    frames = 0;
  }
  if (performance) {
    shadowtext(5, 5, t); 
  }
	    
  glutSwapBuffers();
  glEnable(GL_LIGHTING);
}

void keyboard(unsigned char key, int x, int y)
{
  GLint params[2];
    
  switch (key) {
  case 'h':
    printf("help\n\n");
    printf("w         -  Toggle wireframe/filled\n");
    printf("c         -  Toggle culling\n");
    printf("n         -  Toggle facet/smooth normal\n");
    printf("b         -  Toggle bounding box\n");
    printf("r         -  Reverse polygon winding\n");
    printf("m         -  Toggle color/material/none mode\n");
    printf("p         -  Toggle performance indicator\n");
    printf("s/S       -  Scale model smaller/larger\n");
    printf("t         -  Show model stats\n");
    printf("o         -  Weld vertices in model\n");
    printf("+/-       -  Increase/decrease smoothing angle\n");
    printf("W         -  Write model to file (out.obj)\n");
    printf("q/escape  -  Quit\n\n");
    break;
        
  case 't':
    stats = !stats;
    break;
        
  case 'p':
    performance = !performance;
    break;
        
  case 'm':
    material_mode++;
    if (material_mode > 2)
      material_mode = 0;
    printf("material_mode = %d\n", material_mode);
    lists();
    break;
        
  case 'd':
    glmDelete(model);
    init();
    lists();
    break;
        
  case 'w':
    glGetIntegerv(GL_POLYGON_MODE, params);
    if (params[0] == GL_FILL)
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    break;
        
  case 'c':
    if (glIsEnabled(GL_CULL_FACE))
      glDisable(GL_CULL_FACE);
    else
      glEnable(GL_CULL_FACE);
    break;
        
  case 'b':
    bounding_box = !bounding_box;
    break;
        
  case 'n':
    facet_normal = !facet_normal;
    lists();
    break;
        
  case 'r':
    glmReverseWinding(model);
    lists();
    break;
        
  case 's':
    glmScale(model, 0.8);
    lists();
    break;
        
  case 'S':
    glmScale(model, 1.25);
    lists();
    break;
        
  case 'o':
    //printf("Welded %d\n", glmWeld(model, weld_distance));
    glmVertexNormals(model, smoothing_angle);
    lists();
    break;
        
  case 'O':
    weld_distance += 0.01;
    printf("Weld distance: %.2f\n", weld_distance);
    glmWeld(model, weld_distance);
    glmFacetNormals(model);
    glmVertexNormals(model, smoothing_angle);
    lists();
    break;
        
  case '-':
    smoothing_angle -= 1.0;
    printf("Smoothing angle: %.1f\n", smoothing_angle);
    glmVertexNormals(model, smoothing_angle);
    lists();
    break;
        
  case '+':
    smoothing_angle += 1.0;
    printf("Smoothing angle: %.1f\n", smoothing_angle);
    glmVertexNormals(model, smoothing_angle);
    lists();
    break;
        
  case 'W':
    glmScale(model, 1.0/scale);
    glmWriteOBJ(model, "out.obj", GLM_SMOOTH | GLM_MATERIAL);
    break;
        
  case 'R':
    {
      GLuint i;
      GLfloat swap;
      for (i = 1; i <= model->numvertices; i++) {
	swap = model->vertices[3 * i + 1];
	model->vertices[3 * i + 1] = model->vertices[3 * i + 2];
	model->vertices[3 * i + 2] = -swap;
      }
      glmFacetNormals(model);
      lists();
      break;
    }
  case 'y' : 
  case 'Y' : 
    rasterFlag = 1 - rasterFlag; 
    break;
  case 'f' :
  case 'F' :
    flatonly = 1 - flatonly;
    break;
  case 27:
    exit(0);
    break;
  }
    
  glutPostRedisplay();
}

void menu(int item)
{
  int i = 0;
  DIR* dirp;
  char* name;
  struct dirent* direntp;
    
  if (item > 0) {
    keyboard((unsigned char)item, 0, 0);
  } else {
    dirp = opendir(DATA_DIR);
    while ((direntp = readdir(dirp)) != NULL) {
      if (strstr(direntp->d_name, ".obj")) {
	i++;
	if (i == -item)
	  break;
      }
    }
    if (!direntp)
      return;
    name = (char*)malloc(strlen(direntp->d_name) + strlen(DATA_DIR) + 1);
    strcpy(name, DATA_DIR);
    strcat(name, direntp->d_name);
    model = glmReadOBJ(name);
    scale = glmUnitize(model);
    glmFacetNormals(model);
    glmVertexNormals(model, smoothing_angle);
        
    if (model->nummaterials > 0)
      material_mode = 2;
    else
      material_mode = 0;
        
    lists();
    free(name);
        
    glutPostRedisplay();
  }
}

static GLint      mouse_state;
static GLint      mouse_button;

void mouse(int button, int state, int x, int y)
{
  GLdouble model[4*4];
  GLdouble proj[4*4];
  GLint view[4];
    
  /* fix for two-button mice -- left mouse + shift = middle mouse */
  if (button == GLUT_LEFT_BUTTON && glutGetModifiers() & GLUT_ACTIVE_SHIFT)
    button = GLUT_MIDDLE_BUTTON;
    
  gltbMouse(button, state, x, y);
    
  mouse_state = state;
  mouse_button = button;
    
  if (state == GLUT_DOWN && button == GLUT_MIDDLE_BUTTON) {
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);
    glGetIntegerv(GL_VIEWPORT, view);
    gluProject((GLdouble)x, (GLdouble)y, 0.0,
	       model, proj, view,
	       &pan_x, &pan_y, &pan_z);
    gluUnProject((GLdouble)x, (GLdouble)y, pan_z,
		 model, proj, view,
		 &pan_x, &pan_y, &pan_z);
    pan_y = -pan_y;
  }
    
  glutPostRedisplay();
}

void motion(int x, int y)
{
  GLdouble model[4*4];
  GLdouble proj[4*4];
  GLint view[4];
    
  gltbMotion(x, y);
    
  if (mouse_state == GLUT_DOWN && mouse_button == GLUT_MIDDLE_BUTTON) {
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);
    glGetIntegerv(GL_VIEWPORT, view);
    gluProject((GLdouble)x, (GLdouble)y, 0.0,
	       model, proj, view,
	       &pan_x, &pan_y, &pan_z);
    gluUnProject((GLdouble)x, (GLdouble)y, pan_z,
		 model, proj, view,
		 &pan_x, &pan_y, &pan_z);
    pan_y = -pan_y;
  }
    
  glutPostRedisplay();
}

int main(int argc, char** argv)
{
  int buffering = GLUT_DOUBLE;
  struct dirent* direntp;
  DIR* dirp;
  int models;
    
  glutInitWindowSize(512, 512);
  glutInit(&argc, argv);
    
  while (--argc) {
    if (strcmp(argv[argc], "-sb") == 0)
      buffering = GLUT_SINGLE;
    else
      model_file = argv[argc];
  }
    
  if (!model_file) {
    model_file = "data/cube.obj";
  }
    
  glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | buffering);
  glutCreateWindow("Smooth");
    
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
    
  models = glutCreateMenu(menu);
  dirp = opendir(DATA_DIR);
  if (!dirp) {
    fprintf(stderr, "%s: can't open data directory.\n", argv[0]);
  } else {
    while ((direntp = readdir(dirp)) != NULL) {
      if (strstr(direntp->d_name, ".obj")) {
	entries++;
	glutAddMenuEntry(direntp->d_name, -entries);
      }
    }
    closedir(dirp);
  }
    
  glutCreateMenu(menu);
  glutAddMenuEntry("Smooth", 0);
  glutAddMenuEntry("", 0);
  glutAddSubMenu("Models", models);
  glutAddMenuEntry("", 0);
  glutAddMenuEntry("[w]   Toggle wireframe/filled", 'w');
  glutAddMenuEntry("[c]   Toggle culling on/off", 'c');
  glutAddMenuEntry("[n]   Toggle face/smooth normals", 'n');
  glutAddMenuEntry("[b]   Toggle bounding box on/off", 'b');
  glutAddMenuEntry("[p]   Toggle frame rate on/off", 'p');
  glutAddMenuEntry("[t]   Toggle model statistics", 't');
  glutAddMenuEntry("[m]   Toggle color/material/none mode", 'm');
  glutAddMenuEntry("[r]   Reverse polygon winding", 'r');
  glutAddMenuEntry("[s]   Scale model smaller", 's');
  glutAddMenuEntry("[S]   Scale model larger", 'S');
  glutAddMenuEntry("[o]   Weld redundant vertices", 'o');
  glutAddMenuEntry("[+]   Increase smoothing angle", '+');
  glutAddMenuEntry("[-]   Decrease smoothing angle", '-');
  glutAddMenuEntry("[W]   Write model to file (out.obj)", 'W');
  glutAddMenuEntry("", 0);
  glutAddMenuEntry("[Esc] Quit", 27);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
    
  init();
   
  glutMainLoop();
  return 0;
}
