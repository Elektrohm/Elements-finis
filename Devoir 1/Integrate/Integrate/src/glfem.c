/*
 *  glfem.c
 *  Library for EPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  GLFW  http://www.glfw.org/ (version utilisée 3.3.2)
 *  BOV   https://git.immc.ucl.ac.be/hextreme/NGP/-/tree/master/deps/BOV
 *
 */
 
 
#include "glfem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



double   zoom_init;
double   translate_init[2];
float    current_text_color[4] = {1.0,0.0,0.0,1.0};
float    current_color[4]      = {0.0,0.0,0.0,1.0};


typedef bov_window_t glfemWindow;
glfemWindow* theCurrentWindow = NULL;


// ======================================================================================
// ======================= Generic Functions for homework ===============================
// ======================================================================================

void glfemSetColor(float color[4]) 
{
    for(int i = 0; i < 4; ++i) {
        current_color[i] = color[i]; }
}

void glfemSetTextColor(float color[4]) 
{
    for(int i = 0; i < 4; ++i) {
        current_text_color[i] = color[i]; }
}



static void glfemKeyCallback(GLFWwindow* self,
                         int key, int scancode, int action,int mods) {
  bov_window_t* window = (bov_window_t*) glfwGetWindowUserPointer(self);
  if(action==GLFW_PRESS || action==GLFW_REPEAT) {
	switch(key) {
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose(self,GL_TRUE);
			break;
		case GLFW_KEY_H:
		case GLFW_KEY_K:
			if(window->help_needed==0) {
				window->help_needed = 1;
			}
			else {
				window->help_needed = 0;
			}
			break;
		case GLFW_KEY_R:
			window->param.zoom = zoom_init;
			window->param.translate[0] = translate_init[0];
			window->param.translate[1] = translate_init[1];
			break;
	}
  }
  if(key==GLFW_KEY_ESCAPE)
	glfwSetWindowShouldClose(self,GL_TRUE);
}
	

void glfemWindowCreate(const char *windowName,int w, int h,int n,double *x, double *y){
  // Initialize Window
  bov_window_t *window = bov_window_new(w,h, windowName);
  theCurrentWindow = window;
  bov_window_set_color(window, (GLfloat[4]){0.9, 0.9, 0.8, 0.0});
  
  // Set Viewport and Initial Zoom & Translate based on nodes
  double minX  = femMin(x,n);
  double maxX  = femMax(x,n);
  double minY  = femMin(y,n);
  double maxY  = femMax(y,n);
  double sizeX = (maxX-minX)/1.45;
  double meanX = (maxX+minX)/2.0; 
  double sizeY = (maxY-minY)/1.45;
  double meanY = (maxY+minY)/2.0;
    
  double ratio = (GLfloat) h / (GLfloat) w;
  double size = fmax(sizeX,sizeY);
  double left,right,top,bottom;
  if (ratio > 1.0) {
      left = meanX - size;
      right = meanX + size;
      bottom = meanY - size*ratio;
      top = meanY + size*ratio;  }   
  else {
      left = meanX - size/ratio;
      right = meanX + size/ratio;
      bottom = meanY - size;
      top = meanY + size;  }   
      
  if ((fabs(top) + fabs(bottom)) > (fabs(left) + fabs(right))) {
      window->param.zoom = 1/(0.5*(fabs(top) + fabs(bottom))); }   
  else {
      window->param.zoom = 1/(0.5*(fabs(left) + fabs(right))); }   
  
  window->param.translate[0] = -(left + right) / 2.0;
  window->param.translate[1] = -(bottom + top) / 2.0;
  
  zoom_init = window->param.zoom;
  translate_init[0] = window->param.translate[0];
  translate_init[1] = window->param.translate[1];
  
  // Set KeyCallBack and Help menu
  glfwSetKeyCallback(window->self, glfemKeyCallback);
 	
  window->help = bov_text_new((GLubyte[]) {
	"   [esc]   exit\n"
	"     r     reset zoom and translation\n"
	"    h k    display/hide keyboard shortcuts\n"
  }, GL_STATIC_DRAW);
  bov_text_set_space_type(window->help, PIXEL_SPACE);
  bov_text_set_fontsize(window->help, 20.0f); 
  bov_text_set_boldness(window->help, 0.1f);
  bov_text_set_outline_width(window->help, 0.5f);
  bov_text_set_color(window->help,current_text_color);
  
}



void glfemWindowUpdate() 
{
    bov_text_set_pos(theCurrentWindow->help, 
            (GLfloat[2]){-20.0f, theCurrentWindow->size[1] - 30.0f} );
    bov_window_update(theCurrentWindow);
}



void glfemWindowFree() 
{
    bov_window_delete(theCurrentWindow);
}



int glfemWindowShouldClose() 
{
    return bov_window_should_close(theCurrentWindow);
}




void glfemDrawMessage(char *message, double pos[2]){

  
  bov_text_t* text = bov_text_new((const GLubyte *)message, GL_STATIC_DRAW);
  text->param =  (bov_text_param_t) {
		.fillColor = {current_text_color[0],current_text_color[1],current_text_color[2],current_text_color[3]},
		.outlineColor = {1.0f ,1.0f, 1.0f, 2.0f},
		.pos = {0.0f, 0.0f},
		.shift = {0.0f, 0.0f},
		.fontSize = 20.0f,
		.boldness = 0.0f,
		.outlineWidth = 0.0f,
		.spaceType = PIXEL_SPACE};
  bov_text_set_pos(text, (GLfloat[2]) {pos[0], pos[1]});
  
  bov_text_draw(theCurrentWindow, text);
  bov_text_delete(text);
}


void glfemDrawNodes(double *x, double *y, int n) {
  GLfloat (*coord)[2] = malloc(sizeof(coord[0])*n);
  for(int i = 0; i < n; ++i){
    coord[i][0] = x[i];
    coord[i][1] = y[i];
  }
  bov_points_t* points = bov_points_new(coord,n,GL_STATIC_DRAW);
  bov_points_set_color(points,current_color);
  bov_points_set_width(points,0.0125/zoom_init);
  bov_points_draw(theCurrentWindow,points, 0, BOV_TILL_END);
  bov_points_delete(points);
  free(coord);
}

void glfemDrawElement(double *x, double *y, int n){
  GLfloat (*coord)[2] = malloc(sizeof(coord[0])*n);
  for(int i = 0; i < n; ++i){
    coord[i][0] = x[i];
    coord[i][1] = y[i];
  }
  bov_points_t* points = bov_points_new(coord,n,GL_STATIC_DRAW);
  bov_points_set_color(points,current_color);
  bov_points_set_width(points,0.0025/zoom_init);

  bov_line_loop_draw(theCurrentWindow, points, 0, BOV_TILL_END);
  bov_points_delete(points);
  free(coord);
}


