#include <iostream>
#include <GL/GLee.h>         // No need to link to GL/gl.h
#include <GL/glfw.h>      // Include OpenGL Framework library
#include <GL/freeglut.h>  // Include FreeGLUT so we can easily draw spheres and calculate our viewing frustrum
#include <math.h>         // Used only for sin() and cos() functions
#include <cstdio>
#include <stdlib.h>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <vector>
#include <map>
using namespace std;

GLint max_contador = 0;
extern GLfloat max_distance_hydrogen_bond;
extern GLint forca_externa_contador_max_hb;

extern GLint forca_externa_contador_max;
GLint forca_externa_contador_max_best = 0;

extern GLint forca_externa_contador_max_t;
GLint forca_externa_contador_max_best_t = 0;

extern GLint calibration_minimal_distance;

extern GLint contador_restart_life;

void change_properties(){

	// Procurando valores pra tracao e colisao ------------------------------------------------
	if(contador_restart_life > max_contador){
		printf("------------- Encontrado menor ------------- %d %d\n",contador_restart_life,max_contador);
		forca_externa_contador_max_best = forca_externa_contador_max;
		forca_externa_contador_max_best_t = forca_externa_contador_max_t;
		max_contador = contador_restart_life;
	}
	forca_externa_contador_max++;
	if(forca_externa_contador_max > 30){
		forca_externa_contador_max = 0;
		forca_externa_contador_max_t++;
		if(forca_externa_contador_max_t > 30){
				forca_externa_contador_max_t = 0;
		}
	}
	printf("forca_externa_contador_max: %d t: %d best: %d bestt: %d max_contador_best: %d\n",forca_externa_contador_max,forca_externa_contador_max_t,forca_externa_contador_max_best,forca_externa_contador_max_best_t,max_contador);
	// ----------------------------------------------------------------------------------------
}
