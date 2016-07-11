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
extern map<string, GLfloat> vs;
extern map<string, map<string, map<string, map<string, map<string, map<string, GLfloat> > > > > > vp;

//extern GLfloat max_distance_hydrogen_bond;
extern GLint forca_externa_contador_max_hb;

extern GLint forca_externa_contador_max;
GLint forca_externa_contador_max_best = 0;

extern GLint forca_externa_contador_max_t;
GLint forca_externa_contador_max_best_t = 0;

extern GLint calibration_minimal_distance;

extern GLint contador_restart_life;

GLint max_life = 0;
GLfloat best_valueF = 0;
GLint best_valueI = 0;

extern void sV(string s1, string s2, GLfloat valor);
extern GLfloat gV(string s1, string s2);

extern map<string, map<string, GLfloat*> > vF;
extern map<string, map<string, GLint*> > vI;

void change_properties() {
//	printf("V1==> %f\n",*vF["tuning"]["max_distance_hydrogen_bond"]);
//	if(contador_restart_life > max_life){
//		best_valueF = *vF["tuning"]["max_distance_hydrogen_bond"];
//		max_life = contador_restart_life;
//	}
//	*vF["tuning"]["max_distance_hydrogen_bond"] += 0.1;
//	printf("V2==> %f (%f / %d)\n",*vF["tuning"]["max_distance_hydrogen_bond"], best_valueF, max_life);


//	printf("V1==> %d\n", *vI["tuning"]["forca_externa_contador_max"]);
	if (contador_restart_life > max_life) {
		best_valueI = *vI["tuning"]["forca_externa_contador_max"];
		max_life = contador_restart_life;
	}
	*vI["tuning"]["forca_externa_contador_max"] += 1;
	printf("V2==> %d (%d / %d)\n", *vI["tuning"]["forca_externa_contador_max"], best_valueI, max_life);

//	for (map<string, map<string, GLfloat*> >::iterator i = vF.begin(); i != vF.end(); ++i) {
//			for (map<string, GLfloat*>::iterator ii = i->second.begin(); ii != i->second.end(); ++ii) {
//				printf("Variavel %s %f\n",(*ii).first,(*ii).second);
////				string campochave = (*ii).first;
////				tw[campochave] = Convert((*ii).second);
//			}
//	}

// Procurando valores pra tracao e colisao ------------------------------------------------
//	if(contador_restart_life > max_contador){
//		printf("------------- Encontrado menor ------------- %d %d\n",contador_restart_life,max_contador);
//		forca_externa_contador_max_best = forca_externa_contador_max;
//		forca_externa_contador_max_best_t = forca_externa_contador_max_t;
//		max_contador = contador_restart_life;
//	}
//	forca_externa_contador_max++;
//	if(forca_externa_contador_max > 30){
//		forca_externa_contador_max = 0;
//		forca_externa_contador_max_t++;
//		if(forca_externa_contador_max_t > 30){
//			forca_externa_contador_max_t = 0;
//		}
//	}
//	printf("forca_externa_contador_max: %d %f: %d best: %d bestt: %d max_contador_best: %d\n",forca_externa_contador_max,forca_externa_contador_max_t,forca_externa_contador_max_best,forca_externa_contador_max_best_t,max_contador);
	// ----------------------------------------------------------------------------------------
}
