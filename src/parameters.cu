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

int contadorll = 0;
//extern bool mousefree;
extern GLint windowWidth;
extern GLint windowHeight;
///extern string fps;
//map<int, map<string, map<string, GLfloat> > > atom_statistic;
map<string, string> tw;

//class GUI: public QMainWindow {
//public:
//	GUI();
//}
void displayMe(void) {
	contadorll++;
	printf("A %d\n", contadorll);
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_POLYGON);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.5, 0.0, 0.0);
	glVertex3f(0.5, 0.5, 0.0);

//	glMatrixMode(GL_PROJECTION);
//	glPushMatrix();
//	glLoadIdentity();
//	gluOrtho2D(0.0, 300, 0.0, 300);
//	glMatrixMode(GL_MODELVIEW);
//	glPushMatrix();
//	glLoadIdentity();
//	glColor3f(0.0, 1.0, 0.0); // Green
//	glRasterPos2i(10, 10);
//	void * font = GLUT_BITMAP_9_BY_15;
//	stringstream ss;
//	ss << contadorll;
//	string texto = ss.str();
//	for (string::iterator i = texto.begin(); i != texto.end(); ++i) {
//		char c = *i;
//		glutBitmapCharacter(font, c);
//	}
//	glMatrixMode(GL_MODELVIEW);
//	glPopMatrix();
//
//	glMatrixMode(GL_PROJECTION);
//	glPopMatrix();
//	// ----- Stop Drawing Stuff! ------
//	glfwSwapBuffers();

	glEnd();
	glFlush();
}

void help(){
	tw["Campo1"] = "Teste1";
	tw["Campo2"] = "Teste2";
	tw["Campo3"] = "Teste3";
	tw["Campo4"] = "Teste4";
}


void parameters_window() {
	help();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0, windowWidth, 0.0, windowHeight);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glColor3f(0.0, 1.0, 0.0);				// Green
	GLint distancia_topo = windowHeight;
	string texto = "";
	void * font = GLUT_BITMAP_9_BY_15;
	typedef map<string, string>::iterator it_type;
	for(it_type iterator = tw.begin(); iterator != tw.end(); iterator++) {
		distancia_topo -= 14;
	    glRasterPos2i(10, distancia_topo);
		texto = iterator->first+" = "+iterator->second;
		for (string::iterator i = texto.begin(); i != texto.end(); ++i) {
			char c = *i;
			glutBitmapCharacter(font, c);
		}
	}
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
}

void handleKeypress(int theKey, int theAction) {
/*
	if (theAction == GLFW_PRESS) {
		if (pressionando_control) {
			switch (theKey) {
			case 'C':
				add_atom('C');
				break;

			case 'H':
				add_atom('H');
				break;

			case 'O':
				add_atom('O');
				break;

			case 'N':
				add_atom('N');
				break;

			case '1':
				ativa_desativa_colisao_tensao();
				break;

			case '2':
				ativa_desativa_tensao_hb();
				break;

			case '3':
				compare_protein_build_MD(0);
				break;

			case '4':
				ativa_desativa_comparation();
				break;

			case '5':
				restaura_posicoes();
				break;

			case '6':
				ativa_desativa_comparation_speed();
				break;

			case '.':
				calibration_precision++;
				distance_calibration();
				break;

			case ',':
				calibration_precision--;
				distance_calibration();
				break;
			}
		} else {
			switch (theKey) {
			case 289:
				pressionando_control = true;
				break;

			case 'W':
				holdingForward = true;
				break;

			case 'S':
				holdingBackward = true;
				break;

			case 'A':
				holdingLeftStrafe = true;
				break;

			case 'D':
				holdingRightStrafe = true;
				break;

			case 'K':
				pressionando_k = true;
				break;

			case 'J':
				pressionando_j = true;
				break;

			case 'I':
				pressionando_i = true;
				break;

			case 'M':
				pressionando_m = true;
				break;

			case 'Y':
				posz[0] -= 0.5;
				break;

			case 'H':
				posz[0] += 0.5;
				break;

			case 'O':
				pressionando_o = true;
				break;

			case 'P':
				pressionando_p = true;
				break;

			case '0':
				pressionando_0 = true;
				break;

			case 'L':
				pressionando_l = true;
				break;

			case '8':
				pressionando_8 = true;
				break;

			case '9':
				pressionando_9 = true;
				break;

			case 'B':
				show_variables();
				break;

			case 'G':
				add_chain();
				break;

			case '1':
				camera_position(0.0, 0.0, caixa_tamanho, 0.0, 0.0, 0.0);
				break;

			case '2':
				camera_position(caixa_tamanho * 6.0, 0.0, -90.0, 0.0, -90.0, 0.0);
				break;

			case '3':
				camera_position(0.0, caixa_tamanho * 6.0, -90.0, 90.0, 180.0, -90.0);
				break;

			case '4':
				camera_position(caixa_tamanho * 3.0, caixa_tamanho * 3.0, caixa_tamanho, 27.0, -30.0, 0.0);
				break;

			case 'E':
				ativa_desativa_perimetro();
				break;

			case 'R':
				ativa_desativa_rastreio();
				break;

			case 'F':
				ativa_desativa_forca();
				break;

			case 'Z':
				rem_energy();
				break;

			case 'X':
				add_energy();
				break;

			case 'C':
				continua();
				break;

			case ',':
				caixa_tamanho -= 0.5;
				break;

			case '.':
				caixa_tamanho += 0.5;
				break;

			case 'T':
				ativa_desativa_base();
				break;

			case '5':
				ativa_desativa_base_line();
				break;

			default:
				// Do nothing...
				break;
			}
		}
	} else // If a key is released, toggle the relevant key-release flag
	{
		switch (theKey) {
		case 289:
			pressionando_control = false;
			break;

		case 'W':
			holdingForward = false;
			break;

		case 'S':
			holdingBackward = false;
			break;

		case 'A':
			holdingLeftStrafe = false;
			break;

		case 'D':
			holdingRightStrafe = false;
			break;

		case 'K':
			pressionando_k = false;
			break;

		case 'J':
			pressionando_j = false;
			break;

		case 'I':
			pressionando_i = false;
			break;

		case 'M':
			pressionando_m = false;
			break;

		case 'O':
			pressionando_o = false;
			break;

		case 'P':
			pressionando_p = false;
			break;

		case '0':
			pressionando_0 = false;
			break;

		case 'L':
			pressionando_l = false;
			break;

		case '8':
			pressionando_8 = false;
			break;

		case '9':
			pressionando_9 = false;
			break;

		default:
			// Do nothing...
			break;
		}
	}
	*/
}
