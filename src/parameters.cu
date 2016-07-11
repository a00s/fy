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
#include "gvariaveis.h"

using namespace std;
GLint menu_page = 0;

extern map<string, map<string, GLfloat*> > vF;
extern map<string, map<string, GLint*> > vI;
extern map<string, map<string, bool*> > vB;

//extern map<string, GLfloat> vs;
//extern map<string, map<string, map<string, map<string, map<string, map<string, GLfloat> > > > > > vp;

int contadorll = 0;
map<string, string> tw;

// Variavel de sistema
//void setVS(string s1, GLfloat valor) {
//	vs[s1] = valor;
//}

//GLfloat getVS(string s1) {
//	return 0.0;
//	return vs[s1];
//}

void parameters_window() {

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0, windowWidth, 0.0, windowHeight);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	if (menu_page == 0) {
		glColor3f(0.0, 1.0, 0.0);
	} else if (menu_page == 1) {
		glColor3f(1.0, 1.0, 0.0);
	} else if (menu_page == 2) {
		glColor3f(0.0, 1.0, 1.0);
	} else if (menu_page == 3) {
		glColor3f(1.0, 1.0, 1.0);
	} else if (menu_page == 4) {
		glColor3f(0.0f, 1.0f, 0.0f);
	} else {
		glColor3f(0.0, 1.0, 0.0);				// Green
	}
	GLint distancia_topo = windowHeight;
	string texto = "";
	void * font = GLUT_BITMAP_9_BY_15;

	typedef map<string, string>::iterator it_type;
	for (it_type iterator = tw.begin(); iterator != tw.end(); iterator++) {
		distancia_topo -= 14;
		glRasterPos2i(10, distancia_topo);
		texto = iterator->first + " = " + iterator->second;
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

//string Convert(GLfloat number) {
//	ostringstream buff;
//	buff << number;
//	return buff.str();
//}

string Convert(GLfloat *number) {
	ostringstream buff;
	buff << *number;
	return buff.str();
}

string Convert(GLint *number) {
	ostringstream buff;
	buff << *number;
	return buff.str();
}

string Convert(bool *number) {
	ostringstream buff;
	buff << *number;
	return buff.str();
}

void imprime_menu(string tipomenu) {
	tw.clear();
//	tw[""] = "Variaveis de Screen";
	for (map<string, map<string, GLfloat*> >::iterator i = vF.begin(); i != vF.end(); ++i) {
		if ((*i).first == tipomenu) {
			for (map<string, GLfloat*>::iterator ii = i->second.begin(); ii != i->second.end(); ++ii) {
				string campochave = (*ii).first;
				tw[campochave] = Convert((*ii).second);
			}
		}
	}
	for (map<string, map<string, GLint*> >::iterator i = vI.begin(); i != vI.end(); ++i) {
		if ((*i).first == tipomenu) {
			for (map<string, GLint*>::iterator ii = i->second.begin(); ii != i->second.end(); ++ii) {
				string campochave = (*ii).first;
				tw[campochave] = Convert((*ii).second);
			}
		}
	}
	for (map<string, map<string, bool*> >::iterator i = vB.begin(); i != vB.end(); ++i) {
		if ((*i).first == tipomenu) {
			for (map<string, bool*>::iterator ii = i->second.begin(); ii != i->second.end(); ++ii) {
				string campochave = (*ii).first;
				tw[campochave] = Convert((*ii).second);
				if (Convert((*ii).second) == "1") {
					tw[campochave] = "true";
				} else {
					tw[campochave] = "false";
				}
			}
		}
	}
}

void handleKeypress(int theKey, int theAction) {
	if (theAction == GLFW_PRESS) {
		// ------------- Padrao ------------------
		switch (theKey) {
		case 289:
			pressionando_control = true;
			break;
		case '0':
			menu_page = 0;
			break;
		case '1':
			menu_page = 1;
			break;
		case '2':
			menu_page = 2;
			break;
		case '3':
			menu_page = 3;
			break;
		case '4':
			menu_page = 4;
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
		default:
			break;
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
		}
	}

	if (menu_page == 0) {
		tw.clear();
		tw["0"] = "Menu";
		tw["1"] = "Navegacao";
		tw["2"] = "Alteracoes";
		tw["3"] = "Variaveis de tela";
		tw["4"] = "Variaveis";
		tw["-"] = "-------------";
		tw["W"] = "Sobe";
		tw["S"] = "Desce";
		tw["A"] = "Esquerda";
		tw["D"] = "Direita";
	} else if (menu_page == 1) {
		// Navegacao
		tw.clear();
		tw["P"] = "Restaura posicoes";
		tw["C"] = "Compara Protein build MD";
		tw["V"] = "A/D comparacao speed";
		tw["X"] = "A/D comparacao";
		tw["E"] = "Camera 1";
		tw["R"] = "Camera 2";
		tw["T"] = "Camera 3";
		tw["Y"] = "Camera 4";
		tw[","] = "Diminui caixa";
		tw["."] = "Aumenta caixa";
		tw["B"] = "A/D base";
		tw["N"] = "A/D linha base";
		tw["M"] = "A/D perimetro";
		tw["J"] = "A/D rastreio";
		tw["F"] = "A/D forca";
		tw["G"] = "A/D ghost protein";
		if (theAction == GLFW_PRESS) {
			switch (theKey) {
			case 'P':
				restaura_posicoes();
				break;
			case 'C':
				compare_protein_build_MD(0);
				break;
			case 'V':
				ativa_desativa_comparation_speed();
				break;
			case 'X':
				ativa_desativa_comparation();
				break;
			case 'E':
				camera_position(0.0, 0.0, -15.0, 0.0, 0.0, 0.0);
				break;
			case 'R':
				camera_position(-40.0, 0.0, -50.0, 0.0, 90.0, 0.0);
				break;
			case 'T':
				camera_position(0.0, 40.0, -50.0, 90.0, 0.0, 0.0);
				break;
			case 'Y':
				camera_position(caixa_tamanho * 3.0, caixa_tamanho * 3.0, caixa_tamanho, 27.0, -30.0, 0.0);
				break;
			case ',':
				caixa_tamanho -= 0.5;
				break;
			case '.':
				caixa_tamanho += 0.5;
				break;
			case 'B':
				ativa_desativa_base();
				break;
			case 'N':
				ativa_desativa_base_line();
				break;
			case 'M':
				ativa_desativa_perimetro();
				break;
			case 'J':
				ativa_desativa_rastreio();
				break;
			case 'F':
				ativa_desativa_forca();
				break;
			case 'G':
				ativa_desativa_ghost_protein();
				break;
			case 'I':
				angulo_adicional_teste -= 1.0;
//				posx[36] -= 0.1;
				break;
			case 'O':
				angulo_adicional_teste += 1.0;
//				posx[36] += 0.1;
				break;
			case '9':
				posy[36] += 0.1;
				break;
			case 'K':
				posy[36] -= 0.1;
				break;
			case 'L':
				posz[36] -= 0.1;
				break;
			case ';':
				posz[36] += 0.1;
				break;
			case '[':
				sequencial_mostra--;
				break;
			case ']':
				sequencial_mostra++;
				break;
			}
		} else {

		}
	} else if (menu_page == 2) {
		// Alteracoes
		tw.clear();
		tw["T"] = "A/D tensao de colisao";
		tw["H"] = "A/D tensao de HB";
		tw["."] = "+ tune distancia";
		tw[","] = "- tune distancia";
		tw["Z"] = "+ energia";
		tw["X"] = "- energia";
		tw["P"] = "Pausa/Continua";
		tw["G"] = "Adiciona aminoacido";
		tw["K"] = "+X";
		tw["J"] = "-X";
		tw["I"] = "+Y";
		tw["M"] = "-Y";

		if (theAction == GLFW_PRESS) {
			switch (theKey) {
			case 'T':
				ativa_desativa_colisao_tensao();
				break;
			case 'H':
				ativa_desativa_tensao_hb();
				break;
			case '.':
				calibration_precision++;
				distance_calibration();
				break;
			case ',':
				calibration_precision--;
				distance_calibration();
				break;
			case 'Z':
				rem_energy();
				break;
			case 'X':
				add_energy();
				break;
			case 'P':
				continua();
				break;
			case 'G':
				add_chain();
				break;
				// X+
			case 'K':
				pressionando_k = false;
				break;
				// X-
			case 'J':
				pressionando_j = false;
				break;
				// Y+
			case 'I':
				pressionando_i = false;
				break;
				// Y-
			case 'M':
				pressionando_m = false;
				break;
			}
		} else {

		}
	} else if (menu_page == 3) {
		imprime_menu("screen");
	} else if (menu_page == 4) {
		// Variaveis
		imprime_menu("tuning");
//		tw.clear();
//		tw[""] = "Propriedades";

//		for (map<string, map<string, map<string, map<string, map<string, map<string, GLfloat> > > > > >::iterator i = vp.begin(); i != vp.end(); ++i) {
//			if ((*i).first != "PDBcalibrationMin" && (*i).first != "PDBcalibrationMax") {
//				for (map<string, map<string, map<string, map<string, map<string, GLfloat> > > > >::iterator ii = i->second.begin(); ii != i->second.end(); ++ii) {
//					for (map<string, map<string, map<string, map<string, GLfloat> > > >::iterator iii = ii->second.begin(); iii != ii->second.end(); ++iii) {
//						for (map<string, map<string, map<string, GLfloat> > >::iterator iiii = iii->second.begin(); iiii != iii->second.end(); ++iiii) {
//							for (map<string, map<string, GLfloat> >::iterator iiiii = iiii->second.begin(); iiiii != iiii->second.end(); ++iiiii) {
//								for (map<string, GLfloat>::iterator iiiiii = iiiii->second.begin(); iiiiii != iiiii->second.end(); ++iiiiii) {
//									string campochave = (*i).first + ' ' + (*ii).first + ' ' + (*iii).first + ' ' + (*iiii).first + ' ' + (*iiiii).first + ' ' + (*iiiiii).first;
//									tw[campochave] = Convert((*iiiiii).second);
//								}
//							}
//						}
//					}
//				}
//			}
//		}
	}
//	if (pressionando_control) {
//		switch (theKey) {
//		case 'C':
//			add_atom('C');
//			break;
//
//		case 'H':
//			add_atom('H');
//			break;
//
//		case 'O':
//			add_atom('O');
//			break;
//
//		case 'N':
//			add_atom('N');
//			break;
//

}
