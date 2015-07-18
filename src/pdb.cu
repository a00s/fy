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
using namespace std;

vector<GLfloat> pdb_coordenadas_x;
vector<GLfloat> pdb_coordenadas_y;
vector<GLfloat> pdb_coordenadas_z;
GLint contador = 0;

void read_pdb() {
	printf("Lendo arquivo PDB\n");
	string line = "";
	ifstream myfile("files/Tunnel_4UG0.pdb");
	if (myfile.is_open()) {
		string tipo;
		GLfloat fix_x = 0;
		GLfloat fix_y = 0;
		GLfloat fix_z = 0;
		while (getline(myfile, line)) {
			tipo = line.substr(0, 6);
			if (tipo.erase(tipo.find_last_not_of(" ") + 1) == "ATOM" || tipo == "HETATM") {
				contador++;
				pdb_coordenadas_x.resize(contador);
				pdb_coordenadas_y.resize(contador);
				pdb_coordenadas_z.resize(contador);

				pdb_coordenadas_x[contador - 1] = strtof(line.substr(30, 8).c_str(), NULL) - fix_x;
				pdb_coordenadas_y[contador - 1] = strtof(line.substr(38, 8).c_str(), NULL) - fix_y;
				pdb_coordenadas_z[contador - 1] = strtof(line.substr(46, 8).c_str(), NULL) - fix_z;

				if (contador == 1) {
					fix_x = pdb_coordenadas_x[contador - 1];
					fix_y = pdb_coordenadas_y[contador - 1];
					fix_z = pdb_coordenadas_z[contador - 1];
					pdb_coordenadas_x[contador - 1] = 0.0;
					pdb_coordenadas_y[contador - 1] = 0.0;
					pdb_coordenadas_z[contador - 1] = 0.0;
				}

//				GLfloat cart_x = strtof(line.substr(30, 8).c_str(), NULL);
//				GLfloat cart_y = strtof(line.substr(38, 8).c_str(), NULL);
//				GLfloat cart_z = strtof(line.substr(46, 8).c_str(), NULL);
////				printf("%s\n",line.c_str());
//				printf("%f %f %f\n", cart_x, cart_y, cart_z);

//				glColor3ub(255, 255, 255);
//				glTranslatef(0.0, 0.0, 0.0);
//				glutSolidSphere(0.1f, 6, 6);
//				glTranslatef(0.0, 0.0, 0.0);
			}
		}
		myfile.close();
	} else {
		cout << "Unable to open file\n";
	}
}

void add_pdb() {
	for (int i = 0; i < contador; i++) {
		glColor3ub(255, 255, 255);
		glTranslatef(pdb_coordenadas_x[i], pdb_coordenadas_y[i], pdb_coordenadas_z[i]);
		glutSolidCube(0.1f);
		glTranslatef(-pdb_coordenadas_x[i], -pdb_coordenadas_y[i], -pdb_coordenadas_z[i]);
	}
}

void grava_pdb() {
	printf("Lendo arquivo PDB\n");
	//Grava
	ofstream myfile("example.txt");
	if (myfile.is_open()) {
		myfile << "This is a line.\n";
		myfile << "This is another line.\n";
		myfile.close();
	} else
		cout << "Unable to open file";
}
