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
#include "mysql_connection.h"
#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>
#include <cppconn/prepared_statement.h>
using namespace std;

vector<GLfloat> pdbamin_coordenadas_x;
vector<GLfloat> pdbamin_coordenadas_y;
vector<GLfloat> pdbamin_coordenadas_z;

struct VectorProtein {
	GLint aminoseq;
	string label;
	GLfloat x;
	GLfloat y;
	GLfloat z;
};
vector<VectorProtein> pdb_analise;

struct StructAmino {
	GLfloat x;
	GLfloat y;
	GLfloat z;
};
vector<StructAmino> AminoTemp1;
vector<StructAmino> AminoTemp2;

vector<StructAmino> AminoTemp1MD;
vector<StructAmino> AminoTemp2MD;

GLint cont_amin = 0;
extern int atomos_quantidade;
extern char atomo_letra[1000];
extern GLfloat nucleo_proximity[1000];
extern GLfloat nucleo_proximity_free[1000];
extern GLfloat nucleo_proximity_HB[1000];
extern GLint electron_quantidade[1000];
extern GLfloat electron_raio[1000];
extern GLfloat electron_raio_HB[1000];
extern GLint electron_arested[1000][4][2];
extern GLfloat massa[1000];
extern GLfloat posx[1000];
extern GLfloat posy[1000];
extern GLfloat posz[1000];
extern GLfloat velocidade_x[1000]; // velocidade
extern GLfloat velocidade_y[1000]; // velocidade
extern GLfloat velocidade_z[1000]; // velocidade
extern GLfloat electron_y[1000][4];
extern GLfloat electron_z[1000][4];
extern GLint amino[1000];
extern GLint amino_sequencial[1000];
extern GLint atomo_letraN[1000];
extern GLint contador_amino;
extern int get_amino_number(const char *amino_sigla);
extern int get_atom_number(const char *atomo_letra_local);
extern bool atomo_base[1000];
extern GLint chain_ultimo_atomo;
map<string, int> mapa_atomo;
map<string, int> electron_contador;

void connect_electron(string atomo1, string atomo2) {
	electron_arested[mapa_atomo[atomo1]][electron_contador[atomo1]][0] = mapa_atomo[atomo2];
	electron_arested[mapa_atomo[atomo1]][electron_contador[atomo1]][1] = electron_contador[atomo2];
	electron_arested[mapa_atomo[atomo2]][electron_contador[atomo2]][0] = mapa_atomo[atomo1];
	electron_arested[mapa_atomo[atomo2]][electron_contador[atomo2]][1] = electron_contador[atomo1];
	electron_contador[atomo1]++;
	electron_contador[atomo2]++;
}

void add_atoml(string atomo_ll) {
	if (atomos_quantidade == 998) {
		printf("Quantidade de atomos chegou no limite \n");
		return;
	}

	if (atomo_ll == "C") {
		atomo_letra[atomos_quantidade] = 'C';
		nucleo_proximity[atomos_quantidade] = 1.7;
		nucleo_proximity_free[atomos_quantidade] = 3.4; // Van der waals angstron
		electron_quantidade[atomos_quantidade] = 4;
		electron_raio[atomos_quantidade] = 3.4;
		massa[atomos_quantidade] = 12.0107;
	} else if (atomo_ll == "H") {
		atomo_letra[atomos_quantidade] = 'H';
		nucleo_proximity[atomos_quantidade] = 1.2;
		nucleo_proximity_free[atomos_quantidade] = 2.4;  // Van der waals angstron
		electron_quantidade[atomos_quantidade] = 1;
		electron_raio[atomos_quantidade] = 2.4;
		electron_raio_HB[atomos_quantidade] = 2.4;
		massa[atomos_quantidade] = 1.0079;
		nucleo_proximity_HB[atomos_quantidade] = 1.2;

//		electron_y_HB[atomos_quantidade] = 270;
	} else if (atomo_ll == "O") {
		atomo_letra[atomos_quantidade] = 'O';
		nucleo_proximity[atomos_quantidade] = 1.52;
		nucleo_proximity_free[atomos_quantidade] = 3.04;  // Van der waals angstron
		electron_quantidade[atomos_quantidade] = 2;
		electron_raio[atomos_quantidade] = 3.04;
		electron_raio_HB[atomos_quantidade] = 3.04;
		massa[atomos_quantidade] = 15.099;
		nucleo_proximity_HB[atomos_quantidade] = 1.52;

//		electron_y_HB[atomos_quantidade] = 90;
	} else if (atomo_ll == "N") {
		atomo_letra[atomos_quantidade] = 'N';
//		nucleo_proximity[atomos_quantidade] = 1.55; // original
		nucleo_proximity[atomos_quantidade] = 1.46; // average collected
		nucleo_proximity_free[atomos_quantidade] = 3.1;  // Van der waals angstron
		electron_quantidade[atomos_quantidade] = 4;
		electron_raio[atomos_quantidade] = 2.92;
		electron_raio_HB[atomos_quantidade] = 3.1;
		massa[atomos_quantidade] = 14.0067;
		nucleo_proximity_HB[atomos_quantidade] = 1.55;
	} else if (atomo_ll == "S") {
		atomo_letra[atomos_quantidade] = 'S';
		nucleo_proximity[atomos_quantidade] = 1.8;
		nucleo_proximity_free[atomos_quantidade] = 3.6;  // Van der waals angstron
		electron_quantidade[atomos_quantidade] = 4;
		electron_raio[atomos_quantidade] = 3.6;
		massa[atomos_quantidade] = 32.06;
	}
	posx[atomos_quantidade] = 1.0;
	posy[atomos_quantidade] = 3.0;
	posz[atomos_quantidade] = 2.0;
	velocidade_x[atomos_quantidade] = 0.0;
	velocidade_y[atomos_quantidade] = 0.0;
	velocidade_z[atomos_quantidade] = 0.0;
	electron_y[atomos_quantidade][0] = 345.0;
	electron_z[atomos_quantidade][0] = 350.0;

	electron_y[atomos_quantidade][1] = 100.0;
	electron_z[atomos_quantidade][1] = 150.0;

	electron_y[atomos_quantidade][2] = 33.0;
	electron_z[atomos_quantidade][2] = 200.0;

	electron_y[atomos_quantidade][3] = 300.0;
	electron_z[atomos_quantidade][3] = 1.0;

	atomos_quantidade++;
}

void conecta_eletrons(string procura_amino) {
	connect_electron("N1", "C1");
	connect_electron("C1", "C2");
	connect_electron("C2", "O1");

	if (procura_amino == "ILE") {
		connect_electron("C1", "C3");
		connect_electron("C3", "C4");
		connect_electron("C3", "C5");
		connect_electron("C5", "C6");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C4");
		connect_electron("H5", "C4");
		connect_electron("H6", "C4");
		connect_electron("H7", "C5");
		connect_electron("H8", "C5");
		connect_electron("H9", "C6");
		connect_electron("H10", "C6");
		connect_electron("H11", "C6");
	} else if (procura_amino == "GLU") {
		connect_electron("C1", "C3");
		connect_electron("C3", "C4");
		connect_electron("C4", "C5");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "C4");
		connect_electron("H6", "C4");
	} else if (procura_amino == "ALA") {
		connect_electron("C1", "C3");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "C3");
	} else if (procura_amino == "ARG") {
		connect_electron("C1", "C3");
		connect_electron("C3", "C4");
		connect_electron("C4", "C5");
		connect_electron("N2", "C5");
		connect_electron("C6", "N2");
		connect_electron("N3", "C6");
		connect_electron("N4", "C6");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "C4");
		connect_electron("H6", "C4");
		connect_electron("H7", "C5");
		connect_electron("H8", "C5");
		connect_electron("H9", "N2");
		connect_electron("H10", "N3");
		connect_electron("H11", "N3");
		connect_electron("H12", "N4");
		connect_electron("H13", "N4");
	} else if (procura_amino == "LYS") {
		connect_electron("C1", "C3");
		connect_electron("C3", "C4");
		connect_electron("C4", "C5");
		connect_electron("C5", "C6");
		connect_electron("N2", "C6");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "C4");
		connect_electron("H6", "C4");
		connect_electron("H7", "C5");
		connect_electron("H8", "C5");
		connect_electron("H9", "C6");
		connect_electron("H10", "C6");
		connect_electron("H11", "N2");
		connect_electron("H12", "N2");
		connect_electron("H13", "N2");
	} else if (procura_amino == "ASP") {
		connect_electron("C1", "C3");
		connect_electron("C4", "C3");
		connect_electron("O2", "C4");
		connect_electron("O3", "C4");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
	} else if (procura_amino == "TYR") {
		connect_electron("C1", "C3");
		connect_electron("C4", "C3");
		connect_electron("C5", "C4");
		connect_electron("C6", "C5");
		connect_electron("C7", "C6");
		connect_electron("C8", "C7");
		connect_electron("C9", "C8");
		connect_electron("C9", "C4");
		connect_electron("O2", "C7");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "C9");
		connect_electron("H6", "C8");
		connect_electron("H7", "O2");
		connect_electron("H8", "C6");
		connect_electron("H9", "C5");
	} else if (procura_amino == "TRP") {
		connect_electron("C1", "C3");
		connect_electron("C4", "C3");
		connect_electron("C5", "C4");
		connect_electron("C6", "C4");
		connect_electron("C7", "C6");
		connect_electron("C8", "C7");
		connect_electron("C9", "C8");
		connect_electron("C10", "C9");
		connect_electron("C11", "C10");
		connect_electron("C11", "C6");
		connect_electron("N2", "C5");
		connect_electron("N2", "C11");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "C5");
		connect_electron("H6", "N2");
		connect_electron("H7", "C7");
		connect_electron("H8", "C8");
		connect_electron("H9", "C9");
		connect_electron("H10", "C10");
	} else if (procura_amino == "SER") {
		connect_electron("C1", "C3");
		connect_electron("O2", "C3");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
	} else if (procura_amino == "THR") {
		connect_electron("C1", "C3");
		connect_electron("C4", "C3");
		connect_electron("O2", "C3");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C4");
		connect_electron("H5", "C4");
		connect_electron("H6", "C4");
	} else if (procura_amino == "GLY") {
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C1");
	} else if (procura_amino == "MET") {
		connect_electron("C3", "C1");
		connect_electron("C4", "C3");
		connect_electron("S1", "C4");
		connect_electron("C5", "S1");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "C4");
		connect_electron("H6", "C4");
		connect_electron("H7", "C5");
		connect_electron("H8", "C5");
		connect_electron("H9", "C5");
	} else if (procura_amino == "PHE") {
		connect_electron("C1", "C3");
		connect_electron("C4", "C3");
		connect_electron("C5", "C4");
		connect_electron("C6", "C5");
		connect_electron("C7", "C6");
		connect_electron("C8", "C7");
		connect_electron("C9", "C8");
		connect_electron("C9", "C4");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "C5");
		connect_electron("H6", "C6");
		connect_electron("H7", "C7");
		connect_electron("H8", "C8");
		connect_electron("H9", "C9");
	} else if (procura_amino == "LEU") {
		connect_electron("C1", "C3");
		connect_electron("C3", "C4");
		connect_electron("C4", "C5");
		connect_electron("C4", "C6");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "C4");
		connect_electron("H6", "C5");
		connect_electron("H7", "C5");
		connect_electron("H8", "C5");
		connect_electron("H9", "C6");
		connect_electron("H10", "C6");
		connect_electron("H11", "C6");
	} else if (procura_amino == "VAL") {
		connect_electron("C1", "C3");
		connect_electron("C3", "C4");
		connect_electron("C3", "C5");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C4");
		connect_electron("H5", "C4");
		connect_electron("H6", "C4");
		connect_electron("H7", "C5");
		connect_electron("H8", "C5");
		connect_electron("H9", "C5");
	} else if (procura_amino == "ASN") {
		connect_electron("C1", "C3");
		connect_electron("C4", "C3");
		connect_electron("N2", "C4");
		connect_electron("O2", "C4");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "N2");
		connect_electron("H6", "N2");
	} else if (procura_amino == "GLN") {
		connect_electron("C3", "C1");
		connect_electron("C4", "C3");
		connect_electron("C5", "C4");
		connect_electron("N2", "C5");
		connect_electron("O2", "C5");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "C4");
		connect_electron("H6", "C4");
		connect_electron("H7", "N2");
		connect_electron("H8", "N2");
	} else if (procura_amino == "CYS") {
		connect_electron("C3", "C1");
		connect_electron("S1", "C3");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
	} else if (procura_amino == "HIS") {
		connect_electron("C3", "C1");
		connect_electron("C4", "C3");
		connect_electron("C5", "C4");
		connect_electron("N2", "C5");
		connect_electron("C6", "N2");
		connect_electron("N3", "C6");
		connect_electron("N3", "C4");
		connect_electron("H1", "N1");
		connect_electron("H2", "C1");
		connect_electron("H3", "C3");
		connect_electron("H4", "C3");
		connect_electron("H5", "C5");
		connect_electron("H6", "C6");
		connect_electron("H7", "N3");
	} else if (procura_amino == "PRO") {
		connect_electron("C3", "C1");
		connect_electron("C4", "C3");
		connect_electron("C5", "C4");
		connect_electron("N1", "C5");
		connect_electron("H1", "C1");
		connect_electron("H2", "C3");
		connect_electron("H3", "C3");
		connect_electron("H4", "C4");
		connect_electron("H5", "C4");
		connect_electron("H6", "C5");
		connect_electron("H7", "C5");
	}
}

void read_pdb_amino(string procura_amino) {
	cont_amin = 0;
	printf("Lendo arquivo PDB dos aminoacidos\n");
	string line = "";
	ifstream myfile("files/aminoacidos.pdb");
	if (myfile.is_open()) {
		string tipo;
		string nome_amino;
		string atomo_sigla;
		string atomo_label;
		GLfloat fix_x = 0;
		GLfloat fix_y = 0;
		GLfloat fix_z = 0;
		while (getline(myfile, line)) {
			tipo = line.substr(0, 6);
			if (tipo.erase(tipo.find_last_not_of(" ") + 1) == "ATOM" || tipo == "HETATM") {
				nome_amino = line.substr(17, 3);

				if (nome_amino == procura_amino) {
					cont_amin++;
					atomo_label = line.substr(77, 3);
					atomo_label = atomo_label.erase(atomo_label.find_last_not_of(" ") + 1);

					atomo_sigla = line.substr(77, 1);
					pdbamin_coordenadas_x.resize(cont_amin);
					pdbamin_coordenadas_y.resize(cont_amin);
					pdbamin_coordenadas_z.resize(cont_amin);

					pdbamin_coordenadas_x[cont_amin - 1] = strtof(line.substr(30, 8).c_str(), NULL) - fix_x;
					pdbamin_coordenadas_y[cont_amin - 1] = strtof(line.substr(38, 8).c_str(), NULL) - fix_y;
					pdbamin_coordenadas_z[cont_amin - 1] = strtof(line.substr(46, 8).c_str(), NULL) - fix_z;

					if (cont_amin == 1) {
						printf("Zerando coordenadas\n");
						fix_x = pdbamin_coordenadas_x[cont_amin - 1];
						fix_y = pdbamin_coordenadas_y[cont_amin - 1];
						fix_z = pdbamin_coordenadas_z[cont_amin - 1];
						pdbamin_coordenadas_x[cont_amin - 1] = 0.0;
						pdbamin_coordenadas_y[cont_amin - 1] = 0.0;
						pdbamin_coordenadas_z[cont_amin - 1] = 0.0;
					}
//					printf("%s %s %s %f %f %f\n", atomo_label.c_str(), atomo_sigla.c_str(), nome_amino.c_str(), pdbamin_coordenadas_x[cont_amin - 1], pdbamin_coordenadas_y[cont_amin - 1], pdbamin_coordenadas_z[cont_amin - 1]);
					add_atoml(atomo_sigla);
					mapa_atomo[atomo_label] = atomos_quantidade - 1;
					posx[mapa_atomo[atomo_label]] = pdbamin_coordenadas_x[cont_amin - 1];
					posy[mapa_atomo[atomo_label]] = pdbamin_coordenadas_y[cont_amin - 1];
					posz[mapa_atomo[atomo_label]] = pdbamin_coordenadas_z[cont_amin - 1];
					amino[mapa_atomo[atomo_label]] = contador_amino;
					amino[mapa_atomo[atomo_label]] = get_amino_number(nome_amino.c_str());
					amino_sequencial[mapa_atomo[atomo_label]] = contador_amino;
					atomo_letraN[mapa_atomo[atomo_label]] = get_atom_number(atomo_label.c_str());
					if (atomo_label == "N1" || atomo_label == "C1" || atomo_label == "C2") {
						atomo_base[mapa_atomo[atomo_label]] = true;
					}
					if (atomo_label == "N1" && contador_amino > 0) {
						electron_arested[mapa_atomo[atomo_label]][0][0] = chain_ultimo_atomo;
						electron_arested[mapa_atomo[atomo_label]][0][1] = 1;
						electron_arested[chain_ultimo_atomo][1][0] = mapa_atomo[atomo_label];
						electron_arested[chain_ultimo_atomo][1][1] = 0;
						electron_contador[atomo_label]++;
					} else if (atomo_label == "C2") {
						chain_ultimo_atomo = mapa_atomo[atomo_label];
					}

//					printf("---> %d %d %s %d\n", atomo_letraN[mapa_atomo[atomo_label]], mapa_atomo[atomo_label], atomo_label.c_str(), get_amino_number(nome_amino.c_str()));

				}

//				printf("%s\n",nome_amino.c_str());

			}
		}
//		velocidade_z[mapa_atomo["C1"]] = 0.14;
//		velocidade_z[mapa_atomo["N1"]] = 0.13;
//		velocidade_z[mapa_atomo["C1"]] = 0.04;
//		velocidade_z[mapa_atomo["N1"]] = 0.03;
		velocidade_z[mapa_atomo["C1"]] = 0.001;
		velocidade_z[mapa_atomo["N1"]] = 0.003;
		conecta_eletrons(procura_amino);

		myfile.close();
		mapa_atomo.clear();
		electron_contador.clear();
		contador_amino++;
	} else {
		cout << "Unable to open file\n";
	}
}

void add_pdb_amino() {
	for (int i = 0; i < cont_amin; i++) {
		glColor3ub(255, 255, 255);
		glTranslatef(pdbamin_coordenadas_x[i], pdbamin_coordenadas_y[i], pdbamin_coordenadas_z[i]);
		glutSolidCube(0.1f);
//		glutSolidSphere(0.5f, 30, 30);
		glTranslatef(-pdbamin_coordenadas_x[i], -pdbamin_coordenadas_y[i], -pdbamin_coordenadas_z[i]);
	}
}

void load_protein(string PDBID) {
	printf("Lendo proteina a partir do banco: %s\n", PDBID.c_str());
	printf("Calibrando a partir do mysql\n");
	try {
		sql::Driver *driver;
		sql::Connection *con;
		sql::PreparedStatement *pstmt;
		sql::ResultSet *res;
		driver = get_driver_instance();
		con = driver->connect("tcp://127.0.0.1:3306", "a00s_230", "testando");
		con->setSchema("a00s_230");

//		pstmt = con->prepareStatement("SELECT i_306337 amino, i_306344 aminoseq, i_306299 x, i_306307 y, i_306315 z, i_306408 atom, i_331770 atomlabel FROM a_306280 WHERE i_307676=? AND i_306401=1 AND i_331770 IS NOT NULL HAVING aminoseq IN(1) AND atom IN ('C') ORDER BY i_306344,atomlabel DESC");
//		pstmt = con->prepareStatement("SELECT i_306337 amino, i_306344 aminoseq, i_306299 x, i_306307 y, i_306315 z, i_306408 atom, i_331770 atomlabel FROM a_306280 WHERE i_307676=? AND i_306401=1 AND i_331770 IS NOT NULL ORDER BY i_306344,atomlabel DESC");
		pstmt = con->prepareStatement("SELECT i_306337 amino, i_306344 aminoseq, i_306299 x, i_306307 y, i_306315 z, i_306408 atom, i_331770 atomlabel FROM a_306280 WHERE i_307676=? AND i_306401=1 AND i_331770 IS NOT NULL AND i_306344=1 ORDER BY i_306344,atomlabel DESC");
		pstmt->setString(1, PDBID);
		res = pstmt->executeQuery();
		string nome_amino;
		string atomo_sigla;
		string atomo_label;
		int aminoseq = -1;
//		contador_amino = -1;
		while (res->next()) {
//			printf("%s %s\n", res->getString("amino").c_str(), res->getString("atomlabel").c_str());
			if (res->getInt("aminoseq") != aminoseq) {
				if (aminoseq != -1) {
					contador_amino++;
					conecta_eletrons(nome_amino);
				}
				aminoseq = res->getInt("aminoseq");
				mapa_atomo.clear();
				electron_contador.clear();
			}
			nome_amino = res->getString("amino");
			aminoseq = res->getInt("aminoseq");
			cont_amin++;
			atomo_label = res->getString("atomlabel");
			atomo_sigla = res->getString("atom");
			pdbamin_coordenadas_x.resize(cont_amin);
			pdbamin_coordenadas_y.resize(cont_amin);
			pdbamin_coordenadas_z.resize(cont_amin);
			pdbamin_coordenadas_x[cont_amin - 1] = strtof(res->getString("x").c_str(), NULL);
			pdbamin_coordenadas_y[cont_amin - 1] = strtof(res->getString("y").c_str(), NULL);
			pdbamin_coordenadas_z[cont_amin - 1] = strtof(res->getString("z").c_str(), NULL);
			add_atoml(atomo_sigla);
			mapa_atomo[atomo_label] = atomos_quantidade - 1;
			posx[mapa_atomo[atomo_label]] = pdbamin_coordenadas_x[cont_amin - 1];
			posy[mapa_atomo[atomo_label]] = pdbamin_coordenadas_y[cont_amin - 1];
			posz[mapa_atomo[atomo_label]] = pdbamin_coordenadas_z[cont_amin - 1];
			amino[mapa_atomo[atomo_label]] = contador_amino;
			amino[mapa_atomo[atomo_label]] = get_amino_number(nome_amino.c_str());
			amino_sequencial[mapa_atomo[atomo_label]] = contador_amino;
			atomo_letraN[mapa_atomo[atomo_label]] = get_atom_number(atomo_label.c_str());
			if (atomo_label == "N1" || atomo_label == "C1" || atomo_label == "C2") {
				atomo_base[mapa_atomo[atomo_label]] = true;
			}
			if (atomo_label == "N1" && contador_amino > 0) {
//				printf("Interligando com o ultimo atomo %d\n",chain_ultimo_atomo);
				electron_arested[mapa_atomo[atomo_label]][0][0] = chain_ultimo_atomo;
				electron_arested[mapa_atomo[atomo_label]][0][1] = 1;
				electron_arested[chain_ultimo_atomo][1][0] = mapa_atomo[atomo_label];
				electron_arested[chain_ultimo_atomo][1][1] = 0;
				electron_contador[atomo_label]++;
			} else if (atomo_label == "C2") {
//				printf("Setando ultimo atomo %d\n",mapa_atomo[atomo_label]);
				chain_ultimo_atomo = mapa_atomo[atomo_label];
			}
		}
//		velocidade_z[mapa_atomo["C1"]] = 0.001;
//		velocidade_z[mapa_atomo["N1"]] = 0.003;
		conecta_eletrons(nome_amino);
		delete res;
		delete pstmt;
		delete con;
		printf("PDBid Loaded\n");
	} catch (sql::SQLException &e) {
		printf("%d\n", e.getErrorCode());
	}
}

GLfloat compare_protein_calculate() {
//	printf("Comparando\n");
	uint tamanho_vetor_1 = AminoTemp1.size();
	uint tamanho_vetor_2 = AminoTemp2.size();
	GLfloat distance_local = 0.0;
	for (uint i = 0; i < tamanho_vetor_1; ++i) {
		for (uint ii = 0; ii < tamanho_vetor_2; ++ii) {
//			printf("%f %f %f / %f %f %f\n",AminoTemp1[i].x,AminoTemp1[i].y,AminoTemp1[i].z,AminoTemp2[ii].x,AminoTemp2[ii].y,AminoTemp2[ii].z);
//			printf("D: (%f %f %f | %f %f %f)   %f        %f\n", AminoTemp2[i].x, AminoTemp2[i].y, AminoTemp2[i].z, AminoTemp2[ii].x, AminoTemp2[ii].y, AminoTemp2[ii].z, abs(sqrt(pow((AminoTemp2[ii].x - AminoTemp1[i].x), 2.0) + pow((AminoTemp2[ii].y - AminoTemp1[i].y), 2.0) + pow((AminoTemp2[ii].z - AminoTemp1[i].z), 2.0))), distance_local);
			distance_local += abs(sqrt(pow((AminoTemp2[ii].x - AminoTemp1[i].x), 2.0) + pow((AminoTemp2[ii].y - AminoTemp1[i].y), 2.0) + pow((AminoTemp2[ii].z - AminoTemp1[i].z), 2.0)));
		}
	}
//	printf("%f\n",distance_local);
	return distance_local;
}

GLfloat compare_protein_calculateMD() {
	uint tamanho_vetor_1 = AminoTemp1MD.size();
	uint tamanho_vetor_2 = AminoTemp2MD.size();
	GLfloat distance_local = 0.0;
	for (uint i = 0; i < tamanho_vetor_1; ++i) {
		for (uint ii = 0; ii < tamanho_vetor_2; ++ii) {
//			printf("%f %f %f / %f %f %f\n",AminoTemp1[i].x,AminoTemp1[i].y,AminoTemp1[i].z,AminoTemp2[ii].x,AminoTemp2[ii].y,AminoTemp2[ii].z);
			distance_local += abs(sqrt(pow((AminoTemp2MD[ii].x - AminoTemp1MD[i].x), 2.0) + pow((AminoTemp2MD[ii].y - AminoTemp1MD[i].y), 2.0) + pow((AminoTemp2MD[ii].z - AminoTemp1MD[i].z), 2.0)));
		}
	}
//	printf("%f\n",distance_local);
	return distance_local;
}

void compare_protein_build_MD(int c_distancia) {
	uint tamanho_vetor = pdb_analise.size();
	int last_seq = -1;
	int cont_seq = 0;
	int cont_vector = -1;
	int i1c = 0;
	int i2c = c_distancia;
	GLfloat soma_total = 0.0;

	for (GLint i = 0; i < atomos_quantidade; i++) {
//		printf("3d\n");
		if (last_seq >= 0) {
			if (last_seq != amino_sequencial[i]) {
				cont_seq++;
				cont_vector = -1;
				if (cont_seq > i2c) {
					soma_total += compare_protein_calculateMD();
					cont_seq = 1; // antes tava 0
					i1c++;
					i2c++;
					i = 0;
				}
			}
		} else {
			last_seq = amino_sequencial[i];
		}
//		printf("Last seq: %d\n", last_seq);
		if (i1c == cont_seq) {
			cont_vector++;
			AminoTemp1MD.resize(cont_vector + 1);
			AminoTemp1MD[cont_vector].x = posx[i];
			AminoTemp1MD[cont_vector].y = posy[i];
			AminoTemp1MD[cont_vector].z = posz[i];
		}
		if (i2c == cont_seq) {
			if (i2c != i1c) {
				cont_vector++;
			}
			AminoTemp2MD.resize(cont_vector + 1);
			AminoTemp2MD[cont_vector].x = posx[i];
			AminoTemp2MD[cont_vector].y = posy[i];
			AminoTemp2MD[cont_vector].z = posz[i];
		}
		last_seq = amino_sequencial[i];
	}
	soma_total += compare_protein_calculate();
	printf("Total final1: %f\n", soma_total);
}

void compare_protein_build(int c_distancia) {
//	c_distancia++;
	uint tamanho_vetor = pdb_analise.size();
	printf("%d\n", tamanho_vetor);
	int last_seq = -1;
	int cont_seq = 0;
	int cont_vector = -1;
	int i1c = 0;
	int i2c = c_distancia;
	GLfloat soma_total = 0.0;
	for (uint i = 0; i < tamanho_vetor; ++i) {
		if (last_seq > 0) {
			if (last_seq != pdb_analise[i].aminoseq) {
				cont_seq++;
				cont_vector = -1;
				if (cont_seq > i2c) {
					soma_total += compare_protein_calculate();
					cont_seq = 1;
					i1c++;
					i2c++;
					i = 0;
				}
			}
		} else {
			last_seq = pdb_analise[i].aminoseq;
		}
//		printf("ContsetA %d\n",cont_seq);
		if (i1c == cont_seq) {
//			printf("Aqui 1: %d\n", cont_vector);
			cont_vector++;
			AminoTemp1.resize(cont_vector + 1);
			AminoTemp1[cont_vector].x = pdb_analise[i].x;
			AminoTemp1[cont_vector].y = pdb_analise[i].y;
			AminoTemp1[cont_vector].z = pdb_analise[i].z;
		}
//		printf("ContsetB %d\n",cont_seq);
		if (i2c == cont_seq) {
//			printf("Aqui 2\n");
			if (i2c != i1c) {
				cont_vector++;
			}
			AminoTemp2.resize(cont_vector + 1);
			AminoTemp2[cont_vector].x = pdb_analise[i].x;
			AminoTemp2[cont_vector].y = pdb_analise[i].y;
			AminoTemp2[cont_vector].z = pdb_analise[i].z;
		}
//		printf("ContsetC %d\n",cont_seq);
		last_seq = pdb_analise[i].aminoseq;
	}
	soma_total += compare_protein_calculate();
	printf("Total final2: %f\n", soma_total);
}

void load_protein_position(string PDBID) {
	printf("Lendo proteina para analise: %s\n", PDBID.c_str());
	try {
		sql::Driver *driver;
		sql::Connection *con;
		sql::PreparedStatement *pstmt;
		sql::ResultSet *res;
		driver = get_driver_instance();
		con = driver->connect("tcp://127.0.0.1:3306", "a00s_230", "testando");
		con->setSchema("a00s_230");
//		pstmt = con->prepareStatement("SELECT i_306344 aminoseq, i_306299 x, i_306307 y, i_306315 z, i_331770 atomlabel FROM a_306280 WHERE i_307676=? AND i_306401=1 AND i_331770 IS NOT NULL ORDER BY i_306344,atomlabel DESC");
		pstmt = con->prepareStatement("SELECT i_306344 aminoseq, i_306299 x, i_306307 y, i_306315 z, i_331770 atomlabel FROM a_306280 WHERE i_307676=? AND i_306401=1 AND i_331770 IS NOT NULL AND i_306344=1 ORDER BY i_306344,atomlabel DESC");
		pstmt->setString(1, PDBID);
		res = pstmt->executeQuery();
		int cont_atom = -1;
		while (res->next()) {
			cont_atom++;
			pdb_analise.resize(cont_atom + 1);
			pdb_analise[cont_atom].aminoseq = res->getInt("aminoseq");
			pdb_analise[cont_atom].label = res->getString("atomlabel").c_str();
			pdb_analise[cont_atom].x = res->getDouble("x");
			pdb_analise[cont_atom].y = res->getDouble("y");
			pdb_analise[cont_atom].z = res->getDouble("z");
		}
		delete res;
		delete pstmt;
		delete con;
		printf("PDBid Loaded\n");
//		compare_protein_build_MD(0);
		compare_protein_build(0);
	} catch (sql::SQLException &e) {
		printf("%d\n", e.getErrorCode());
	}
}

