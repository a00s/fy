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
#include "mysql_connection.h"
#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>
#include <cppconn/prepared_statement.h>
#include "gvariaveis.h"
#include "struct.h"
using namespace std;

extern vector<VectorProtein> ProteinOriginal;
extern GLint calibration_precision;
extern GLint calibration_precision_out;
extern GLfloat calibrationMin[20][34][2][20][34]; // [Amino1][Atom1][Same|Another][Amino2][Atom2]
extern GLfloat calibrationMax[20][34][2][20][34]; // [Amino1][Atom1][Same|Another][Amino2][Atom2]
extern map<string, map<string, map<string, map<string, map<string, map<string, GLfloat> > > > > > vs;
extern int atomos_quantidade;

int get_amino_number(const char *amino_sigla) {
	//  Glu = 0
	//  Tyr = 1
	//  Gln = 2
	//  Ile = 3
	//  Trp = 4
	//  Lys = 5
	//  Leu = 6
	//	Ala = 7
	//	Arg = 8
	//	Asp = 9
	//	Ser = 10
	//	Thr = 11
	//	Gly = 12
	//	Met = 13
	//	Phe = 14
	//	Val = 15
	//	Asn = 16
	//	Cys = 17
	//	His = 18
	//	Pro = 19

	if (!strcmp(amino_sigla, "GLU")) {
		return 0;
	} else if (!strcmp(amino_sigla, "TYR")) {
		return 1;
	} else if (!strcmp(amino_sigla, "GLN")) {
		return 2;
	} else if (!strcmp(amino_sigla, "ILE")) {
		return 3;
	} else if (!strcmp(amino_sigla, "TRP")) {
		return 4;
	} else if (!strcmp(amino_sigla, "LYS")) {
		return 5;
	} else if (!strcmp(amino_sigla, "LEU")) {
		return 6;
	} else if (!strcmp(amino_sigla, "ALA")) {
		return 7;
	} else if (!strcmp(amino_sigla, "ARG")) {
		return 8;
	} else if (!strcmp(amino_sigla, "ASP")) {
		return 9;
	} else if (!strcmp(amino_sigla, "SER")) {
		return 10;
	} else if (!strcmp(amino_sigla, "THR")) {
		return 11;
	} else if (!strcmp(amino_sigla, "GLY")) {
		return 12;
	} else if (!strcmp(amino_sigla, "MET")) {
		return 13;
	} else if (!strcmp(amino_sigla, "PHE")) {
		return 14;
	} else if (!strcmp(amino_sigla, "VAL")) {
		return 15;
	} else if (!strcmp(amino_sigla, "ASN")) {
		return 16;
	} else if (!strcmp(amino_sigla, "CYS")) {
		return 17;
	} else if (!strcmp(amino_sigla, "HIS")) {
		return 18;
	} else if (!strcmp(amino_sigla, "PRO")) {
		return 19;
	}
	return 99;
}

int get_atom_number(const char *atomo_letra_local) {
	if (!strcmp(atomo_letra_local, "C1")) {
		return 0;
	} else if (!strcmp(atomo_letra_local, "C2")) {
		return 1;
	} else if (!strcmp(atomo_letra_local, "C3")) {
		return 2;
	} else if (!strcmp(atomo_letra_local, "C4")) {
		return 3;
	} else if (!strcmp(atomo_letra_local, "C5")) {
		return 4;
	} else if (!strcmp(atomo_letra_local, "C6")) {
		return 5;
	} else if (!strcmp(atomo_letra_local, "C7")) {
		return 6;
	} else if (!strcmp(atomo_letra_local, "C8")) {
		return 7;
	} else if (!strcmp(atomo_letra_local, "C9")) {
		return 8;
	} else if (!strcmp(atomo_letra_local, "C10")) {
		return 9;
	} else if (!strcmp(atomo_letra_local, "C11")) {
		return 10;
	} else if (!strcmp(atomo_letra_local, "N1")) {
		return 11;
	} else if (!strcmp(atomo_letra_local, "N2")) {
		return 12;
	} else if (!strcmp(atomo_letra_local, "N3")) {
		return 13;
	} else if (!strcmp(atomo_letra_local, "N4")) {
		return 14;
	} else if (!strcmp(atomo_letra_local, "O1")) {
		return 15;
	} else if (!strcmp(atomo_letra_local, "O2")) {
		return 16;
	} else if (!strcmp(atomo_letra_local, "O3")) {
		return 17;
	} else if (!strcmp(atomo_letra_local, "O4")) {
		return 18;
	} else if (!strcmp(atomo_letra_local, "H1")) {
		return 19;
	} else if (!strcmp(atomo_letra_local, "H2")) {
		return 20;
	} else if (!strcmp(atomo_letra_local, "H3")) {
		return 21;
	} else if (!strcmp(atomo_letra_local, "H4")) {
		return 22;
	} else if (!strcmp(atomo_letra_local, "H5")) {
		return 23;
	} else if (!strcmp(atomo_letra_local, "H6")) {
		return 24;
	} else if (!strcmp(atomo_letra_local, "H7")) {
		return 25;
	} else if (!strcmp(atomo_letra_local, "H8")) {
		return 26;
	} else if (!strcmp(atomo_letra_local, "H9")) {
		return 27;
	} else if (!strcmp(atomo_letra_local, "H10")) {
		return 28;
	} else if (!strcmp(atomo_letra_local, "H11")) {
		return 29;
	} else if (!strcmp(atomo_letra_local, "H12")) {
		return 30;
	} else if (!strcmp(atomo_letra_local, "H13")) {
		return 31;
	} else if (!strcmp(atomo_letra_local, "H14")) {
		return 32;
	} else if (!strcmp(atomo_letra_local, "S1")) {
		return 33;
	}
	return 0;
}

void distance_calibration() {
	printf("Calibrando a partir do mysql\n");
	try {
		sql::Driver *driver;
		sql::Connection *con;
		sql::PreparedStatement *pstmt;
		sql::ResultSet *res;
		driver = get_driver_instance();
		con = driver->connect("tcp://127.0.0.1:3306", "a00s_230", "testando");
		con->setSchema("a00s_230");
		printf("Calibration precision %d\n", calibration_precision);
//		printf("Aqui 0.3\n");
		// ------ Dictionary -------
		// Same = 0 / Different = 1
		// Glu = 0
		// Tyr = 1
		// Gln = 2
		// Ile = 3
		// Trp = 4
		// Lys = 5

		// C1 = 0
		// C2 = 1
		// C3 = 2
		// C4 = 3
		// C5 = 4
		// C6 = 5
		// C7 = 6
		// C8 = 7
		// C9 = 8
		// C10 = 9
		// C11 = 10

		// N1 = 11
		// N2 = 12
		// N3 = 13
		// N4 = 14

		// O1 = 15
		// O2 = 16
		// O3 = 17
		// O4 = 18

		// H1 = 19
		// H2 = 20
		// H3 = 21
		// H4 = 22
		// H5 = 23
		// H6 = 24
		// H7 = 25
		// H8 = 26
		// H9 = 27
		// H10 = 28
		// H11 = 29
		// H12 = 30
		// H13 = 31
		// H14 = 32

		// S1 = 33

		// [Amino1][Atom1][Same|Different][Amino2][Atom2]
		// -------- Same amino acid ---------
//		SELECT * FROM (SELECT * FROM a_380484 WHERE i_380512 = 1 ORDER BY i_380529 DESC) tabtemp GROUP BY i_380488,i_380500,i_380494,i_380506
//		pstmt = con->prepareStatement("SELECT i_380488 amino, i_380494 atom1, i_380506 atom2, i_380517 min_distance, i_380523 max_distance FROM a_380484 WHERE i_380529=? AND i_380512=1");
		pstmt = con->prepareStatement("SELECT i_380488 amino, i_380494 atom1, i_380506 atom2, i_380517 min_distance, i_380523 max_distance FROM (SELECT * FROM a_380484 WHERE i_380512 = 1 AND i_393242/*pdbid*/ IS NULL ORDER BY i_380529 DESC) tabtemp GROUP BY i_380488,i_380500,i_380494,i_380506");
//		pstmt->setInt(1, calibration_precision);
		res = pstmt->executeQuery();
		while (res->next()) {
			calibrationMin[get_amino_number(res->getString("amino").c_str())][get_atom_number(res->getString("atom1").c_str())][0][get_amino_number(res->getString("amino").c_str())][get_atom_number(res->getString("atom2").c_str())] = res->getDouble(4);
			calibrationMax[get_amino_number(res->getString("amino").c_str())][get_atom_number(res->getString("atom1").c_str())][0][get_amino_number(res->getString("amino").c_str())][get_atom_number(res->getString("atom2").c_str())] = res->getDouble(5);
			//vs["calibrationMin"][res->getString("amino").c_str()][res->getString("atom1").c_str()]["Same"][res->getString("amino").c_str()][res->getString("atom2").c_str()] = res->getDouble(4);
			//vs["calibrationMax"][res->getString("amino").c_str()][res->getString("atom1").c_str()]["Same"][res->getString("amino").c_str()][res->getString("atom2").c_str()] = res->getDouble(5);
		}

		// -------- Different amino acid ---------
//		pstmt = con->prepareStatement("SELECT i_380488 aminofrom, i_380500 aminoto, i_380494 atom1, i_380506 atom2, i_380517 min_distance, i_380523 max_distance FROM a_380484 WHERE i_380529=? AND i_380512 IS NULL");
		pstmt = con->prepareStatement("SELECT i_380488 aminofrom, i_380500 aminoto, i_380494 atom1, i_380506 atom2, i_380517 min_distance, i_380523 max_distance FROM (SELECT * FROM a_380484 WHERE i_380512 IS NULL AND i_393242/*pdbid*/ IS NULL ORDER BY i_380529 DESC) tabtemp GROUP BY i_380488,i_380500,i_380494,i_380506");
//		pstmt->setInt(1, calibration_precision_out);
		res = pstmt->executeQuery();
		while (res->next()) {
			calibrationMin[get_amino_number(res->getString("aminofrom").c_str())][get_atom_number(res->getString("atom1").c_str())][1][get_amino_number(res->getString("aminoto").c_str())][get_atom_number(res->getString("atom2").c_str())] = res->getDouble(5);
		}
		delete res;
		delete pstmt;
		delete con;
		printf("System calibrated\n");
	} catch (sql::SQLException &e) {
		printf("Error %d\n", e.getErrorCode());
	}
}

void distance_calibration_pdb(string PDBID) {
	printf("Calibrando a partir de um pdb no mysql\n");
	try {
		sql::Driver *driver;
		sql::Connection *con;
		sql::PreparedStatement *pstmt;
		sql::ResultSet *res;
		driver = get_driver_instance();
		con = driver->connect("tcp://127.0.0.1:3306", "a00s_230", "testando");
		con->setSchema("a00s_230");
		// [Amino1][Atom1][Same|Different][Amino2][Atom2]
		// -------- Same amino acid ---------
		pstmt = con->prepareStatement("SELECT i_380488 amino, i_380494 atom1, i_380506 atom2, i_380517 min_distance, i_380523 max_distance FROM a_380484 WHERE i_380512 = 1 AND i_393242/*pdbid*/ = ? ORDER BY i_380529 DESC");
		pstmt->setString(1, PDBID);
		res = pstmt->executeQuery();
		while (res->next()) {
//			printf("Recalibrando %s %s %s\n",res->getString("amino").c_str(),res->getString("atom1").c_str(),res->getString("atom2").c_str());
			calibrationMin[get_amino_number(res->getString("amino").c_str())][get_atom_number(res->getString("atom1").c_str())][0][get_amino_number(res->getString("amino").c_str())][get_atom_number(res->getString("atom2").c_str())] = res->getDouble(4);
			calibrationMax[get_amino_number(res->getString("amino").c_str())][get_atom_number(res->getString("atom1").c_str())][0][get_amino_number(res->getString("amino").c_str())][get_atom_number(res->getString("atom2").c_str())] = res->getDouble(5);
			vs["PDBcalibrationMin"][res->getString("amino").c_str()][res->getString("atom1").c_str()]["Same"][res->getString("amino").c_str()][res->getString("atom2").c_str()] = res->getDouble(4);
			vs["PDBcalibrationMax"][res->getString("amino").c_str()][res->getString("atom1").c_str()]["Same"][res->getString("amino").c_str()][res->getString("atom2").c_str()] = res->getDouble(5);
		}
		// -------- Different amino acid ---------
		pstmt = con->prepareStatement("SELECT i_380488 amino, i_380494 atom1, i_380500 amino2, i_380506 atom2, i_380517 min_distance FROM a_380484 WHERE i_380512 IS NULL AND i_393242/*pdbid*/ = ? ORDER BY i_380529 DESC");
		pstmt->setString(1, PDBID);
		res = pstmt->executeQuery();
		while (res->next()) {
//			if (get_amino_number(res->getString("amino").c_str()) == 19 && get_amino_number(res->getString("amino2").c_str()) == 2 && get_atom_number(res->getString("atom1").c_str()) == 0 && get_atom_number(res->getString("atom2").c_str()) == 0) {
//				printf("Recalibrando %s %s %s %s\n", res->getString("amino").c_str(), res->getString("amino2").c_str(), res->getString("atom1").c_str(), res->getString("atom2").c_str());
//			}
			calibrationMin[get_amino_number(res->getString("amino").c_str())][get_atom_number(res->getString("atom1").c_str())][1][get_amino_number(res->getString("amino2").c_str())][get_atom_number(res->getString("atom2").c_str())] = res->getDouble(5);
			vs["PDBcalibrationMin"][res->getString("amino").c_str()][res->getString("atom1").c_str()]["Different"][res->getString("amino2").c_str()][res->getString("atom2").c_str()] = res->getDouble(5);
//			if (get_amino_number(res->getString("amino").c_str()) == 19 && get_amino_number(res->getString("amino2").c_str()) == 2 && get_atom_number(res->getString("atom1").c_str()) == 0 && get_atom_number(res->getString("atom2").c_str()) == 0) {
//				printf("Resultado %f %f %s\n", calibrationMin[19][0][1][2][0], res->getDouble(5),res->getString(5).c_str());
//			}
		}
//		printf("Procurando %f %d %d %d\n", calibrationMin[19][0][1][2][0], get_amino_number("PRO"), get_amino_number("GLN"), get_atom_number("C1"));
		delete res;
		delete pstmt;
		delete con;
		printf("System calibrated PDB\n");
//		typedef map<string, map<string, map<string, map<string, map<string, GLfloat> > > > > inner_map;
//		typedef map<string, inner_map> outer_map;

//		for (outer_map::iterator i = ivs.begin(), iend = ivs.end(); i != iend; ++i)
//		{
//		    inner_map &innerMap = i->second;
//		    for (inner_map::iterator j = innerMap.begin(), jend = innerMap.end(); j != jend; ++j)
//		    {
//		        /* ... */
//		    }
//		}
//		for (map<string, map<string, map<string, map<string, map<string, map<string, GLfloat> > > > > >::iterator i = vs.begin(); i != vs.end(); ++i) {
//			for (map<string, map<string, map<string, map<string, map<string, GLfloat> > > > >::iterator ii = i->second.begin(); ii != i->second.end(); ++ii) {
//				for (map<string, map<string, map<string, map<string, GLfloat> > > >::iterator iii = ii->second.begin(); iii != ii->second.end(); ++iii) {
//					for (map<string, map<string, map<string, GLfloat> > >::iterator iiii = iii->second.begin(); iiii != iii->second.end(); ++iiii) {
//						for (map<string, map<string, GLfloat> >::iterator iiiii = iiii->second.begin(); iiiii != iiii->second.end(); ++iiiii) {
//							for (map<string, GLfloat>::iterator iiiiii = iiiii->second.begin(); iiiiii != iiiii->second.end(); ++iiiiii) {
//								printf("%s %s %s %s %s %s %f\n", (*i).first.c_str(), (*ii).first.c_str(), (*iii).first.c_str(), (*iiii).first.c_str(), (*iiiii).first.c_str(), (*iiiiii).first.c_str(), (*iiiiii).second);
//							}
//						}
//					}
//				}
//			}
//		}
//		map<string, map<string, map<string, map<string, map<string, map<string, GLfloat> > > > > >::iterator ivs;
//		for (ivs = vs.begin(); ivs != vs.end(); ivs++) {
//			printf("%s\n", ivs->first.c_str());
//			map<string, map<string, map<string, map<string, map<string, GLfloat> > > > >::iterator ivs2;
//			//			string tvar = ivs->first.c_str();
//			for (ivs2 = ivs->second; ivs2 != ivs->second.end(); ivs2++) {
//				//				string  = ivs2->first.c_str();
//				//				tw[tvar] = tvar2;
//			}
//		}
	} catch (sql::SQLException &e) {
		printf("%d\n", e.getErrorCode());
	}
}

void hsvtorgb(unsigned char *r, unsigned char *g, unsigned char *b, unsigned char h, unsigned char s, unsigned char v) {
	unsigned char region, fpart, p, q, t;

	if (s == 0) {
		/* color is grayscale */
		*r = *g = *b = v;
		return;
	}

	/* make hue 0-5 */
	region = h / 43;
	/* find remainder part, make it from 0-255 */
	fpart = (h - (region * 43)) * 6;

	/* calculate temp vars, doing integer multiplication */
	p = (v * (255 - s)) >> 8;
	q = (v * (255 - ((s * fpart) >> 8))) >> 8;
	t = (v * (255 - ((s * (255 - fpart)) >> 8))) >> 8;

	/* assign temp vars based on color cone region */
	switch (region) {
	case 0:
		*r = v;
		*g = t;
		*b = p;
		break;
	case 1:
		*r = q;
		*g = v;
		*b = p;
		break;
	case 2:
		*r = p;
		*g = v;
		*b = t;
		break;
	case 3:
		*r = p;
		*g = q;
		*b = v;
		break;
	case 4:
		*r = t;
		*g = p;
		*b = v;
		break;
	default:
		*r = v;
		*g = p;
		*b = q;
		break;
	}

	return;
}

GLfloat CompareWithGhost() {
	GLfloat distance_error = 0;
	for (GLint i = 0; i < atomos_quantidade; i++) {
		if (ProteinOriginal[i].show) {
			distance_error += show_distance(posx[i], posxGhost[i], posy[i], posyGhost[i], posz[i], poszGhost[i]);
		}
	}
//	printf("Distance Error: %f\n",distance_error);
	return distance_error;
}
