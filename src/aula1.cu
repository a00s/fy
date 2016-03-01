#include <iostream>

#include <GL/GLee.h>         // No need to link to GL/gl.h
#include <GL/glfw.h>      // Include OpenGL Framework library
//#include <GLFW/glfw3.h>
#include <GL/freeglut.h>  // Include FreeGLUT so we can easily draw spheres and calculate our viewing frustrum
#include <math.h>         // Used only for sin() and cos() functions
#include <cstdio>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "gvariaveis.h"
#include "struct.h"
#include "matrix3d.h"
//#include <unordered_map>

using namespace std;
// Variaveis de tela
map<string, GLfloat> vs;
// Variaveis propriedades
map<string, map<string, map<string, map<string, map<string, map<string, GLfloat> > > > > > vp;
map<string, map<string, GLfloat*> > vF;
map<string, map<string, GLint*> > vI;
map<string, map<string, bool*> > vB;

const float TO_RADS = 3.141592654f / 180.0f; // The value of 1 degree in radians

GLint windowWidth = 1500;                    // Width of our window
GLint windowHeight = 768;                    // Heightof our window

GLint midWindowX = windowWidth / 2;         // Middle of the window horizontally
GLint midWindowY = windowHeight / 2;         // Middle of the window vertically

GLfloat fieldOfView = 15.0f; // Define our field of view (i.e. how quickly foreshortening occurs)
GLfloat near = 1.0f; // The near (Z Axis) point of our viewing frustrum (default 1.0f)
GLfloat far = 1500.0f; // The far  (Z Axis) point of our viewing frustrum (default 1500.0f)

GLint frameCount = 0; // How many frames have we drawn?
// Location of the sun (i.e. how far deep into the screen is it?)
GLfloat sunZLocation = -40.0f;

GLfloat lightPos[] = { 0.0f, 0.0f, 0.0f, 1.0f };
//GLfloat lightPos[] = { 0.0f, 0.0f, 300.0f, 1.0f };

// How fast we move (higher values mean we move and strafe faster)
GLfloat movementSpeedFactor = 0.2f;
//GLfloat movementSpeedFactor = 1.0f;
int atomos_quantidade = 0;
GLfloat caixa_tamanho = 2000;
int cont_loop_electron = 0;
int cont_loop_electron_time = 1;
GLint ghost_angulo = 0;

GLfloat max_distance_hydrogen_bond = 4.5;
GLfloat max_distance_disulfide_bond = 2.8;
GLint forca_externa_contador_max = 4;
GLint forca_externa_contador_max_t = 10;
GLint forca_externa_contador_max_hb = 5;
GLint forca_externa_contador_max_db = 1;
GLint calibration_precision = 10;
GLint calibration_precision_out = 10;
GLint calibration_minimal_distance = 4;
GLint contador_restart = 0;
GLint contador_restart_max_error = 600;
GLint contador_restart_life = 0;
GLint atomo_selecionado = -1;

GLfloat last_x;
GLfloat last_y;

// Camera rotation
GLfloat camXRot = 0.0f;
GLfloat camYRot = 0.0f;
GLfloat camZRot = 0.0f;

// Camera position
GLfloat camXPos = 0.0f;
GLfloat camYPos = 0.0f;
GLfloat camZPos = 0.0f;

// Camera movement speed
GLfloat camXSpeed = 0.0f;
GLfloat camYSpeed = 0.0f;
GLfloat camZSpeed = 0.0f;

// Atomo 1
GLfloat posx[1000];
GLfloat posy[1000];
GLfloat posz[1000];

GLfloat posx_backup[1000];
GLfloat posy_backup[1000];
GLfloat posz_backup[1000];

// Electron positions
GLfloat posEx[1000][4];
GLfloat posEy[1000][4];
GLfloat posEz[1000][4];

GLfloat posEx_HB[1000];
GLfloat posEy_HB[1000];
GLfloat posEz_HB[1000];

GLfloat electron_raio[1000];
GLfloat electron_raio_HB[1000];

GLfloat electron_y[1000][4];
GLfloat electron_z[1000][4];
GLfloat electron_y_HB[1000];
GLfloat electron_z_HB[1000];

GLint electron_arested[1000][4][2]; // Primeiro eh o atomo e segundo o numero do electron
GLint electron_arested_HB[1000]; // Primeiro no qual tem Hydrogen Bond
GLint electron_quantidade[1000];
GLfloat electron_arested_min_distance[1000][4]; // Primeiro eh o atomo e segundo o numero do electron

GLfloat velocidade_x[1000]; // velocidade
GLfloat velocidade_y[1000]; // velocidade
GLfloat velocidade_z[1000]; // velocidade
GLfloat velocidade_x_backup[1000]; // velocidade
GLfloat velocidade_y_backup[1000]; // velocidade
GLfloat velocidade_z_backup[1000]; // velocidade

GLfloat velocidade_x_original[1000]; // velocidade
GLfloat velocidade_y_original[1000]; // velocidade
GLfloat velocidade_z_original[1000]; // velocidade

char atomo_letra[1000]; // define que atomo que eh C H O...
GLint atomo_letraN[1000]; // C1=0 C2=1...
string atomo_letraL[1000]; // C1 C2...

GLfloat massa[1000];
bool atomo_base[1000];
GLint amino[1000];
GLint amino_sequencial[1000];
GLint contador_amino = 0;

GLfloat nucleo_proximity[1000];
GLfloat nucleo_proximity_free[1000];
GLfloat nucleo_proximity_HB[1000];

GLfloat percurso[500][3];
GLint percurso_contador = -1;
GLint percurso_contador_loop = 0;
GLint chain_ultimo_atomo = -1;

GLfloat calibrationMin[20][34][2][20][34]; // [Amino1][Atom1][Same|Another][Amino2][Atom2]
GLfloat calibrationMax[20][34][2][20][34]; // [Amino1][Atom1][Same|Another][Amino2][Atom2]

GLint ghost_protein_atom[3];
GLfloat ghost_x[3];
GLfloat ghost_y[3];
GLfloat ghost_z[3];

double lastTime = glfwGetTime();
string fps = "";
int contadorFrames = 0;
GLfloat pH = 7.0;

bool mousefree = false;

// Hoding any keys down?
bool holdingForward = false;
bool holdingBackward = false;
bool holdingLeftStrafe = false;
bool holdingRightStrafe = false;
bool pressionando_k = false;
bool pressionando_j = false;
bool pressionando_i = false;
bool pressionando_m = false;
//bool pressionando_o = false;
//bool pressionando_p = false;
//bool pressionando_0 = false;
//bool pressionando_l = false;
//bool pressionando_8 = false;
//bool pressionando_9 = false;
bool pressionando_control = false;

bool show_power_sphere = false;
bool show_perimetro = false;
bool show_base_line = false;
bool show_base = true;
bool rastreio = false;
bool paused = false;
bool show_colisao_tensao = false;
bool show_tensao_hb = true;
bool show_comparation = false;
bool show_comparation_speed = false;
bool show_ghost_protein = true;

GLint forca_externa_contador = 0;
GLint forca_externa_contador_t = 0;
GLint forca_externa_contador_hb = 0;
GLint forca_externa_contador_db = 0; // Dissulfite bond

// ----------------- teste sphera ----------------------
GLdouble testesphera[50][3];
GLint testesphera_contador = 0;

//GLfloat C_nucleo_proximity = 1.7;
//GLfloat C_nucleo_proximity_free = 1.7; // Van der waals angstron
//GLfloat C_electron_raio = 3.4;
//GLfloat C_massa = 12.0107;
//
//GLfloat H_nucleo_proximity = 1.2;
//GLfloat H_nucleo_proximity_free = 1.2; // Van der waals angstron
//GLfloat H_electron_raio = 2.4;
//GLfloat H_massa = 1.0079;
//
//GLfloat O_nucleo_proximity = 1.52;
//GLfloat O_nucleo_proximity_free = 1.52; // Van der waals angstron
//GLfloat O_electron_raio = 3.04;
//GLfloat O_massa = 15.099;
//
//GLfloat N_nucleo_proximity = 1.55;
//GLfloat N_nucleo_proximity_free = 1.55; // Van der waals angstron
//GLfloat N_electron_raio = 3.1;
//GLfloat N_massa = 14.0067;

GLint contador_sleep = 0;
map<int, map<string, map<string, GLfloat> > > atom_statistic;
vector<VectorProtein> ProteinOriginal;

// --------- External functions -------------

// pdb.cu
void read_pdb();
void add_pdb();

// aminos.cu
void read_pdb_amino(string procura_amino);
void add_pdb_amino();
void add_atoml(char atomo_ll);
void connect_electron(string atomo1, string atomo2);
void load_protein(string PDBID);
void conecta_eletrons(string procura_amino);
GLfloat load_protein_position(string PDBID, string model, int distance);
void compare_protein_build_MD(int c_distancia);
GLfloat compare_protein_build(int c_distancia);
void load_protein_models(string PDBID);

// util.cu
void distance_calibration();
void distance_calibration_pdb(string PDBID);
int get_amino_number(const char *amino_sigla);
int get_atom_number(const char *atomo_letra_local);
void hsvtorgb(unsigned char *r, unsigned char *g, unsigned char *b, unsigned char h, unsigned char s, unsigned char v);

//parameters.cu
void parameters_window();
//void seta_variaveis();
void handleKeypress(int theKey, int theAction);
//void setVS(string s1, GLfloat valor);
//GLfloat getVS(string s1);
//void sV(string s1, string s2, GLfloat);
//GLfloat gV(string s1, string s2);

//learn.cu
void change_properties();

// -------------------------------------------
//class vsc {
//public:
//	string s1;
//	string s2;
//	string s3;
//	string s4;
//	GLfloat f1;
//	GLfloat f2;
//	GLfloat f3;
//	GLfloat f4;
//};
//map<string, string, string> vs;

// Function to convert degrees to radians
float toRads(const float &theAngleInDegrees) {
	return theAngleInDegrees * TO_RADS;
}

// Function to convert radians to degree
float toAngle(const float &theAngleInDegrees) {
	return theAngleInDegrees * 180 / 3.141592654f;
}

// Fix negative angles
float fixNegativeAngle(const float &Angle) {
	if (Angle >= 360) {
		return Angle - 360;
	} else if (Angle < 0) {
		return 360 + Angle;
	} else {
		return Angle;
	}
}

// Function to check if OpenGL is having issues - pass it a unique string of some kind to track down where in the code it's moaning
void checkGLError(const char * errorLocation) {
	unsigned int gle = glGetError();

	if (gle != GL_NO_ERROR) {
		cout << "GL Error discovered from caller " << errorLocation << ": ";

		switch (gle) {
		case GL_INVALID_ENUM:
			cout << "Invalid enum." << endl;
			break;

		case GL_INVALID_VALUE:
			cout << "Invalid value.\n";
			break;

		case GL_INVALID_OPERATION:
			cout << "Invalid Operation.\n";
			break;

		case GL_STACK_OVERFLOW:
			cout << "Stack overflow.\n";
			break;

		case GL_STACK_UNDERFLOW:
			cout << "Stack underflow.\n";
			break;

		case GL_OUT_OF_MEMORY:
			cout << "Out of memory.\n";
			break;
		default:
			cout << "Unknown error! Enum code is: " << gle << endl;
			break;

		} // End of switch

	} // End of if error detected

} // End of chechGLError function

void add_energy() {
	for (GLint i = 0; i < atomos_quantidade; i++) {
		if (velocidade_x[i] > 0) {
			velocidade_x[i] += 0.01;
		} else {
			velocidade_x[i] -= 0.01;
		}
		if (velocidade_y[i] > 0) {
			velocidade_y[i] += 0.01;
		} else {
			velocidade_y[i] -= 0.01;
		}
		if (velocidade_z[i] > 0) {
			velocidade_z[i] += 0.01;
		} else {
			velocidade_z[i] -= 0.01;
		}
	}
}

void seta_variaveis() {
	// -------- Definindo variaveis -----------
	vF["screen"]["movementSpeedFactor"] = &movementSpeedFactor;
	vF["screen"]["caixa_tamanho"] = &caixa_tamanho;
	vF["tuning"]["max_distance_hydrogen_bond"] = &max_distance_hydrogen_bond;
	vF["tuning"]["max_distance_disulfide_bond"] = &max_distance_disulfide_bond;
	vF["contador"]["last_x"] = &last_x;
	vF["contador"]["last_y"] = &last_y;

	vI["fixo"]["atomos_quantidade"] = &atomos_quantidade;
	vI["contador"]["cont_loop_electron"] = &cont_loop_electron;
	vI["contador"]["cont_loop_electron_time"] = &cont_loop_electron_time;
	vI["tuning"]["forca_externa_contador_max"] = &forca_externa_contador_max;
	vI["tuning"]["forca_externa_contador_max_t"] = &forca_externa_contador_max_t;
	vI["tuning"]["forca_externa_contador_max_hb"] = &forca_externa_contador_max_hb;
	vI["tuning"]["forca_externa_contador_max_hb"] = &forca_externa_contador_max_db;
	vI["tuning"]["calibration_precision"] = &calibration_precision;
	vI["tuning"]["calibration_precision_out"] = &calibration_precision_out;
	vI["tuning"]["calibration_minimal_distance"] = &calibration_minimal_distance;
	vI["contador"]["contador_restart"] = &contador_restart;
	vI["tuning"]["contador_restart_max_error"] = &contador_restart_max_error;
	vI["contador"]["contador_restart_life"] = &contador_restart_life;
	vI["screen"]["atomo_selecionado"] = &atomo_selecionado;
	vI["contador"]["contador_amino"] = &contador_amino;
	vI["fixo"]["percurso_contador"] = &percurso_contador;
	vI["contador"]["percurso_contador_loop"] = &percurso_contador_loop;
	vI["contador"]["chain_ultimo_atomo"] = &chain_ultimo_atomo;
	vI["contador"]["forca_externa_contador"] = &forca_externa_contador;
	vI["contador"]["forca_externa_contador_t"] = &forca_externa_contador_t;
	vI["contador"]["forca_externa_contador_hb"] = &forca_externa_contador_hb;
	vI["contador"]["forca_externa_contador_db"] = &forca_externa_contador_db;
	vI["contador"]["testesphera_contador"] = &testesphera_contador;
	vI["contador"]["contador_sleep"] = &contador_sleep;

	vB["screen"]["mousefree"] = &mousefree;
	vB["screen"]["show_power_sphere"] = &show_power_sphere;
	vB["screen"]["show_perimetro"] = &show_perimetro;
	vB["screen"]["show_base_line"] = &show_base_line;
	vB["screen"]["show_base"] = &show_base;
	vB["screen"]["rastreio"] = &rastreio;
	vB["screen"]["paused"] = &paused;
	vB["screen"]["show_colisao_tensao"] = &show_colisao_tensao;
	vB["screen"]["show_tensao_hb"] = &show_tensao_hb;
	vB["screen"]["show_comparation"] = &show_comparation;
	vB["screen"]["show_comparation_speed"] = &show_comparation_speed;
	vB["screen"]["show_ghost_protein"] = &show_ghost_protein;

}

void initGL() {
	seta_variaveis();
	// ----- GLFW Settings -----

//	glfwDisable(GLFW_MOUSE_CURSOR); // Hide the mouse cursor

	// ----- Window and Projection Settings -----

	// Set the window title
	glfwSetWindowTitle("FY Project | Thiago Benazzi Maia");

	// Setup our viewport to be the entire size of the window
	glViewport(0, 0, (GLsizei) windowWidth, (GLsizei) windowHeight);

	// Change to the projection matrix, reset the matrix and set up our projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// The following code is a fancy bit of math that is eqivilant to calling:
	// gluPerspective(fieldOfView/2.0f, width/height , near, far);
	// We do it this way simply to avoid requiring glu.h
	GLfloat aspectRatio = (windowWidth > windowHeight) ? float(windowWidth) / float(windowHeight) : float(windowHeight) / float(windowWidth);
	GLfloat fH = tan(float(fieldOfView / 360.0f * 3.14159f)) * near;
	GLfloat fW = fH * aspectRatio;
	glFrustum(-fW, fW, -fH, fH, near, far);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// ----- OpenGL settings -----
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set out clear colour to black, full alpha
	glfwSwapInterval(1); // Lock to vertical sync of monitor (normally 60Hz, so 60fps)
	glShadeModel(GL_SMOOTH);    // Enable (gouraud) shading
	glEnable(GL_DEPTH_TEST);    // Enable depth testing
	glClearDepth(1.0f);         // Clear the entire depth of the depth buffer
	glDepthFunc(GL_LEQUAL);	// Set our depth function to overwrite if new value less than or equal to current value
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST); // Ask for nicest perspective correction
	glEnable(GL_CULL_FACE); // Do not draw polygons facing away from us
	glLineWidth(2.0f);			// Set a 'chunky' line width

	// ---- Set up OpenGL lighting -----

	// Enable lighting
	glEnable(GL_LIGHTING);

	// Ambient, diffuse and specular lighting values (note that these are ALL FOUR COMPONENT VECTORS!)
	//GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	GLfloat diffuseLight[] = { 0.7f, 0.7f, 0.7f, 1.0f };
	GLfloat specularLight[] = { 1.0f, 1.0f, 1.0f, 1.0f };

	GLint specularMagnitude = 64; // Define how "tight" our specular highlights are (larger number = smaller specular highlight). The valid range is is 1 to 128

	// Setup and enable light 0
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos); // Specify the position of the light
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight); // Specify ambient light properties
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight); // Specify diffuse light properties
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight); // Specify specular light properties
	glEnable(GL_LIGHT0);

	// Enable colour tracking of materials
	glEnable(GL_COLOR_MATERIAL);

	// Define the shininess of the material we'll use to draw things
	GLfloat materialSpecularReflectance[] = { 1.0f, 1.0f, 1.0f, 1.0f };

	// Set Material properties to follow glColor values
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	// Use our shiny material and magnitude
	glMaterialfv(GL_FRONT, GL_SPECULAR, materialSpecularReflectance);
	glMateriali(GL_FRONT, GL_SHININESS, specularMagnitude);

	// Check for any OpenGL errors (providing the location we called the function from)
	checkGLError("InitGL");

	// ---------- Aqui vai
	for (int i = 0; i < 1000; i++) {
		electron_arested_HB[i] = -1;
		for (int ii = 0; ii < 4; ii++) {
			for (int iii = 0; iii < 4; iii++) {
				electron_arested[i][ii][iii] = -1;
			}
		}
	}
	last_x = 0 - 10.0;
	last_y = 0 - 10.0;
	srand((unsigned) time(0));

	// Definindo atomos base do ghost
	ghost_protein_atom[0] = 3;
//	ghost_protein_atom[2] = 16;
	ghost_protein_atom[2] = 99;

//	distance_calibration();
//	distance_calibration_pdb("1wqc");
//	load_protein("1wqc");

	distance_calibration_pdb("1oeh");
	load_protein("1oeh");
//	load_protein_models("1oeh");

//	posx_backup = posx;
	memcpy(posx_backup, posx, sizeof(posx));
	memcpy(posy_backup, posy, sizeof(posx));
	memcpy(posz_backup, posz, sizeof(posx));
	memcpy(velocidade_x_original, velocidade_x, sizeof(velocidade_x));
	memcpy(velocidade_y_original, velocidade_y, sizeof(velocidade_y));
	memcpy(velocidade_z_original, velocidade_z, sizeof(velocidade_z));
//	add_energy();
//	compare_protein_build_MD(0);

//	 vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
//	    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
//	    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);

	//------------ teste
//int colorA[3] = [0, 0, 255];
//int colorB[3] = [255, 0, 0];

// -----------------
	// ------------ descomentar pra iniciar pausado -------------
	paused = true;
	for (GLint i = 0; i < atomos_quantidade; i++) {
		velocidade_x_backup[i] = velocidade_x[i];
		velocidade_y_backup[i] = velocidade_y[i];
		velocidade_z_backup[i] = velocidade_z[i];
		velocidade_x[i] = 0.0;
		velocidade_y[i] = 0.0;
		velocidade_z[i] = 0.0;
		cont_loop_electron_time = 999999999;
	}
	// ---------------------------------------------------------

//	load_protein("8RXN");
// ---------------------
//	read_pdb_amino("TRP");
//	read_pdb_amino("LEU");
//	read_pdb_amino("TYR");
//	read_pdb_amino("ILE");
//	read_pdb_amino("GLN");
//	read_pdb_amino("TRP");
//	read_pdb_amino("LYS");
//	read_pdb_amino("ASP");
//	read_pdb_amino("GLY");
//	read_pdb_amino("PRO");
//	read_pdb_amino("SER");
//	read_pdb_amino("ARG");

//	add_pdb_amino();
//	read_pdb();

}

// Function to move the camera the amount we've calculated in the calculateCameraMovement function
void moveCamera() {
	camXPos += camXSpeed;
	camYPos += camYSpeed;
	camZPos += camZSpeed;
}

// Function to deal with mouse position changes, called whenever the mouse cursorm moves
void handleMouseMove(int mouseX, int mouseY) {
	if (mousefree == true) {
		return;
	}
	GLfloat vertMouseSensitivity = 40.0f;
	GLfloat horizMouseSensitivity = 40.0f;

//cout << "Mouse cursor is at position (" << mouseX << ", " << mouseY << endl;

	int horizMovement = mouseX - midWindowX;
	int vertMovement = mouseY - midWindowY;

//	camXRot += vertMovement / vertMouseSensitivity;
//	camYRot += horizMovement / horizMouseSensitivity;
	camXRot += (vertMovement / vertMouseSensitivity);
	camYRot += (horizMovement / horizMouseSensitivity);

// Control looking up and down with the mouse forward/back movement
// Limit loking up to vertically up
//	if (camXRot < -90.0f) {
//		camXRot = -90.0f;
//	}
//
//// Limit looking down to vertically down
//	if (camXRot > 90.0f) {
//		camXRot = 90.0f;
//	}
//
//// Looking left and right. Keep the angles in the range -180.0f (anticlockwise turn looking behind) to 180.0f (clockwise turn looking behind)
//	if (camYRot < -180.0f) {
//		camYRot += 360.0f;
//	}
//
//	if (camYRot > 180.0f) {
//		camYRot -= 360.0f;
//	}
	// Control looking up and down with the mouse forward/back movement
	// Limit loking up to vertically up
	if (camXRot < -90.0f) {
		camXRot = -90.0f;
	}

	// Limit looking down to vertically down
	if (camXRot > 90.0f) {
		camXRot = 90.0f;
	}

	// Looking left and right. Keep the angles in the range -180.0f (anticlockwise turn looking behind) to 180.0f (clockwise turn looking behind)
	if (camYRot < -180.0f) {
		camYRot += 360.0f;
	}

	if (camYRot > 180.0f) {
		camYRot -= 360.0f;
	}

// Reset the mouse position to the centre of the window each frame
	glfwSetMousePos(midWindowX, midWindowY);
}

void restaura_posicoes() {
	printf("Restaurando posicoes\n");
	memcpy(posx, posx_backup, sizeof(posx));
	memcpy(posy, posy_backup, sizeof(posx));
	memcpy(posz, posz_backup, sizeof(posx));
	memcpy(velocidade_x, velocidade_x_original, sizeof(velocidade_x));
	memcpy(velocidade_y, velocidade_y_original, sizeof(velocidade_y));
	memcpy(velocidade_z, velocidade_z_original, sizeof(velocidade_z));
}

void pause_local() {
	paused = true;
	for (GLint i = 0; i < atomos_quantidade; i++) {
		velocidade_x_backup[i] = velocidade_x[i];
		velocidade_y_backup[i] = velocidade_y[i];
		velocidade_z_backup[i] = velocidade_z[i];
		velocidade_x[i] = 0.0;
		velocidade_y[i] = 0.0;
		velocidade_z[i] = 0.0;
		cont_loop_electron_time = 999999999;
	}
}

void continua() {
	if (paused == false) {
		pause_local();
		return;
	}
	paused = false;
	cont_loop_electron_time = 1;
	cont_loop_electron = 0;
	for (GLint i = 0; i < atomos_quantidade; i++) {
		velocidade_x[i] = velocidade_x_backup[i];
		velocidade_y[i] = velocidade_y_backup[i];
		velocidade_z[i] = velocidade_z_backup[i];
	}
}

void calculaposicoes() {
	if (pressionando_k == true) {
		posx[0] = posx[0] + 0.1f;
	}

	if (pressionando_j == true) {
		posx[0] = posx[0] - 0.1f;
	}

	if (pressionando_i == true) {
		posy[0] = posy[0] + 0.1f;
	}

	if (pressionando_m == true) {
		posy[0] = posy[0] - 0.1f;
	}
}

// Function to calculate which direction we need to move the camera and by what amount
void calculateCameraMovement() {
// Break up our movement into components along the X, Y and Z axis
	float camMovementXComponent = 0.0f;
	float camMovementYComponent = 0.0f;
	float camMovementZComponent = 0.0f;

	if (holdingForward == true) {
		// Control X-Axis movement
		float pitchFactor = cos(toRads(camXRot));
		camMovementXComponent += (movementSpeedFactor * float(sin(toRads(camYRot)))) * pitchFactor;

		// Control Y-Axis movement
		camMovementYComponent += movementSpeedFactor * float(sin(toRads(camXRot))) * -1.0f;

		// Control Z-Axis movement
		float yawFactor = float(cos(toRads(camXRot)));
		camMovementZComponent += (movementSpeedFactor * float(cos(toRads(camYRot))) * -1.0f) * yawFactor;
	}

	if (holdingBackward == true) {
		// Control X-Axis movement
		float pitchFactor = cos(toRads(camXRot));
		camMovementXComponent += (movementSpeedFactor * float(sin(toRads(camYRot))) * -1.0f) * pitchFactor;

		// Control Y-Axis movement
		camMovementYComponent += movementSpeedFactor * float(sin(toRads(camXRot)));

		// Control Z-Axis movement
		float yawFactor = float(cos(toRads(camXRot)));
		camMovementZComponent += (movementSpeedFactor * float(cos(toRads(camYRot)))) * yawFactor;
	}

	if (holdingLeftStrafe == true) {
		// Calculate our Y-Axis rotation in radians once here because we use it twice
		float yRotRad = toRads(camYRot);

		camMovementXComponent += -movementSpeedFactor * float(cos(yRotRad));
		camMovementZComponent += -movementSpeedFactor * float(sin(yRotRad));
	}

	if (holdingRightStrafe == true) {
		// Calculate our Y-Axis rotation in radians once here because we use it twice
		float yRotRad = toRads(camYRot);

		camMovementXComponent += movementSpeedFactor * float(cos(yRotRad));
		camMovementZComponent += movementSpeedFactor * float(sin(yRotRad));
	}

// After combining our movements for any & all keys pressed, assign them to our camera speed along the given axis
	camXSpeed = camMovementXComponent;
	camYSpeed = camMovementYComponent;
	camZSpeed = camMovementZComponent;

// Cap the speeds to our movementSpeedFactor (otherwise going forward and strafing at an angle is twice as fast as just going forward!)
	if (camXSpeed > movementSpeedFactor) {
		camXSpeed = movementSpeedFactor;
	}
	if (camXSpeed < -movementSpeedFactor) {
		camXSpeed = -movementSpeedFactor;
	}

	if (camYSpeed > movementSpeedFactor) {
		camYSpeed = movementSpeedFactor;
	}
	if (camYSpeed < -movementSpeedFactor) {
		camYSpeed = -movementSpeedFactor;
	}

	if (camZSpeed > movementSpeedFactor) {
		camZSpeed = movementSpeedFactor;
	}
	if (camZSpeed < -movementSpeedFactor) {
		camZSpeed = -movementSpeedFactor;
	}
}
void show_variables() {
	for (GLint i = 0; i < atomos_quantidade; i++) {
		//printf("%d)->%d A(X:%f Y:%f Z:%f)  E(X:%f Y:%f Z:%f)\n", i, electron_arested[i][0][0], posx[i], posy[i], posz[i], posEx[i][0], posEy[i][0], posEz[i][0]);
		printf("%d-%c)->%d %d %d %d\n", i, atomo_letra[i], electron_arested[i][0][0], electron_arested[i][1][0], electron_arested[i][2][0], electron_arested[i][3][0]);
	}
	printf("----------------------------------\n");
}

void ativa_desativa_comparation() {
	if (show_comparation) {
		show_comparation = false;
	} else {
		show_comparation = true;
	}
}

void ativa_desativa_comparation_speed() {
	if (show_comparation_speed) {
		show_comparation_speed = false;
	} else {
		show_comparation_speed = true;
	}
}

void ativa_desativa_perimetro() {
	if (show_perimetro) {
		show_perimetro = false;
	} else {
		show_perimetro = true;
	}
}

void ativa_desativa_rastreio() {
	if (rastreio) {
		rastreio = false;
	} else {
		rastreio = true;
	}
}

void ativa_desativa_forca() {
	if (show_power_sphere) {
		show_power_sphere = false;
	} else {
		show_power_sphere = true;
	}
}

void ativa_desativa_ghost_protein() {
	if (show_ghost_protein) {
		show_ghost_protein = false;
	} else {
		show_ghost_protein = true;
	}
}

void ativa_desativa_base() {
	if (show_base) {
		show_base = false;
	} else {
		show_base = true;
	}
}

void ativa_desativa_base_line() {
	if (show_base_line) {
		show_base_line = false;
	} else {
		show_base_line = true;
	}
}

void ativa_desativa_colisao_tensao() {
	if (show_colisao_tensao) {
		show_colisao_tensao = false;
	} else {
		show_colisao_tensao = true;
	}
}

void ativa_desativa_tensao_hb() {
	if (show_tensao_hb) {
		show_tensao_hb = false;
	} else {
		show_tensao_hb = true;
	}
}

void camera_position(GLfloat px, GLfloat py, GLfloat pz, GLfloat rx, GLfloat ry, GLfloat rz) {
//	camXRot = rx;
//	camYRot = ry;
//	camZRot = rz;
	camXRot = rx;
	camYRot = ry;
	camZRot = rz;

	camXPos = px;
	camYPos = py;
	camZPos = pz + 10;
}

void rem_energy() {
	for (GLint i = 0; i < atomos_quantidade; i++) {
		if (velocidade_x[i] > 0) {
			velocidade_x[i] -= 0.01;
		}
		if (velocidade_x[i] < 0) {
			velocidade_x[i] += 0.01;
		}
		if (velocidade_y[i] > 0) {
			velocidade_y[i] -= 0.01;
		}
		if (velocidade_y[i] < 0) {
			velocidade_y[i] += 0.01;
		}
		if (velocidade_z[i] > 0) {
			velocidade_z[i] -= 0.01;
		}
		if (velocidade_z[i] < 0) {
			velocidade_z[i] += 0.01;
		}
	}
}

void add_atom(char atomo_ll) {
	if (atomos_quantidade == 998) {
		printf("Quantidade de atomos chegou no limite \n");
		return;
	}

	if (atomo_ll == 'C') {
		atomo_letra[atomos_quantidade] = 'C';
		nucleo_proximity[atomos_quantidade] = 1.7;
		nucleo_proximity_free[atomos_quantidade] = 3.4; // Van der waals angstron
		electron_quantidade[atomos_quantidade] = 4;
		electron_raio[atomos_quantidade] = 3.4;
		massa[atomos_quantidade] = 12.0107;
	} else if (atomo_ll == 'H') {
		atomo_letra[atomos_quantidade] = 'H';
		nucleo_proximity[atomos_quantidade] = 1.2;
		nucleo_proximity_free[atomos_quantidade] = 2.4;  // Van der waals angstron
		electron_quantidade[atomos_quantidade] = 1;
		electron_raio[atomos_quantidade] = 2.4;
		electron_raio_HB[atomos_quantidade] = 2.4;
		massa[atomos_quantidade] = 1.0079;
		nucleo_proximity_HB[atomos_quantidade] = 1.2;

//		electron_y_HB[atomos_quantidade] = 270;
	} else if (atomo_ll == 'O') {
		atomo_letra[atomos_quantidade] = 'O';
		nucleo_proximity[atomos_quantidade] = 1.52;
		nucleo_proximity_free[atomos_quantidade] = 3.04;  // Van der waals angstron
		electron_quantidade[atomos_quantidade] = 2;
		electron_raio[atomos_quantidade] = 3.04;
		electron_raio_HB[atomos_quantidade] = 3.04;
		massa[atomos_quantidade] = 15.099;
		nucleo_proximity_HB[atomos_quantidade] = 1.52;

//		electron_y_HB[atomos_quantidade] = 90;
	} else if (atomo_ll == 'N') {
		atomo_letra[atomos_quantidade] = 'N';
//		nucleo_proximity[atomos_quantidade] = 1.55; // original
		nucleo_proximity[atomos_quantidade] = 1.46; // average collected
		nucleo_proximity_free[atomos_quantidade] = 3.1;  // Van der waals angstron
		electron_quantidade[atomos_quantidade] = 4;
		electron_raio[atomos_quantidade] = 2.92;
		electron_raio_HB[atomos_quantidade] = 3.1;
		massa[atomos_quantidade] = 14.0067;
		nucleo_proximity_HB[atomos_quantidade] = 1.55;
	} else if (atomo_ll == 'S') {
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

void add_atom_teste() {
	add_atom('C');

	posx[atomos_quantidade - 1] = 0.0;
	posy[atomos_quantidade - 1] = 0.0;

	posz[atomos_quantidade - 1] = 0.0;

	add_atom('O');
	posx[atomos_quantidade - 1] = posx[atomos_quantidade - 2] - 5;
	posy[atomos_quantidade - 1] = posy[atomos_quantidade - 2] - 3;
	posz[atomos_quantidade - 1] = 0.0;
	electron_arested[atomos_quantidade - 1][0][0] = atomos_quantidade - 2;
	electron_arested[atomos_quantidade - 1][0][1] = 1;
	electron_arested[atomos_quantidade - 2][2][0] = atomos_quantidade - 1;
	electron_arested[atomos_quantidade - 2][2][1] = 0;

	add_atom('H');
	posx[atomos_quantidade - 1] = posx[atomos_quantidade - 2];
	posy[atomos_quantidade - 1] = posy[atomos_quantidade - 2] + 5;
	posz[atomos_quantidade - 1] = 0.0;
	electron_arested[atomos_quantidade - 1][0][0] = atomos_quantidade - 3;
	electron_arested[atomos_quantidade - 1][0][1] = 3;
	electron_arested[atomos_quantidade - 3][3][0] = atomos_quantidade - 1;
	electron_arested[atomos_quantidade - 3][3][1] = 0;
	velocidade_y[atomos_quantidade - 1] = 0.05;

}

void add_chain() {
//	if (contador_sleep == 0) {
//		read_pdb_amino("TRP");
//	} else if (contador_sleep == 1) {
//		read_pdb_amino("TRP");
//	}
// GLN TRP
// ---------------------------------------------------
	if (contador_sleep == 0) {
		read_pdb_amino("ASN");
	} else if (contador_sleep == 1) {
		read_pdb_amino("LEU");
	} else if (contador_sleep == 2) {
		read_pdb_amino("TYR");
	} else if (contador_sleep == 3) {
		read_pdb_amino("ILE");
	} else if (contador_sleep == 4) {
		read_pdb_amino("GLN");
	} else if (contador_sleep == 5) {
		read_pdb_amino("TRP");
	} else if (contador_sleep == 6) {
		read_pdb_amino("LYS");
	} else if (contador_sleep == 7) {
		read_pdb_amino("ASP");
	} else if (contador_sleep == 8) {
		read_pdb_amino("GLY");
	} else if (contador_sleep == 9) {
		read_pdb_amino("PRO");
	} else if (contador_sleep == 10) {
		read_pdb_amino("SER");
	} else if (contador_sleep == 11) {
		read_pdb_amino("ARG");
	}
	contador_sleep++;
//	if(contador_sleep == 1){
//		read_pdb_amino("TYR");
//	} if(contador_sleep == 2){
//		read_pdb_amino("TYR");
//	} if(contador_sleep == 3){
//		read_pdb_amino("TYR");
//	}
//	contador_sleep++;
//	add_atom_teste();
// Unico -------------------
//	add_glu(true, true);
// -------------------------

// sequencia que precisa
//	add_leu(false, false);
//	add_tyr(false, false);
//	add_ile(false, false);
//	add_gln(false, false);
//	add_trp(false, false);
//	add_leu(false, false);
//	add_lys(false, false);
// --------------------

//	add_gln(false, false);
//	add_gln(false, false);
//	add_gln(false, false);
//	add_gln(false, false);
//	add_gln(false, false);
//	add_gln(false, false);
//	add_gln(false, false);
//	add_gln(false, false);

//	add_tyr(false, false);
//	add_tyr(false, false);
//	add_ile(false, false);
//	add_gln(false, false);
//	add_trp(false, false);
//	add_trp(false, false);
//	add_leu(false, false);
//	add_lys(false, false);

//	add_tyr(false, false);
//	add_gln(false, false);

//	add_ile(false, false);
//	add_trp(false, false);

//	add_glu(false, false);
//	add_glu(false, false);
//	add_glu(false, false);
//	add_glu(false, false);
//	add_glu(false, false);
//	add_glu(false, false);
//	add_glu(false, false);
//	add_glu(false, false);
//	add_glu(false, false);
//	add_glu(false, true);
}

// Function to draw a grid of lines
void drawGround() {
	GLfloat extent = 600.0f; // How far on the Z-Axis and X-Axis the ground extends
	GLfloat stepSize = 50.0f;  // The size of the separation between points
	GLfloat groundLevel = -50.0f;   // Where on the Y-Axis the ground is drawn

// Set colour to white
	glColor3ub(255, 255, 255);

// Draw our ground grid
	glBegin(GL_LINES);
	for (GLint loop = -extent; loop < extent; loop += stepSize) {
// Draw lines along Z-Axis
		glVertex3f(loop, groundLevel, extent);
		glVertex3f(loop, groundLevel, -extent);

// Draw lines across X-Axis
		glVertex3f(-extent, groundLevel, loop);
		glVertex3f(extent, groundLevel, loop);
	}
	glEnd();
}

GLfloat show_distance(GLfloat x1, GLfloat x2, GLfloat y1, GLfloat y2, GLfloat z1, GLfloat z2) {
	return abs(sqrt(pow((x2 - x1), 2.0) + pow((y2 - y1), 2.0) + pow((z2 - z1), 2.0)));
}

//bool angle_colision(GLint atomo1, GLint atomo2) {
//	GLfloat d_atomo_atomo = show_distance(posx[atomo1], posx[atomo2], posy[atomo1], posy[atomo2], posz[atomo1], posz[atomo2]);
//	GLfloat atomo_electron = electron_raio[atomo1] + electron_raio[atomo2];
//	if (abs(d_atomo_atomo - atomo_electron) < collision_angleE) {
//		return true;
//	}
//	return false;
//}

bool checa_between(GLfloat numero, GLfloat comeco, GLfloat fim, GLfloat tensao) {
	return false;
}

void collision3D(GLint i, GLint ii, string tipo) {
//	printf("C1 %s %f\n", tipo.c_str(), show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]),calibrationMax[amino[i]][atomo_letraN[i]][0][amino[ii]][atomo_letraN[ii]]);
//	printf("C1 %s Dist:%f Cal:%f i:%d ii:%d am1:%d am2:%d | N1C1max:%f\n", tipo.c_str(), show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]), calibrationMax[amino[i]][atomo_letraN[i]][0][amino[ii]][atomo_letraN[ii]], i, ii, amino[i], amino[ii],calibrationMax[0][11][0][0][0]);
//	if (tipo == "tensao") {
//		GLint contador_contrario = 0;
//		if ((velocidade_x[i] > 0 && velocidade_x[ii] < 0) || (velocidade_x[ii] > 0 && velocidade_x[i] < 0)) {
////			printf("X contrario\n");
//			contador_contrario++;
//		}
//		if ((velocidade_y[i] > 0 && velocidade_y[ii] < 0) || (velocidade_y[ii] > 0 && velocidade_y[i] < 0)) {
////			printf("Y contrario\n");
//			contador_contrario++;
//		}
//		if ((velocidade_z[i] > 0 && velocidade_z[ii] < 0) || (velocidade_z[ii] > 0 && velocidade_z[i] < 0)) {
////			printf("Z contrario\n");
//			contador_contrario++;
//		}
//		if (contador_contrario > 1) {
//			printf("Tensao %d %f %f %f %f %f %f\n", contador_contrario, velocidade_x[i], velocidade_x[ii], velocidade_y[i], velocidade_y[ii], velocidade_z[i], velocidade_z[ii]);
//		} else {
//			return;
//		}
//	} else if (tipo == "colisao") {
//		printf("Colisao %f %f %f %f %f %f\n");
//	} else if (tipo == "electron_nucleo") {
//		printf("Electron nucleo %f %f %f %f %f %f\n");
//	}

// In this part for 3d colision we are mostly using as a base the code from Thomas Smid (M.Sc. Physics, Ph.D. Astronomy)
//	GLfloat R = 1;   //(restitution coefficient)  between 0 and 1 (1=perfectly elastic collision)
	GLfloat m1 = massa[i];  // (mass of ball 1)
	GLfloat m2 = massa[ii];  // (mass of ball 2)
	GLfloat r1 = nucleo_proximity_free[i];  // (radius of ball 1)
	GLfloat r2 = nucleo_proximity_free[ii];  // (radius of ball 2)
	GLfloat x1 = posx[i];  // (x-coordinate of ball 1)
	GLfloat y1 = posy[i];  // (y-coordinate of ball 1)
	GLfloat z1 = posz[i];  // (z-coordinate of ball 1)
	GLfloat x2 = posx[ii];  // (x-coordinate of ball 2)
	GLfloat y2 = posy[ii];  // (y-coordinate of ball 2)
	GLfloat z2 = posz[ii]; // (z-coordinate of ball 2)
	GLfloat vx1 = velocidade_x[i]; // (velocity x-component of ball 1)
	GLfloat vy1 = velocidade_y[i]; // (velocity y-component of ball 1)
	GLfloat vz1 = velocidade_z[i]; // (velocity z-component of ball 1)
	GLfloat vx2 = velocidade_x[ii]; // (velocity x-component of ball 2)
	GLfloat vy2 = velocidade_y[ii]; // (velocity y-component of ball 2)
	GLfloat vz2 = velocidade_z[ii]; // (velocity z-component of ball 2)
//	GLint error = 0;    // (0: no error: balls do not collide 2: initial positions impossible (balls overlap))
//	GLfloat pi, r12, m21, d, v, theta2, phi2, st, ct, sp, cp, vx1r, vy1r, vz1r, fvz1r, thetav, phiv, dr, alpha, beta, sbeta, cbeta, dc, sqs, t, a, dvz2, vx2r, vy2r, vz2r, x21, y21, z21, vx21, vy21, vz21, vx_cm, vy_cm, vz_cm;
	GLfloat pi, r12, m21, d, v, theta2, phi2, st, ct, sp, cp, vx1r, vy1r, vz1r, fvz1r, thetav, phiv, dr, alpha, beta, sbeta, cbeta, t, a, dvz2, vx2r, vy2r, vz2r, x21, y21, z21, vx21, vy21, vz21;

//	if (show_distance(posx[0], posx[1], posy[0], posy[1], posz[0], posz[2]) < calibrationMin[0][0][0][0][11]) {
//		printf("C1-N1 %f %fg\n", show_distance(posx[0], posx[1], posy[0], posy[1], posz[0], posz[2]), calibrationMin[0][0][0][0][11]);
//	}

	if (tipo == "tensao") {
//		printf("Tensao\n");
		x1 = posx[i] - (2 * (posx[i] - posx[ii]));
		y1 = posy[i] - (2 * (posy[i] - posy[ii]));
		z1 = posz[i] - (2 * (posz[i] - posz[ii]));
		r1 = show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) / 2;
		r2 = r1;

		if (show_colisao_tensao) {
			glColor3ub(0, 255, 0);
			glBegin(GL_LINES);
			glVertex3f(posx[i], posy[i], posz[i]);
			glVertex3f(posx[ii], posy[ii], posz[ii]);
			glEnd();
		}
	} else if (tipo == "tensaoHB") {
//		printf("Tensao HB\n");
		x1 = posx[i] - (2 * (posx[i] - posx[ii]));
		y1 = posy[i] - (2 * (posy[i] - posy[ii]));
		z1 = posz[i] - (2 * (posz[i] - posz[ii]));
		r1 = show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) / 2;
		r2 = r1;

		if (show_tensao_hb) {
//			glColor3ub(255, 255, 0);
			glColor3ub(0, 128, 255);
			glBegin(GL_LINES);
			glVertex3f(posx[i], posy[i], posz[i]);
			glVertex3f(posx[ii], posy[ii], posz[ii]);
			glEnd();
		}
	} else if (tipo == "tensaoDB") {
		//		printf("Tensao HB\n");
		x1 = posx[i] - (2 * (posx[i] - posx[ii]));
		y1 = posy[i] - (2 * (posy[i] - posy[ii]));
		z1 = posz[i] - (2 * (posz[i] - posz[ii]));
		r1 = show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) / 2;
		r2 = r1;

		if (show_tensao_hb) {
			glColor3ub(255, 128, 0);
			glBegin(GL_LINES);
			glVertex3f(posx[i], posy[i], posz[i]);
			glVertex3f(posx[ii], posy[ii], posz[ii]);
			glEnd();
		}
	} else if (tipo == "electron_nucleo") {
//		printf("Colisao\n");
		r1 = show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) - r2;

	} else {
//		printf("Fucao Colisao\n");
		if (show_colisao_tensao) {
			glColor3ub(255, 0, 0);
			glBegin(GL_LINES);
			glVertex3f(posx[i], posy[i], posz[i]);
			glVertex3f(posx[ii], posy[ii], posz[ii]);
			glEnd();
		}
	}
//     **** initialize some variables ****
	pi = acos(-1.0E0);
//	error = 0;
	r12 = r1 + r2;
	m21 = m2 / m1;
	x21 = x2 - x1;
	y21 = y2 - y1;
	z21 = z2 - z1;
	vx21 = vx2 - vx1;
	vy21 = vy2 - vy1;
	vz21 = vz2 - vz1;

//	vx_cm = (m1 * vx1 + m2 * vx2) / (m1 + m2);
//	vy_cm = (m1 * vy1 + m2 * vy2) / (m1 + m2);
//	vz_cm = (m1 * vz1 + m2 * vz2) / (m1 + m2);

//     **** calculate relative distance and relative speed ***
	d = sqrt(x21 * x21 + y21 * y21 + z21 * z21);
//	if (tipo == "tensao") {
//		d = r1 + r2;
//	}
	v = sqrt(vx21 * vx21 + vy21 * vy21 + vz21 * vz21);

//     **** return if distance between balls smaller than sum of radii ****
	if (d > r12) {
//		error = 2;
//		printf("X2 %f  X1 %f\n",x2,x1);
//		printf("Y2 %f  Y1 %f\n",y2,y1);
//		printf("Z2 %f  Z1 %f\n",z2,z1);
//		printf("X21 %f  Y21 %f  Z21 %f\n", x21, y21, z21);
//		printf("Minha dist pos: %f\n",show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]));
//		printf("Minha dist XYZ: %f\n",show_distance(x1, x2, y1, y2, z1, z2));
//		printf("Saindo 1: d:%f r12:%f r1:%f x21:%f\n", d, r12, r1,x21);
		return;
	}

//     **** return if relative speed = 0 ****
	if (v == 0) {
//		error = 1;
//		printf("Saindo 2\n");
		return;
	}

//     **** shift coordinate system so that ball 1 is at the origin ***
	x2 = x21;
	y2 = y21;
	z2 = z21;

//     **** boost coordinate system so that ball 2 is resting ***
	vx1 = -vx21;
	vy1 = -vy21;
	vz1 = -vz21;

//     **** find the polar coordinates of the location of ball 2 ***
	theta2 = acos(z2 / d);
	if (x2 == 0 && y2 == 0) {
		phi2 = 0;
	} else {
		phi2 = atan2(y2, x2);
	}
	st = sin(theta2);
	ct = cos(theta2);
	sp = sin(phi2);
	cp = cos(phi2);

//     **** express the velocity vector of ball 1 in a rotated coordinate
//          system where ball 2 lies on the z-axis ******
	vx1r = ct * cp * vx1 + ct * sp * vy1 - st * vz1;
	vy1r = cp * vy1 - sp * vx1;
	vz1r = st * cp * vx1 + st * sp * vy1 + ct * vz1;
	fvz1r = vz1r / v;
	if (fvz1r > 1) {
		fvz1r = 1;
	}   // fix for possible rounding errors
	else if (fvz1r < -1) {
		fvz1r = -1;
	}
	thetav = acos(fvz1r);
	if (vx1r == 0 && vy1r == 0) {
		phiv = 0;
	} else {
		phiv = atan2(vy1r, vx1r);
	}

//     **** calculate the normalized impact parameter ***
	dr = d * sin(thetav) / r12;

//     **** return old positions and velocities if balls do not collide ***
	if (thetav > pi / 2 || fabs(dr) > 1) {
		x2 = x2 + x1;
		y2 = y2 + y1;
		z2 = z2 + z1;
		vx1 = vx1 + vx2;
		vy1 = vy1 + vy2;
		vz1 = vz1 + vz2;
//		error = 1;
		return;
	}

//     **** calculate impact angles if balls do collide ***
	alpha = asin(-dr);
	beta = phiv;
	sbeta = sin(beta);
	cbeta = cos(beta);

//     **** calculate time to collision ***
	t = (d * cos(thetav) - r12 * sqrt(1 - dr * dr)) / v;

//     **** update positions and reverse the coordinate shift ***
	x2 = x2 + vx2 * t + x1;
	y2 = y2 + vy2 * t + y1;
	z2 = z2 + vz2 * t + z1;
	x1 = (vx1 + vx2) * t + x1;
	y1 = (vy1 + vy2) * t + y1;
	z1 = (vz1 + vz2) * t + z1;

//  ***  update velocities ***

	a = tan(thetav + alpha);

	dvz2 = 2 * (vz1r + a * (cbeta * vx1r + sbeta * vy1r)) / ((1 + a * a) * (1 + m21));

	vz2r = dvz2;
	vx2r = a * cbeta * dvz2;
	vy2r = a * sbeta * dvz2;
	vz1r = vz1r - m21 * vz2r;
	vx1r = vx1r - m21 * vx2r;
	vy1r = vy1r - m21 * vy2r;

//     **** rotate the velocity vectors back and add the initial velocity
//           vector of ball 2 to retrieve the original coordinate system ****

	velocidade_x[i] = ct * cp * vx1r - sp * vy1r + st * cp * vz1r + vx2;
	velocidade_y[i] = ct * sp * vx1r + cp * vy1r + st * sp * vz1r + vy2;
	velocidade_z[i] = ct * vz1r - st * vx1r + vz2;
	velocidade_x[ii] = ct * cp * vx2r - sp * vy2r + st * cp * vz2r + vx2;
	velocidade_y[ii] = ct * sp * vx2r + cp * vy2r + st * sp * vz2r + vy2;
	velocidade_z[ii] = ct * vz2r - st * vx2r + vz2;

//	if (tipo == "colisao") {
//			posx[i] = posx[i] + 0.05;
//			posy[i] = posy[i] - 0.05;
//			posz[i] = posz[i] + 0.05;
//	}

	return;
}

void drawP(GLfloat oy, GLfloat oz, GLfloat ox) {
}

void drawE(GLfloat oy, GLfloat oz, GLfloat ox, GLint atomN, GLint electronN) {
//	printf("De %f %f %d %d\n", oy, oz, atomN, electronN);
	GLfloat cosseno = cos(toRads(oy)) * ox;
	GLfloat seno = sin(toRads(oy)) * ox;
	GLfloat cossenoZ = cos(toRads(oz)) * cosseno;
	GLfloat senoZ = sin(toRads(oz)) * cosseno;
	glTranslatef(cossenoZ, seno, senoZ);
	posEx[atomN][electronN] = posx[atomN] + cossenoZ;
	posEy[atomN][electronN] = posy[atomN] + seno;
	posEz[atomN][electronN] = posz[atomN] + senoZ;
//	printf("De %f %f\n", posx[0], posEx[0][0]);
	glutSolidSphere(0.05f, 6, 6);
	glTranslatef(-cossenoZ, -seno, -senoZ);
}

void drawE_HB(GLfloat oy, GLfloat oz, GLfloat ox, GLint atomN) {
//	printf("De %f %f %d %d\n", oy, oz, atomN, electronN);
	GLfloat cosseno = cos(toRads(oy)) * ox;
	GLfloat seno = sin(toRads(oy)) * ox;
	GLfloat cossenoZ = cos(toRads(oz)) * cosseno;
	GLfloat senoZ = sin(toRads(oz)) * cosseno;
	glTranslatef(cossenoZ, seno, senoZ);
	posEx_HB[atomN] = posx[atomN] + cossenoZ;
	posEy_HB[atomN] = posy[atomN] + seno;
	posEz_HB[atomN] = posz[atomN] + senoZ;
//	printf("De %f %f\n", posx[0], posEx[0][0]);
	glutSolidSphere(0.2f, 6, 6);
	glTranslatef(-cossenoZ, -seno, -senoZ);
}

// Function to draw our spheres and position the light source
void drawScene() {

// Clear the screen and depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
// Reset the matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
// Move the camera to our location in space
	glRotatef(camXRot, 1.0f, 0.0f, 0.0f); // Rotate our camera on the x-axis (looking up and down)
	glRotatef(camYRot, 0.0f, 1.0f, 0.0f); // Rotate our camera on the  y-axis (looking left and right)
	glTranslatef(-camXPos, -camYPos, -camZPos); // Translate the modelviewm matrix to the position of our camera
// Move everything "into" the screen (i.e. move 300 units along the Z-axis into the screen) so that all positions are now relative to the location of the sun
	glTranslatef(0.0f, 0.0f, sunZLocation);

	glColor3ub(255, 255, 0);
	// ------------ Desenhando cubo
//	glutWireCube(2.0 * caixa_tamanho);

	// Testando sphera
//	for (GLint ii = 0; ii < testesphera_contador; ii++) {
//		glColor3ub(255, 255, 0);
//		glTranslatef(testesphera[ii][0], testesphera[ii][1], testesphera[ii][2]);
//		glutSolidSphere(0.1, 6, 6);
//		glTranslatef(-testesphera[ii][0], -testesphera[ii][1], -testesphera[ii][2]);
//	}

	// --------------------------

// Aidicionando ribosome
//	add_pdb();

// --------- Testando erros --------------
	/*
	 if(contador_restart > contador_restart_max_error){
	 printf("Ultrapassou limite de erro %d\n",contador_restart);
	 restaura_posicoes();
	 change_properties();
	 contador_restart_life = 0;
	 }
	 contador_restart_life++;
	 */

	if (show_comparation) {
		compare_protein_build_MD(1);
	}

	contador_restart = 0;
	for (GLint i = 0; i < atomos_quantidade; i++) {
		if (posy[i] > caixa_tamanho) {
			velocidade_y[i] *= -1;
			if (posy[i] > caixa_tamanho + 0.5) {
				posy[i] = caixa_tamanho;
			}
		} else if (posy[i] < caixa_tamanho * -1) {
			velocidade_y[i] *= -1;
			if (posy[i] < (caixa_tamanho * -1) - 0.5) {
				posy[i] = caixa_tamanho * -1;
			}
		}
		if (posx[i] > caixa_tamanho) {
			velocidade_x[i] *= -1;
			if (posx[i] > caixa_tamanho + 0.5) {
				posx[i] = caixa_tamanho;
			}
		} else if (posx[i] < caixa_tamanho * -1) {
			velocidade_x[i] *= -1;
			if (posx[i] < (caixa_tamanho * -1) - 0.5) {
				posx[i] = caixa_tamanho * -1;
			}
		}
		if (posz[i] > caixa_tamanho) {
			velocidade_z[i] *= -1;
			if (posz[i] > caixa_tamanho + 0.5) {
				posz[i] -= caixa_tamanho;
			}
		} else if (posz[i] < caixa_tamanho * -1) {
			velocidade_z[i] *= -1;
			if (posz[i] < (caixa_tamanho * -1) - 0.5) {
				posz[i] = caixa_tamanho * -1;
			}
		}

// ===================== Movendo ======================
		posx[i] += velocidade_x[i];
		posy[i] += velocidade_y[i];
		posz[i] += velocidade_z[i];
// ================ Checa colisao de nucleo ================
		GLint same_amino = 0; //0 same / 1 different
		for (GLint ii = 0; ii < atomos_quantidade; ii++) {
			if (i != ii) {
				same_amino = 0;
				if (amino_sequencial[i] != amino_sequencial[ii]) {
					same_amino = 1;
				}
				if (calibrationMin[amino[i]][atomo_letraN[i]][same_amino][amino[ii]][atomo_letraN[ii]] == 0.0) {
					printf("Nao existe valor para %d %d %s %d %s %f\n", same_amino, amino[i], atomo_letraL[i].c_str(), amino[ii], atomo_letraL[ii].c_str(), show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]));
//					printf("Nao existe valor para %s\n", amino[i]);
				}
				// Comentado pra resolve o problema de atomos se alinhando
				if (electron_arested[i][0][0] == ii || electron_arested[i][1][0] == ii || electron_arested[i][2][0] == ii || electron_arested[i][3][0] == ii) {
					if (show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) <= calibrationMin[amino[i]][atomo_letraN[i]][same_amino][amino[ii]][atomo_letraN[ii]]) {
						collision3D(i, ii, "colisao"); // Comentado pra resolve o problema de atomos se alinhando
					}
				} else if (show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) <= calibrationMin[amino[i]][atomo_letraN[i]][same_amino][amino[ii]][atomo_letraN[ii]]) {
					if (electron_arested[i][0][0] == ii || electron_arested[i][1][0] == ii || electron_arested[i][2][0] == ii || electron_arested[i][3][0] == ii) {
						collision3D(i, ii, "colisao"); // Comentado pra resolve o problema de atomos se alinhando
					} else if (forca_externa_contador == forca_externa_contador_max && show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) < calibration_minimal_distance) {
//						printf("Aqui %d %d\n",atomo_letraN[i],atomo_letraN[ii]);
						collision3D(i, ii, "colisao");
					}
				}
			}
		}

		// ================ Checa tensao de nucleo ================
		for (GLint ii = 0; ii < atomos_quantidade; ii++) {
			if (i != ii) {
				same_amino = 0;
				if (amino_sequencial[i] == amino_sequencial[ii]) {
					if (show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) >= calibrationMax[amino[i]][atomo_letraN[i]][same_amino][amino[ii]][atomo_letraN[ii]]) {
						if (forca_externa_contador_t == forca_externa_contador_max_t && show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) < calibration_minimal_distance) {
							collision3D(i, ii, "tensao");
						} else if (electron_arested[i][0][0] == ii || electron_arested[i][1][0] == ii || electron_arested[i][2][0] == ii || electron_arested[i][3][0] == ii) {
							collision3D(i, ii, "tensao");
						}
					}
				} else if (electron_arested[i][0][0] == ii || electron_arested[i][1][0] == ii || electron_arested[i][2][0] == ii || electron_arested[i][3][0] == ii) {
					same_amino = 1;
					if (show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) >= calibrationMin[amino[i]][atomo_letraN[i]][same_amino][amino[ii]][atomo_letraN[ii]] + 0.04) {
						collision3D(i, ii, "tensao");
					}
//				} else if (atomo_letra[i] == 'S' && atomo_letra[ii] == 'S' && forca_externa_contador_db == forca_externa_contador_max_db && calibrationMin[amino[i]][atomo_letraN[i]][same_amino][amino[ii]][atomo_letraN[ii]] < max_distance_disulfide_bond) {
				} else if (atomo_letra[i] == 'S' && atomo_letra[ii] == 'S' && forca_externa_contador_db == forca_externa_contador_max_db && show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) < max_distance_disulfide_bond) {
//					printf("Possivel DB %s %s  CalMin %f \n",atomo_letra[i], atomo_letra[ii], calibrationMin[amino[i]][atomo_letraN[i]][same_amino][amino[ii]][atomo_letraN[ii]]);
//					same_amino = 1;
					if (show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) >= calibrationMin[amino[i]][atomo_letraN[i]][same_amino][amino[ii]][atomo_letraN[ii]]) {
//						printf("CMDB %f %d %d\n",calibrationMin[amino[i]][atomo_letraN[i]][same_amino][amino[ii]][atomo_letraN[ii]],atomo_letraN[i],atomo_letraN[ii]);
						collision3D(i, ii, "tensaoDB");
					}
				} else if (forca_externa_contador_hb == forca_externa_contador_max_hb && amino_sequencial[i] != amino_sequencial[ii] && calibrationMin[amino[i]][atomo_letraN[i]][same_amino][amino[ii]][atomo_letraN[ii]] < max_distance_hydrogen_bond && ((atomo_letraN[i] == 11 && atomo_letraN[ii] == 15) || (atomo_letraN[i] == 15 && atomo_letraN[ii] == 11))) {
					// Hrydrogen bonds
					same_amino = 1;
					if (show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) < max_distance_hydrogen_bond) {
						if (show_distance(posx[i], posx[ii], posy[i], posy[ii], posz[i], posz[ii]) >= calibrationMin[amino[i]][atomo_letraN[i]][same_amino][amino[ii]][atomo_letraN[ii]] + 0.04) {
							collision3D(i, ii, "tensaoHB");
						}
					}
				} else if (electron_arested[i][0][0] == ii || electron_arested[i][1][0] == ii || electron_arested[i][2][0] == ii || electron_arested[i][3][0] == ii) {
					same_amino = 1;
				}
			}
		}

// ===================== Nucleo ======================
		if (show_comparation) {
			int temp = 0;
			if (atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MD"] >= atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MIN"] && atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MD"] <= atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MAX"]) {
			} else if (atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MD"] < atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MIN"]) {
				temp = atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MD"] - atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MIN"];
				contador_restart += (temp * -1);
				temp = temp * 7;
				contador_restart++;
			} else if (atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MD"] > atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MAX"]) {
				temp = atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MD"] - atom_statistic[amino_sequencial[i]][atomo_letraL[i]]["MAX"];
				contador_restart += temp;
				temp = temp * 7;
			}
			unsigned char r, g, b;
			hsvtorgb(&r, &g, &b, 120 - temp, 255, 255);
			glColor3ub(r, g, b);
		} else if (show_comparation_speed) {
//			int temp = fabs((velocidade_x_backup[i]) + fabs(velocidade_y_backup[i]) + fabs(velocidade_z_backup[i])) * 3000;
			int temp = fabs((velocidade_x[i]) + fabs(velocidade_y[i]) + fabs(velocidade_z[i])) * 3000;
			unsigned char r, g, b;
			hsvtorgb(&r, &g, &b, 120 - temp, 255, 255);
			glColor3ub(r, g, b);

		} else {
			if (atomo_letra[i] == 'C') {
				glColor3ub(126, 152, 150);
			} else if (atomo_letra[i] == 'H') {
				glColor3ub(233, 233, 233);
			} else if (atomo_letra[i] == 'O') {
				glColor3ub(255, 0, 0);
			} else if (atomo_letra[i] == 'N') {
				glColor3ub(78, 179, 255);
			} else if (atomo_letra[i] == 'S') {
				glColor3ub(255, 128, 0);
			} else {
				glColor3ub(200, 67, 55);
			}
		}
		if (vp["global"]["atomo_selecionado"][""][""][""][""] == i) {
			glColor3ub(255, 255, 255);
		}
		glTranslatef(posx[i], posy[i], posz[i]);
		if (show_base_line == false) {
			if (show_base == false || (show_base == true && atomo_base[i] == true)) {
				if (show_perimetro) {
					glutSolidSphere(nucleo_proximity_free[i], 30, 30);
				} else {
//					glutSolidCube(nucleo_proximity[i] * .50f);
					glutSolidSphere(nucleo_proximity[i] * 0.2f, 20, 20);
				}
			}
		}
		if (rastreio) {
			// Linha teto
			glBegin(GL_LINES);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, caixa_tamanho - posy[i], 0.0);
			glEnd();

			// Linha chao
			glBegin(GL_LINES);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, -caixa_tamanho - posy[i], 0.0);
			glEnd();

			// Linha direita
			glBegin(GL_LINES);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(caixa_tamanho - posx[i], 0.0, 0.0);
			glEnd();

			// Linha esquerda
			glBegin(GL_LINES);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(-caixa_tamanho - posx[i], 0.0, 0.0);
			glEnd();

			// Linha frente
			glBegin(GL_LINES);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, caixa_tamanho - posz[i]);
			glEnd();

			// Linha fundo
			glBegin(GL_LINES);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, -caixa_tamanho - posz[i]);
			glEnd();
		}
// ===================== Circulo do perimetro ==============
//		if (show_perimetro) {
//			glColor3ub(255, 255, 0);
//			glutWireSphere(nucleo_proximity_free[i], 7, 7);
//			glutWireSphere(nucleo_proximity_free[i], 7, 7);
//		}
// ===================== Esfera Forca ======================
		if (show_power_sphere) {
			glColor3ub(255, 255, 0);
			glTranslatef(velocidade_x[i] * 5, velocidade_y[i] * 5, velocidade_z[i] * 5);
			glutSolidSphere(nucleo_proximity[i] * .35, 6, 6);
			glTranslatef(-velocidade_x[i] * 5, -velocidade_y[i] * 5, -velocidade_z[i] * 5);
		}
		// ===================== Esfera Electron ======================
		// Trying to use PKa, but actualy I think I shouldnt
		GLint electron_quantidade_pka = electron_quantidade[i];
		if (atomo_letra[i] == 'N' && i > 3) {
			electron_quantidade_pka = electron_quantidade[i] - 1;
		}

		for (GLint Es = 0; Es < electron_quantidade_pka; Es++) {
			if (electron_arested[i][Es][0] == -1) {
				glColor3ub(255, 0, 0);
			} else {
				glColor3ub(255, 200, 100);
			}
			if (electron_arested[i][Es][0] != -1) {
				if (show_base_line == false || (show_base_line && atomo_base[electron_arested[i][Es][0]] && atomo_base[i])) {
					glColor3ub(255, 255, 255);
					glBegin(GL_LINES);
					glVertex3f(0.0, 0.0, 0.0);
				}
				GLfloat difx = 0.0;
				GLfloat dify = 0.0;
				GLfloat difz = 0.0;

				GLfloat difxE = 0.0;
				GLfloat difyE = 0.0;
				GLfloat difzE = 0.0;
				GLfloat distancef = show_distance(posx[i], posx[electron_arested[i][Es][0]], posy[i], posy[electron_arested[i][Es][0]], posz[i], posz[electron_arested[i][Es][0]]);
				if (posx[i] > posx[electron_arested[i][Es][0]]) {
					difx = 0 - show_distance(posx[i], posx[electron_arested[i][Es][0]], 0.0, 0.0, 0.0, 0.0);
					difxE = electron_raio[i] * difx / distancef;
				} else {
					difx = show_distance(posx[i], posx[electron_arested[i][Es][0]], 0.0, 0.0, 0.0, 0.0);
					difxE = electron_raio[i] * difx / distancef;
				}
				if (posy[i] > posy[electron_arested[i][Es][0]]) {
					dify = 0 - show_distance(posy[i], posy[electron_arested[i][Es][0]], 0.0, 0.0, 0.0, 0.0);
					difyE = electron_raio[i] * dify / distancef;
				} else {
					dify = show_distance(posy[i], posy[electron_arested[i][Es][0]], 0.0, 0.0, 0.0, 0.0);
					difyE = electron_raio[i] * dify / distancef;
				}
				if (posz[i] > posz[electron_arested[i][Es][0]]) {
					difz = 0 - show_distance(posz[i], posz[electron_arested[i][Es][0]], 0.0, 0.0, 0.0, 0.0);
					difzE = electron_raio[i] * difz / distancef;
				} else {
					difz = show_distance(posz[i], posz[electron_arested[i][Es][0]], 0.0, 0.0, 0.0, 0.0);
					difzE = electron_raio[i] * difz / distancef;
				}
				if (show_base_line == false || (show_base_line && atomo_base[electron_arested[i][Es][0]] && atomo_base[i])) {
					glVertex3f(difx, dify, difz);
					glEnd();
				}
				if (difxE != difxE) {
					difxE = 0.0;
				}
				if (difyE != difyE) {
					difyE = 0.0;
				}
				if (difzE != difzE) {
					difzE = 0.0;
				}
				posEx[i][Es] = posx[i] + difxE;
				posEy[i][Es] = posy[i] + difyE;
				posEz[i][Es] = posz[i] + difzE;

			}
		}
// ===================================================
		glTranslatef(-posx[i], -posy[i], -posz[i]);

// ===================================================
	}
	// Adicionando proteina de comparacao
	if (show_ghost_protein) {
		// ---------------------------------------------
		glColor4f(1.0, 0.0, 0.0, 1.0);
		glBegin(GL_LINES);
		glVertex3f(5.0, 0.0, 0.0);
		glVertex3f(-5.0, 0.0, 0.0);
		glEnd();

		glBegin(GL_LINES);
		glColor4f(0.0, 1.0, 0.0, 1.0);
		glVertex3f(0.0, 5.0, 0.0);
		glVertex3f(0.0, -5.0, 0.0);
		glEnd();

		glColor4f(0.0, 0.0, 1.0, 1.0);
		glBegin(GL_LINES);
		glVertex3f(0.0, 0.0, 5.0);
		glVertex3f(0.0, 0.0, -5.0);
		glEnd();

		// Pegar posicao de proteina 0 e 2(equivalente) e colocar 0 na posicao 0,0,0
		// proteina
		GLfloat prot_0_x = 0;
		GLfloat prot_0_y = 0;
		GLfloat prot_0_z = 0;

		GLfloat prot_2_x = posx[ghost_protein_atom[2]] - posx[ghost_protein_atom[0]];
		GLfloat prot_2_y = posy[ghost_protein_atom[2]] - posy[ghost_protein_atom[0]];
		GLfloat prot_2_z = posz[ghost_protein_atom[2]] - posz[ghost_protein_atom[0]];

		GLfloat prot_2t_x = prot_2_x;
		GLfloat prot_2t_y = prot_2_y;
		GLfloat prot_2t_z = prot_2_z;

		// original
		glColor4f(0.0, 0.0, 1.0, 1.0);
		glTranslatef(0, 0, 0);
		glutWireSphere(0.4, 20, 20);
		glTranslatef(prot_2_x, prot_2_y, prot_2_z);
		glutWireSphere(0.5, 20, 20);
		glTranslatef(-prot_2_x, -prot_2_y, -prot_2_z);

		// ghost
		ghost_x[0] = ProteinOriginal[ghost_protein_atom[0]].x;
		ghost_y[0] = ProteinOriginal[ghost_protein_atom[0]].y;
		ghost_z[0] = ProteinOriginal[ghost_protein_atom[0]].z;
		ghost_x[2] = ProteinOriginal[ghost_protein_atom[2]].x;
		ghost_y[2] = ProteinOriginal[ghost_protein_atom[2]].y;
		ghost_z[2] = ProteinOriginal[ghost_protein_atom[2]].z;

		GLfloat ghost_0_x = 0;
		GLfloat ghost_0_y = 0;
		GLfloat ghost_0_z = 0;

		GLfloat ghost_2_x = ghost_x[2] - ghost_x[0];
		GLfloat ghost_2_y = ghost_y[2] - ghost_y[0];
		GLfloat ghost_2_z = ghost_z[2] - ghost_z[0];

		GLfloat ghost_2t_x = ghost_2_x;
		GLfloat ghost_2t_y = ghost_2_y;
		GLfloat ghost_2t_z = ghost_2_z;

//		GLfloat ghost_2_x_final = 0.0;
//		GLfloat ghost_2_y_final = 0.0;
//		GLfloat ghost_2_z_final = 0.0;

		glColor4f(1.0, 0.0, 0.0, 1.0);
		glTranslatef(0, 0, 0);
		glutWireSphere(0.5, 20, 20);
		glTranslatef(ghost_2_x, ghost_2_y, ghost_2_z);
		glutWireSphere(0.5, 20, 20);
		glTranslatef(-ghost_2_x, -ghost_2_y, -ghost_2_z);

		// ------------------- Tentando calcular o angulo ---------------
		// 0, 1, 0 eixo Y
		// 0, 0, 1 eixo Z
		// 1, 0, 0 eixo X

		// A proteina
		// B ghost

		// Achar angulo Y de A e rotacionar A para que fique perpendicular com X (protein_eixo_y)
		// Achar angulo Y de B e rotacionar B para que fique perpendicular com X (ghost_eixo_y)
		// Achar diferenca de angulo Y entre A e B (diferenca_eixo_y)
		// Encontrar angulo Z entre 0 e A, para isso utilizar X e Y
		// Encontrar angulo Z entre 0 e B, para isso utilizar X e Y
		// Achar diferenca de angulo Z entre A e B (diferenca_eixo_z)

		// Parte 2
		// Achar angulo Y de B e rotacionar B para que fique perpendicular com X (ghost_eixo_y)
		// Rotacionar no eixo Z usando variavel diferenca_eixo_z
		// Rotacionar de volta no eixo Y usando ghost_eixo_y + diferenca_eixo_y
// ---------------------------------------------------------------------

		// Achar angulo Y de A e rotacionar A para que fique perpendicular com X (protein_eixo_y)
		GLfloat protein_eixo_y = toAngle(atan(prot_2t_z / prot_2t_x));
		matrix3d(prot_2t_x, prot_2t_y, prot_2t_z, 0, 1, 0, protein_eixo_y);
		prot_2t_x = outputMatrix[0][0];
		prot_2t_y = outputMatrix[1][0];
		prot_2t_z = outputMatrix[2][0];

		// Achar angulo Y de B e rotacionar B para que fique perpendicular com X (ghost_eixo_y)
		GLfloat ghost_eixo_y = toAngle(atan(ghost_2t_z / ghost_2t_x));
		matrix3d(ghost_2t_x, ghost_2t_y, ghost_2t_z, 0, 1, 0, ghost_eixo_y);
		ghost_2t_x = outputMatrix[0][0];
		ghost_2t_y = outputMatrix[1][0];
		ghost_2t_z = outputMatrix[2][0];

//		glColor4f(1.0, 1.0, 0.0, 1.0);
//		glTranslatef(prot_2t_x, prot_2t_y, prot_2t_z);
//		glutWireSphere(0.5, 20, 20);
//		glTranslatef(-prot_2t_x, -prot_2t_y, -prot_2t_z);

		// Achar diferenca de angulo Y entre A e B (diferenca_eixo_y)
		GLfloat diferenca_eixo_y = ghost_eixo_y - protein_eixo_y;
		// Encontrar angulo Z entre 0 e A, para isso utilizar X e Y
		GLfloat protein_eixo_z = toAngle(atan(prot_2t_x / prot_2t_y));
		// Encontrar angulo Z entre 0 e B, para isso utilizar X e Y
		GLfloat ghost_eixo_z = toAngle(atan(ghost_2t_x / ghost_2t_y));
		// Achar diferenca de angulo Z entre A e B (diferenca_eixo_z)
		GLfloat diferenca_eixo_z = ghost_eixo_z - protein_eixo_z;


//		printf("Result: %f %f %f\n", protein_eixo_z, ghost_eixo_z, diferenca_eixo_z);

		// Parte 2
		// Achar angulo Y de B e rotacionar B para que fique perpendicular com X (ghost_eixo_y)
		ghost_eixo_y = toAngle(atan(ghost_2_z / ghost_2_x));
		matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 0, 1, 0, ghost_eixo_y);
		ghost_2_x = outputMatrix[0][0];
		ghost_2_y = outputMatrix[1][0];
		ghost_2_z = outputMatrix[2][0];

		// Rotacionar no eixo Z usando variavel diferenca_eixo_z
		matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 0, 0, 1, diferenca_eixo_z);
		ghost_2_x = outputMatrix[0][0];
		ghost_2_y = outputMatrix[1][0];
		ghost_2_z = outputMatrix[2][0];

		// Rotacionar de volta no eixo Y usando ghost_eixo_y + diferenca_eixo_y
//		matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 0, 1, 0, -(diferenca_eixo_y + ghost_eixo_y));
		matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 0, 1, 0, -(ghost_eixo_y-diferenca_eixo_y));
		ghost_2_x = outputMatrix[0][0];
		ghost_2_y = outputMatrix[1][0];
		ghost_2_z = outputMatrix[2][0];

//		printf("Result: %f %f %f\n", protein_eixo_y, ghost_eixo_y, diferenca_eixo_y);
		printf("Dist ghost: %f %f %f %f\n", show_distance(0, ghost_2_x, 0, ghost_2_y, 0, ghost_2_z), ghost_eixo_y, protein_eixo_y, -(ghost_eixo_y-diferenca_eixo_y));
		// ------ Angulo Z
//		GLfloat ghost_angulo_z = toAngle(atan(ghost_2_x / ghost_2_y));
//		GLfloat proteina_angulo_z = toAngle(atan(prot_2_x / prot_2_y));
//		GLfloat angulo_diferenca_z = ghost_angulo_z - proteina_angulo_z;
//		printf("Angulo Z %f %f %f\n", ghost_angulo, proteina_angulo, angulo_diferenca);
//		matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 0, 0, 1, angulo_diferenca);
//		ghost_2_x = outputMatrix[0][0];
//		ghost_2_y = outputMatrix[1][0];
//		ghost_2_z = outputMatrix[2][0];

		// ------ Angulo X
//		GLfloat ghost_angulo_x = toAngle(atan(ghost_2_z / ghost_2_y));
//		GLfloat proteina_angulo_x = toAngle(atan(prot_2_z / prot_2_y));
//		GLfloat angulo_diferenca_x = proteina_angulo_x - ghost_angulo_x;
//
//		matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 0, 0, 1, angulo_diferenca_z);
//		ghost_2_x = outputMatrix[0][0];
//		ghost_2_y = outputMatrix[1][0];
//		ghost_2_z = outputMatrix[2][0];

//		GLfloat trian_qual_a = show_distance(0, outputMatrix[0][0], 0, outputMatrix[1][0], 0, outputMatrix[2][0]);
//		GLfloat trian_qual_b = show_distance(prot_2_x, outputMatrix[0][0], prot_2_y, outputMatrix[1][0], prot_2_z, outputMatrix[2][0]);
//		GLfloat trian_qual_c = show_distance(0, prot_2_x, 0, prot_2_y, 0, prot_2_z);
//		GLfloat angulo_erro = toAngle(acos((pow(trian_qual_a, 2) + pow(trian_qual_c, 2) - pow(trian_qual_b, 2)) / (2 * trian_qual_a * trian_qual_c)));
//		GLfloat angulo_erro_proteina = toAngle(acos(prot_2_z / trian_qual_c));
//		GLfloat angulo_erro_ghost = toAngle(acos(prot_2_z / trian_qual_c));
//		GLfloat angulo_erro2 = acos(1.0);
//		GLfloat protein_vetor = show_distance(0, prot_2_x, 0, prot_2_y, 0, prot_2_z);
//		GLfloat protein_alfa = toAngle(acos(prot_2_x / protein_vetor));
//
//		GLfloat ghost_vetor = show_distance(0, ghost_2_x, 0, ghost_2_y, 0, ghost_2_z);
//		GLfloat ghost_alfa = toAngle(acos(ghost_2_x / ghost_vetor));
//		GLgloat ghost_vetor = show_distance(0, ghost_2_x, 0, ghost_2_y, 0, ghost_2_z);
//		printf("PV: %f %f %f %f\n", proteina_angulo_x, protein_alfa, ghost_angulo_x, ghost_alfa);
//		printf("Diferenca: %f %f %f %f %f %f\n", ghost_angulo_z, proteina_angulo_z, angulo_diferenca_z, prot_2_x, prot_2_y, angulo_erro);
//		printf("Angulo %f %f %f\n", angulo_diferenca_z, angulo_diferenca_x, angulo_erro);
//		printf("Erro de: %f\n");

//		matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 0, 1, 0, proteina_angulo_x);
//		ghost_2_x = outputMatrix[0][0];
//		ghost_2_y = outputMatrix[1][0];
//		ghost_2_z = outputMatrix[2][0];

		glColor4f(1.0, 1.0, 0.0, 1.0);
		glBegin(GL_LINES);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(prot_2_x, prot_2_y, prot_2_z);
		glEnd();

		glColor4f(0.0, 1.0, 0.0, 1.0);
		glTranslatef(ghost_2_x, ghost_2_y, ghost_2_z);
		glutWireSphere(0.5, 20, 20);
		glTranslatef(-ghost_2_x, -ghost_2_y, -ghost_2_z);
		// ------------------------------- Angulo X -----------------------------------
//		GLfloat distancia = 0.0;
//		GLfloat ultima_distancia = 100.0;
//		GLfloat ghost_angulo_z = 0.0;
//
//		for (GLint lo = 0; lo < 3; lo++) {
//			distancia = 0.0;
//			ultima_distancia = 100.0;
//			ghost_angulo_z = 0.0;
//			for (GLfloat di = 0.0; di < 360; di++) {
//				matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 0, 0, 1, di);
//				distancia = show_distance(prot_2_x, outputMatrix[0][0], prot_2_y, outputMatrix[1][0], prot_2_z, outputMatrix[2][0]);
//				if (distancia < ultima_distancia) {
//					ghost_angulo_z = di;
//				}
//				ultima_distancia = distancia;
//			}
//			matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 0, 0, 1, ghost_angulo_z);
//			ghost_2_x = outputMatrix[0][0];
//			ghost_2_y = outputMatrix[1][0];
//			ghost_2_z = outputMatrix[2][0];
//
//			// ------------------------------- Angulo Z -----------------------------------
//			distancia = 0.0;
//			ultima_distancia = 100.0;
//			GLfloat ghost_angulo_x = 0.0;
//			for (GLfloat di = 0.0; di < 360; di++) {
//				matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 0, 1, 0, di);
//				distancia = show_distance(prot_2_x, outputMatrix[0][0], prot_2_y, outputMatrix[1][0], prot_2_z, outputMatrix[2][0]);
//				if (distancia < ultima_distancia) {
//					ghost_angulo_x = di;
//				}
//				ultima_distancia = distancia;
//			}
//			matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 0, 1, 0, ghost_angulo_x);
//			ghost_2_x = outputMatrix[0][0];
//			ghost_2_y = outputMatrix[1][0];
//			ghost_2_z = outputMatrix[2][0];
//
//			// ------------------------------- Angulo Z -----------------------------------
//			distancia = 0.0;
//			ultima_distancia = 100.0;
//			GLfloat ghost_angulo_y = 0.0;
//			for (GLfloat di = 0.0; di < 360; di++) {
//				matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 1, 0, 0, di);
//				distancia = show_distance(prot_2_x, outputMatrix[0][0], prot_2_y, outputMatrix[1][0], prot_2_z, outputMatrix[2][0]);
//				if (distancia < ultima_distancia) {
//					ghost_angulo_y = di;
//				}
//				ultima_distancia = distancia;
//			}
//			matrix3d(ghost_2_x, ghost_2_y, ghost_2_z, 1, 0, 0, ghost_angulo_y);
//			ghost_2_x = outputMatrix[0][0];
//			ghost_2_y = outputMatrix[1][0];
//			ghost_2_z = outputMatrix[2][0];
//		}
//
//		distancia = show_distance(prot_2_x, ghost_2_x, prot_2_y, ghost_2_y, prot_2_z, ghost_2_z);
//		printf("Distancia final %f\n", distancia);
//		glColor4f(0.0, 1.0, 0.0, 1.0);
//		glTranslatef(ghost_2_x, ghost_2_y, ghost_2_z);
//		glutWireSphere(0.5, 20, 20);
//		glTranslatef(-ghost_2_x, -ghost_2_y, -ghost_2_z);

//		}
		// Mover X ateh q distancia seje menor
		// Mover Y ateh q distancia seje menor
		// Mover Z ateh q distancia seje menor

		// ---------------------------------------------
		// Funcionando -------------------------------------------
//		GLfloat diferenca_x = posx[ghost_protein_atom[0]] - ProteinOriginal[ghost_protein_atom[0]].x;
//		GLfloat diferenca_y = posy[ghost_protein_atom[0]] - ProteinOriginal[ghost_protein_atom[0]].y;
//		GLfloat diferenca_z = posz[ghost_protein_atom[0]] - ProteinOriginal[ghost_protein_atom[0]].z;
//
//		GLfloat ox1 = ProteinOriginal[ghost_protein_atom[0]].x;
//		GLfloat oy1 = ProteinOriginal[ghost_protein_atom[0]].y;
//		GLfloat ox2 = ProteinOriginal[ghost_protein_atom[2]].x;
//		GLfloat oy2 = ProteinOriginal[ghost_protein_atom[2]].y;
//		GLfloat oh = show_distance(ox1, ox2, oy1, oy2, 0, 0);
//		GLfloat oax = abs(show_distance(ox1, ox2, 0, 0, 0, 0));
//		GLfloat oay = abs(show_distance(oy1, oy2, 0, 0, 0, 0));
//		GLfloat oa1 = acos(oax / oh);
//
//		GLfloat fx1 = posx[ghost_protein_atom[0]];
//		GLfloat fy1 = posy[ghost_protein_atom[0]];
//		GLfloat fx2 = posx[ghost_protein_atom[2]];
//		GLfloat fy2 = posy[ghost_protein_atom[2]];
//		GLfloat fh = show_distance(fx1, fx2, fy1, fy2, 0, 0);
//		GLfloat fax = abs(show_distance(fx1, fx2, 0, 0, 0, 0));
//		GLfloat fay = abs(show_distance(fy1, fy2, 0, 0, 0, 0));
//		GLfloat fa1 = acos(fax / fh);
//
//		GLfloat difOF = fa1 - oa1;
//
//		GLfloat oz1, oz2, oaz, fz1, fz2, faz, difOFz;
//
//
//		for (int i = 0; i < ProteinOriginal.size(); i++) {
//			if (ProteinOriginal[i].show) {
//				if (i == 3 || i == 99) {
//					matrix3d(ProteinOriginal[i].x - ProteinOriginal[3].x, ProteinOriginal[i].y - ProteinOriginal[3].y, ProteinOriginal[i].z - ProteinOriginal[3].z, 0, 0, 1, toAngle(difOF));
//					ox1 = (outputMatrix[0][0] + diferenca_x + ProteinOriginal[3].x);
//					oz1 = ProteinOriginal[ghost_protein_atom[0]].z;
//					ox2 = ProteinOriginal[ghost_protein_atom[2]].x;
//					oz2 = ProteinOriginal[ghost_protein_atom[2]].z;
//					oh = show_distance(ox1, ox2, oz1, oz2, 0, 0);
//					oax = abs(show_distance(ox1, ox2, 0, 0, 0, 0));
//					oaz = abs(show_distance(oz1, oz2, 0, 0, 0, 0));
//					oa1 = acos(oax / oh);
//
//					fx1 = posx[ghost_protein_atom[0]];
//					fz1 = posz[ghost_protein_atom[0]];
//					fx2 = posx[ghost_protein_atom[2]];
//					fz2 = posz[ghost_protein_atom[2]];
//					fh = show_distance(fx1, fx2, fz1, fz2, 0, 0);
//					fax = abs(show_distance(fx1, fx2, 0, 0, 0, 0));
//					faz = abs(show_distance(fz1, fz2, 0, 0, 0, 0));
//					fa1 = acos(fax / fh);
//					difOFz = fa1 - oa1;
//					glColor4f(1.0, 0.0, 0.0, 1.0);
//					glTranslatef((outputMatrix[0][0] + diferenca_x + ProteinOriginal[3].x), (outputMatrix[1][0] + diferenca_y + ProteinOriginal[3].y), (outputMatrix[2][0] + diferenca_z + ProteinOriginal[3].z));
//					glutWireSphere(0.5, 20, 20);
//					glTranslatef(-(outputMatrix[0][0] + diferenca_x + ProteinOriginal[3].x), -(outputMatrix[1][0] + diferenca_y + ProteinOriginal[3].y), -(outputMatrix[2][0] + diferenca_z + ProteinOriginal[3].z));
//				}
//			}
//		}
	}

	if (cont_loop_electron == cont_loop_electron_time) {
		cont_loop_electron = 0;
	}
	cont_loop_electron++;
// Monta percurso
	if (massa[0] != 0) {
		for (GLint is = 0; is <= percurso_contador; is++) {
			glColor3ub(255, 255, 255);
			glTranslatef(percurso[is][0], percurso[is][1], percurso[is][2]);
			glutSolidSphere(0.1f, 6, 6);
			glTranslatef(-percurso[is][0], -percurso[is][1], -percurso[is][2]);
		}
// Adiciona posicao atual ao percurso
		if (rastreio) {
			if (percurso_contador_loop == 10 && percurso_contador < 498) {
				percurso_contador++;
				percurso[percurso_contador][0] = posx[0];
				percurso[percurso_contador][1] = posy[0];
				percurso[percurso_contador][2] = posz[0];
				percurso_contador_loop = 0;
			}
			percurso_contador_loop++;
		}
	}
	forca_externa_contador++;
	forca_externa_contador_t++;
	forca_externa_contador_hb++;
	forca_externa_contador_db++;
	if (forca_externa_contador > forca_externa_contador_max) {
		forca_externa_contador = 0;
	}
	if (forca_externa_contador_t > forca_externa_contador_max_t) {
		forca_externa_contador_t = 0;
	}
	if (forca_externa_contador_hb > forca_externa_contador_max_hb) {
		forca_externa_contador_hb = 0;
	}
	if (forca_externa_contador_db > forca_externa_contador_max_db) {
		forca_externa_contador_db = 0;
	}
// -----------------------------------------------------------
	contadorFrames++;
	double currentTime = glfwGetTime();
	if (currentTime - lastTime >= 1.0) {
		stringstream ss;
		ss << (1000.0 / double(contadorFrames));
		fps = ss.str();
		contadorFrames = 0;
		lastTime += 1.0;
	}

// Imprime texto na tela
	parameters_window();
//	glMatrixMode(GL_PROJECTION);
//	glPushMatrix();
//	glLoadIdentity();
//	gluOrtho2D(0.0, windowWidth, 0.0, windowHeight);
//	glMatrixMode(GL_MODELVIEW);
//	glPushMatrix();
//	glLoadIdentity();
//	glColor3f(0.0, 1.0, 0.0);				// Green
//	glRasterPos2i(10, 10);
//	void * font = GLUT_BITMAP_9_BY_15;
//	for (string::iterator i = fps.begin(); i != fps.end(); ++i) {
//		char c = *i;
//		glutBitmapCharacter(font, c);
//	}
//	glMatrixMode(GL_MODELVIEW);
//	glPopMatrix();
//
//	glMatrixMode(GL_PROJECTION);
//	glPopMatrix();
//
//	parameters_window();
// ----- Stop Drawing Stuff! ------
	glfwSwapBuffers();				// Swap the buffers to display the scene (so we don't have to watch it being drawn!)
}

static void mouse_button_callback(int key, int action) {
//	printf("Evento mouse %d %d\n", key, action);
	if (key == 0 && action == 0) {
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);
		GLdouble objX;
		GLdouble objY;
		GLdouble objZ;
		GLdouble objXf;
		GLdouble objYf;
		GLdouble objZf;
		GLdouble modelview[16];
		glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
		GLdouble projection[16];
		glGetDoublev(GL_PROJECTION_MATRIX, projection);
		gluUnProject(viewport[2] / 2, viewport[3] - viewport[3] / 2, 0, modelview, projection, viewport, &objX, &objY, &objZ);
		gluUnProject(viewport[2] / 2, viewport[3] - viewport[3] / 2, 1, modelview, projection, viewport, &objXf, &objYf, &objZf);
		GLint encontrado_i = -1;
		GLdouble encontrado_d = 100;
		GLdouble res_calculo = 0;
		for (GLint i = 0; i < atomos_quantidade; i++) {
			res_calculo = fabs(show_distance(objX, objXf, objY, objYf, objZ, objZf) - show_distance(posx[i], objX, posy[i], objY, posz[i], objZ) - show_distance(posx[i], objXf, posy[i], objYf, posz[i], objZf));
			if (res_calculo < encontrado_d) {
				encontrado_d = res_calculo;
				encontrado_i = i;
			}
		}
		if (encontrado_d <= 0.003) {
			vp["global"]["atomo_selecionado"][""][""][""][""] = encontrado_i;
		} else {
			vp["global"]["atomo_selecionado"][""][""][""][""] = -1;
		}
		printf("Atomo selecionado: %d\n", encontrado_i);
		printf("%s A:%d %d V(%f %f %f)(%f %f %f)\n", atomo_letraL[encontrado_i].c_str(), amino[encontrado_i], amino_sequencial[encontrado_i], velocidade_x[encontrado_i], velocidade_y[encontrado_i], velocidade_z[encontrado_i], velocidade_x_backup[encontrado_i], velocidade_y_backup[encontrado_i], velocidade_z_backup[encontrado_i]);
	}
}

// Fire it up...
int main(int argc, char **argv) {
	cout << "Controls: Use WSAD and the mouse to move around!" << endl;

// Frame counter and window settings variables
	int redBits = 8, greenBits = 8, blueBits = 8;
	int alphaBits = 8, depthBits = 24, stencilBits = 0;

// Flag to keep our main loop running
	bool running = true;

// ----- Intialiase FreeGLUT -----

// Note: We're only using freeGLUT to draw some spheres, so if you modify the code to not include any calls
// to glutSolidSphere, then you don't need this, the header or the lib...
	glutInit(&argc, argv);

// Initialise GLFW
	glfwInit();

// Ask for 4x AntiAliasing (this doesn't mean we'll get it - it'll work only if the GLX_ARB_multisample extension is available)
// Note: Hints must be provided BEFORE the window is opened! But we can't query for it with GLEE until the window is opened! Catch 22!
	glfwOpenWindowHint(GLFW_FSAA_SAMPLES, 4);

// Create a window
	if (!glfwOpenWindow(windowWidth, windowHeight, redBits, greenBits, blueBits, alphaBits, depthBits, stencilBits, GLFW_WINDOW)) {
		cout << "Failed to open window!" << endl;
		glfwTerminate();
		return 0;
	}

// ----- Initialise GLEE -----

// Initialise GLee once we've got a rendering context
// Note: We don't really have to do this because it's called automatically, but if we do it - we KNOW it's been called!
//	GLeeInit();
//
//// Check for the OpenGL extension which allows us to use vsync
//	if (GLEE_GLX_SGI_swap_control) {
//		cout << "Extension found: GLX_SGI_swap_control (vsync can be used)." << endl;
//		glfwSwapInterval(1);
//	} else {
//		cout << "Extension NOT found: GLX_SGI_swap_control (vsync cannot be used)." << endl;
//		glfwSwapInterval(0);
//	}
//
//// Check for the OpenGL extension which allows us to use antialiasing
//	if (GLEE_ARB_multitexture) {
//		cout << "Extension found: GLX_ARB_multitexture (anti-aliasing can be used)." << endl;
//
//// If the extension's available, we likely got anti-aliasing, so disable line smoothing as it comes free with the AA
//		glDisable(GL_LINE_SMOOTH);
//	} else {
//		cout << "Extension NOT found: GLX_ARB_multitexture (anti-aliasing cannot be used)." << endl;
//
//// If the extention's not available, turn on line smoothing
//		glEnable(GL_LINE_SMOOTH);
//	}

// Set the mouse cursor to the centre of our window
	glfwSetMousePos(midWindowX, midWindowY);

// Call our initGL function to set up our OpenGL options
	initGL();

// Specify the function which should execute when a key is pressed or released
	glfwSetKeyCallback(handleKeypress);

// Specify the function which should execute when the mouse is moved
	glfwSetMousePosCallback(handleMouseMove);
	glfwSetMouseButtonCallback(mouse_button_callback);

	while (running == true) {
		calculaposicoes();

// Draw our scene
		drawScene();

// Calculate our camera movement
		calculateCameraMovement();

// Move our camera
		moveCamera();

// Increase our frame counter
		frameCount++;

// Check for any OpenGL errors (providing the location we called the function from)
		checkGLError("Main loop");

// exit if ESC was pressed or window was closed
		running = !glfwGetKey(GLFW_KEY_ESC) && glfwGetWindowParam(GLFW_OPENED);
	}

// Clean up GLFW and exit
	glfwTerminate();

	return 0;
}
