extern GLint windowWidth;
extern GLint windowHeight;
extern GLint calibration_precision;

extern int atomos_quantidade;

extern bool holdingForward;
extern bool holdingBackward;
extern bool holdingLeftStrafe;
extern bool holdingRightStrafe;
extern bool pressionando_k;
extern bool pressionando_j;
extern bool pressionando_i;
extern bool pressionando_m;
extern bool pressionando_o;
extern bool pressionando_p;
extern bool pressionando_0;
extern bool pressionando_l;
extern bool pressionando_8;
extern bool pressionando_9;
extern bool pressionando_control;
extern bool show_power_sphere;
extern bool show_perimetro;
extern bool show_base_line;
extern bool show_base;
extern bool rastreio;
extern bool paused;
extern bool show_colisao_tensao;
extern bool show_tensao_hb;
extern bool show_comparation;
extern bool show_comparation_speed;


extern GLfloat posx[1000];
extern GLfloat posy[1000];
extern GLfloat posz[1000];
extern GLfloat posxGhost[1000];
extern GLfloat posyGhost[1000];
extern GLfloat poszGhost[1000];

extern GLfloat caixa_tamanho;
extern GLint ghost_angulo;

extern GLfloat show_distance(GLfloat x1, GLfloat x2, GLfloat y1, GLfloat y2, GLfloat z1, GLfloat z2);
extern void add_atom(char atomo_ll);
extern void camera_position(GLfloat px, GLfloat py, GLfloat pz, GLfloat rx, GLfloat ry, GLfloat rz);
extern void compare_protein_build_MD(int c_distancia);

extern void ativa_desativa_colisao_tensao();
extern void ativa_desativa_tensao_hb();
extern void ativa_desativa_comparation();
extern void restaura_posicoes();
extern void ativa_desativa_comparation_speed();
extern void distance_calibration();
extern void show_variables();
extern void add_chain();
extern void ativa_desativa_perimetro();
extern void ativa_desativa_rastreio();
extern void ativa_desativa_forca();
extern void rem_energy();
extern void add_energy();
extern void continua();
extern void ativa_desativa_base();
extern void ativa_desativa_base_line();
extern void ativa_desativa_ghost_protein();

extern GLfloat angulo_adicional_teste;
extern GLint sequencial_mostra;
