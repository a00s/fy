#include <iostream>
#include <cmath>

using namespace std;

typedef struct {
    GLfloat x;
    GLfloat y;
    GLfloat z;
}Point;
Point points;

GLfloat rotationMatrix[4][4];
GLfloat inputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};
GLfloat outputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};

void showPoint(){
    cout<<"("<<outputMatrix[0][0]<<","<<outputMatrix[1][0]<<","<<outputMatrix[2][0]<<")"<<endl;
}

void multiplyMatrix()
{
    for(int i = 0; i < 4; i++ ){
        for(int j = 0; j < 1; j++){
            outputMatrix[i][j] = 0;
            for(int k = 0; k < 4; k++){
                outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
            }
        }
    }
}
void setUpRotationMatrix(GLfloat angle, GLfloat u, GLfloat v, GLfloat w)
{
    GLfloat L = (u*u + v * v + w * w);
    angle = angle * M_PI / 180.0; //converting to radian value
    GLfloat u2 = u * u;
    GLfloat v2 = v * v;
    GLfloat w2 = w * w;

    rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
    rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
    rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
    rotationMatrix[0][3] = 0.0;

    rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
    rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
    rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
    rotationMatrix[1][3] = 0.0;

    rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
    rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
    rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
    rotationMatrix[2][3] = 0.0;

    rotationMatrix[3][0] = 0.0;
    rotationMatrix[3][1] = 0.0;
    rotationMatrix[3][2] = 0.0;
    rotationMatrix[3][3] = 1.0;
}

void matrix3d(GLfloat x_local, GLfloat y_local, GLfloat z_local, GLfloat u_axis, GLfloat v_axis, GLfloat w_axis, GLfloat angle_local)
{
//    GLfloat angle;
//    GLfloat u, v, w;
//    cout<<"Enter the initial point you want to transform:";
//    cin>>points.x>>points.y>>points.z;
    inputMatrix[0][0] = x_local;
    inputMatrix[1][0] = y_local;
    inputMatrix[2][0] = z_local;
    inputMatrix[3][0] = 1.0;

//    cout<<"Enter axis vector: ";
//    cin>>u>>v>>w;

//    cout<<"Enter the rotating angle in degree: ";
//    cin>>angle;

    setUpRotationMatrix(angle_local, u_axis, v_axis, w_axis);
    multiplyMatrix();
//    showPoint();

//    return 0;
}
