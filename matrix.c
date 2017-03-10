// Modified by Alex Walker and Jack Wines, 2/11/17

/*** 2 x 2 Matrices ***/
#include <stdio.h>

/* Pretty-prints the given matrix, with one line of text per row of matrix. */
void mat22Print(double m[2][2]) {
    int iMax = 2;
    int jMax = 2;
	for(int i = 0; i < iMax; i++)  {
        printf("|%f %f|\n", m[i][0], m[i][1]);
    }
    printf("\n");
}

/* Returns the determinant of the matrix m. If the determinant is 0.0, then the 
matrix is not invertible, and mInv is untouched. If the determinant is not 0.0, 
then the matrix is invertible, and its inverse is placed into mInv. */
double mat22Invert(double m[2][2], double mInv[2][2]) {
    double determinant = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    if(determinant == 0.0) {
        return determinant;
    }
    mInv[0][0] = m[1][1] / determinant;
    mInv[0][1] = m[0][1] / determinant * (-1);
    mInv[1][0] = m[1][0] / determinant * (-1);
    mInv[1][1] = m[0][0] / determinant;
    return determinant;
}

/* Multiplies a 2x2 matrix m by a 2-column v, storing the result in mTimesV. 
The output should not */
void mat221Multiply(double m[2][2], double v[2], double mTimesV[2]) {
    int matSize = 2;
    for(int i = 0; i < matSize; i++)  {
        double result = 0.0;
        for(int j = 0; j < matSize; j++)  {
            result += m[i][j]*v[j];
        }   
        mTimesV[i] = result;
    }
}

/* Fills the matrix m from its two columns. */
void mat22Columns(double col0[2], double col1[2], double m[2][2]) {
	m[0][0] = col0[0];
    m[0][1] = col1[0];
    m[1][0] = col0[1];
    m[1][1] = col1[1];
}

/* Pretty-prints the given matrix, with one line of text per row of matrix. */
void mat33Print(double m[3][3]) {
    int iMax = 3;
    int jMax = 3;
	for(int i = 0; i < iMax; i++)  {
        printf("|%f %f %f|\n", m[i][0], m[i][1], m[i][2]);
    }
    printf("\n");
}

/* Multiplies the 3x3 matrix m by the 3x3 matrix n. */
void mat333Multiply(double m[3][3], double n[3][3], double mTimesN[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mTimesN[i][j] = (m[i][0] * n[0][j]) + (m[i][1] * n[1][j]) + (m[i][2] * n[2][j]);
        }
    }
}

/* Multiplies the 3x3 matrix m by the 3x1 matrix v. */
void mat331Multiply(double m[3][3], double v[3], double mTimesV[3]) {
    for (int i = 0; i < 3; i++) {
        double result = 0.0;
        for (int j = 0; j < 3; j++) {
            result = result + m[i][j] * v[j];
        }
        mTimesV[i] = result;
    }
}

/* Builds a 3x3 matrix representing 2D rotation and translation in homogeneous 
coordinates. More precisely, the transformation first rotates through the angle 
theta (in radians, counterclockwise), and then translates by the vector (x, y). 
*/
void mat33Isometry(double theta, double x, double y, double isom[3][3]) {
    double rotationMatrix[3][3] = {{cos(theta), -sin(theta), x},
                                   {sin(theta),  cos(theta), y},
                                   {         0,           0, 1}};
    for (int i = 0; i < 3; ++i)  {
        for (int j = 0; j < 3; ++j)  {
            isom[i][j] = rotationMatrix[i][j];
        }
    }
}

// Adds two 3x3 matrices u and v, and stores the result in uPlusV
void mat33Add(double u[3][3], double v[3][3], double uPlusV[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            uPlusV[i][j] = u[i][j] + v[i][j];
        }
    }
}

// Transposes 3x3 matrix m and stores it in transM
void mat33Transpose(double m[3][3], double transM[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            transM[i][j] = m[j][i];
        }
    }
}

/* Given a length-1 3D vector axis and an angle theta (in radians), builds the 
rotation matrix for the rotation about that axis through that angle. Based on 
Rodrigues' rotation formula R = I + (sin theta) U + (1 - cos theta) U^2. */
void mat33AngleAxisRotation(double theta, double axis[3], double rot[3][3]) {
    double u[3][3] = {{       0, -axis[2],  axis[1]},
                      { axis[2],        0, -axis[0]},
                      {-axis[1],  axis[0],        0}};
    
    double uSquared[3][3];
    mat333Multiply(u, u, uSquared);
    
    double identity[3][3] = {{1, 0, 0},
                             {0, 1, 0},
                             {0, 0, 1}};
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rot[i][j] = identity[i][j] + sin(theta) * u[i][j] + (1 - cos(theta)) * uSquared[i][j];
        }
    }
}

/* Given two length-1 3D vectors u, v that are perpendicular to each other. 
Given two length-1 3D vectors a, b that are perpendicular to each other. Builds 
the rotation matrix that rotates u to a and v to b. */
void mat33BasisRotation(double u[3], double v[3], double a[3], double b[3], 
        double rot[3][3]) {
    double w[3];
    vec3Cross(u, v, w);
    double r[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            r[i][0] = u[i];
            r[i][1] = v[i];
            r[i][2] = w[i];
        }
    }
    
    double c[3];
    vec3Cross(a, b, c);
    double s[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            s[i][0] = a[i];
            s[i][1] = b[i];
            s[i][2] = c[i];
        }
    }
    
    double transR[3][3];
    mat33Transpose(r, transR);
    
    mat333Multiply(s, transR, rot);
}

/* Pretty-prints the given matrix, with one line of text per row of matrix. */
void mat44Print(double m[4][4]) {
    int iMax = 4;
    int jMax = 4;
    for(int i = 0; i < iMax; i++)  {
        printf("|%f %f %f %f|\n", m[i][0], m[i][1], m[i][2], m[i][3]);
    }
    printf("\n");
}

/* Multiplies m by n, placing the answer in mTimesN. */
void mat444Multiply(double m[4][4], double n[4][4], double mTimesN[4][4]) {
//    for (int i = 0; i < 4; i++) {
//        for (int j = 0; j < 4; j++) {
//            double result = 0.0;
//            for (int k = 0; k < 4; k++) {
//                result = result + m[i][k] * n[k][j];
//            }
//            mTimesN[i][j] = result;
//        }
//    }
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mTimesN[i][j] = (m[i][0] * n[0][j]) + (m[i][1] * n[1][j]) + (m[i][2] * n[2][j]) + (m[i][3] * n[3][j]);
        }
    }
}

/* Multiplies m by v, placing the answer in mTimesV. */
void mat441Multiply(double m[4][4], double v[4], double mTimesV[4]) {
    for (int i = 0; i < 4; i++) {
        double result = 0.0;
        for (int j = 0; j < 4; j++) {
            result += m[i][j] * v[j];
        }
        mTimesV[i] = result;
    }
}

/* Given a rotation and a translation, forms the 4x4 homogeneous matrix 
representing the rotation followed in time by the translation. */
void mat44Isometry(double rot[3][3], double trans[3], double isom[4][4]) {
    double temp[4][4] = {{rot[0][0], rot[0][1], rot[0][2], trans[0]},
                         {rot[1][0], rot[1][1], rot[1][2], trans[1]},
                         {rot[2][0], rot[2][1], rot[2][2], trans[2]},
                         {        0,         0,         0,        1}};
    
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            isom[i][j] = temp[i][j];
        }
    }
}

/* Given a rotation and translation, forms the 4x4 homogeneous matrix 
representing the inverse translation followed in time by the inverse rotation. 
That is, the isom produced by this function is the inverse to the isom 
produced by mat44Isometry on the same inputs. */
void mat44InverseIsometry(double rot[3][3], double trans[3], 
        double isom[4][4]) {
    double rotInverse[3][3];
    mat33Transpose(rot, rotInverse);
    mat44Isometry(rotInverse, trans, isom);
    
    double transInverse[3] = {-trans[0], -trans[1], -trans[2]};
    double mtInverse[3];
    mat331Multiply(rotInverse, transInverse, mtInverse);
    double temp[4][4] = {{rot[0][0], rot[1][0], rot[2][0], mtInverse[0]},
                         {rot[0][1], rot[1][1], rot[2][1], mtInverse[1]},
                         {rot[0][2], rot[1][2], rot[2][2], mtInverse[2]},
                         {        0,         0,         0,            1}};
    
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            isom[i][j] = temp[i][j];
        }
    }
}

/* Builds a 4x4 matrix representing orthographic projection with a boxy viewing 
volume [left, right] x [bottom, top] x [far, near]. That is, on the near plane 
the box is the rectangle R = [left, right] x [bottom, top], and on the far 
plane the box is the same rectangle R. Keep in mind that 0 > near > far. Maps 
the viewing volume to [-1, 1] x [-1, 1] x [-1, 1]. */
void mat44Orthographic(double left, double right, double bottom, double top, 
        double far, double near, double proj[4][4]) {
    double temp[4][4] = {{2 / (right - left), 0, 0, (-right - left) / (right - left)},
                         {0, 2 / (top - bottom), 0, (-top - bottom) / (top - bottom)},
                         {0, 0, -2 / (near - far), (near + far) / (near - far)},
                         {0, 0, 0, 1}};
    
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            proj[i][j] = temp[i][j];
        }
    }
}

/* Builds a 4x4 matrix that maps a projected viewing volume 
[-1, 1] x [-1, 1] x [-1, 1] to screen [0, w - 1] x [0, h - 1] x [-1, 1]. */
void mat44Viewport(double width, double height, double view[4][4]) {
    double temp[4][4] = {{(width - 1) / 2, 0, 0, (width - 1) / 2},
                         {0, (height - 1) / 2, 0, (height - 1) / 2},
                         {0, 0, 1, 0},
                         {0, 0, 0, 1}};
    
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            view[i][j] = temp[i][j];
        }
    }
}

void mat33Identity(double m[3][3])  {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            m[i][j] = i == j;
        }
    }
}

void mat44Identity(double m[4][4])  {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            m[i][j] = i == j;
        }
    }
}

/* Builds a 4x4 matrix representing perspective projection. The viewing frustum 
is contained between the near and far planes, with 0 > near > far. On the near 
plane, the frustum is the rectangle R = [left, right] x [bottom, top]. On the 
far plane, the frustum is the rectangle (far / near) * R. Maps the viewing 
volume to [-1, 1] x [-1, 1] x [-1, 1]. */
void mat44Perspective(double left, double right, double bottom, double top, 
        double far, double near, double proj[4][4]) {
    double temp[4][4] = {{(-2 * near) / (right - left), 0, (right + left) / (right - left), 0},
                         {0, (-2 * near) / (top - bottom), (top + bottom) / (top - bottom), 0},
                         {0, 0, -(-near - far) / (near - far), -(2 * near * far) / (near - far)},
                         {0, 0, -1, 0}};
    
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            proj[i][j] = temp[i][j];
        }
    }
}