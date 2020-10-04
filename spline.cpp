#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/gl.h>
#include <windows.h>
#include "spline.h"

//эта функция определяет сплайн по заданным параметрам (количество точек, начальное время, струтуры (описаны выше))
void* spline(int num, double t, data_input *dataIn)
{
    FILE *f_phi;
    FILE *f_lambda;
    FILE *f_height;
    FILE *f_accel;
    FILE *f_angul;
    double omega = 7.292115e-5;//угловая скорость вращения Земли
    double a0 = 6378245.0;// экваториальный радиус Земли
    double e = 0.08181336987204;//эксантриситет Земного эллипсоида
    int g = 0;
    //коэффициенты сплайна
    double A, dir_A, dir2_A, B, dir_B, dir2_B;
    //сплайн и его производные для широты
    double phi_S, phi_dirS, phi_dir2S;
    //сплайн и его производные для долготы
    double lambda_S, lambda_dirS, lambda_dir2S;
    //сплайн и его производные для высоты
    double h_S, h_dirS, h_dir2S;
    int i;//счетчик цикла
    double t_old;//время старое
    double t_new;//время новое
    double hag;//сглаживание
    double teta, dir_teta;//тангаж
    double gamma, dir_gamma;//крен
    //высота
    double H, dir_H, dir2_H;
    //величина траекторной скорости
    double W, dir_W, dir2_W;
    //истинный курс
    double Psi, dir_Psi, dir2_Psi;
    //аппроксимированные значения широты и долготы
    double phi, lambda;
    //относительное ускорение и скорость геодезической системы
    //линейные составляющие
    double dir_We1, dir_Wn1, dir_Wh1, We, Wn, Wh;
    //вращательные составляющие
    double w_wave_E, w_wave_N, w_wave_H;
    double dir_We2, dir_Wn2, dir_Wh2;
    //суммарные составляющие
    double dir_Ve, dir_Vn, dir_Vh, Ve, Vn, Vh;
    //проекции векторов относительного ускорения и скорости на оси опорного трехгранника
    double v_ksi, v_eta, v_zeta, dir_v_ksi, dir_v_eta, dir_v_zeta;
    //проекции вектора относительной угловой скорости вращения опорного трехгранника
    double w_ksi, w_eta, w_zeta;
    //проекции вектора угловой скорости вращения Земли на оси опорного трехгранника
    double omega_ksi, omega_eta, omega_zeta;
    //проекции вектора угловой скорости вращения объекта вокруг центра масс на оси опорного треугольника
    double w_wave_ksi, w_wave_eta, w_wave_zeta;
    //проекции вектора абсолютной угловой скорости вращения объекта на оси опорного трехгранника
    double w_abs_ksi, w_abs_eta, w_abs_zeta;
    //проекции вектора консервативных сил на оси опорного трехгранника
    double b_ksi, b_eta, b_zeta;
    //проекции вектора ускорения от активных сил на оси опорного трехгранника
    double a_ksi, a_eta, a_zeta;
    double MNK[3][3];//МНК
    //данные с акселерометра
    double a_x, a_y, a_z;
    //данные c датчиков угловой скорости ДУС
    double w_x, w_y, w_z;
    printf("Enter smoothing step:");
    scanf("%lf", &hag);
    if(num > 1)
    {
        f_phi = fopen("D:/Praktika/spline/data_Phi", "w+t");
        f_lambda = fopen("D:/Praktika/spline/data_Lambda", "w+t");
        f_height = fopen("D:/Praktika/spline/data_Height", "w+t");
        f_accel = fopen("D:/Praktika/spline/data_Accel", "w+t");
        f_angul = fopen("D:/Praktika/spline/data_Angul", "w+t");
        fprintf(f_phi,"----S--------S'--------S''--------t--------t_i--------t_i+1---|\n");
        fprintf(f_lambda,"----S--------S'--------S''--------t--------t_i--------t_i+1---|\n");
        fprintf(f_height,"----S--------S'--------S''--------t--------t_i--------t_i+1---|\n");
        fprintf(f_accel,"----a_x--------a_y--------a_z----|\n");
        fprintf(f_angul,"----w_x--------w_y--------w_z----|\n");
        for(i = 0; i < num - 1; i++)
        {
            do{
                A =(dataIn[i+1].t - t)*(dataIn[i+1].t - t)*(2*(t - dataIn[i].t) + (dataIn[i+1].t - dataIn[i].t));
                dir_A = 6*t*t + 2*(dataIn[i+1].t - dataIn[i].t - 2*dataIn[i].t - 4*dataIn[i+1].t)*t + 4*dataIn[i+1].t*dataIn[i].t - 2*dataIn[i+1].t*(dataIn[i+1].t - dataIn[i].t) + 2*dataIn[i+1].t*dataIn[i+1].t;
                dir2_A = 12*t + 2*((dataIn[i+1].t - dataIn[i].t) - 2*dataIn[i].t - 4*dataIn[i+1].t);
                B =(t - dataIn[i].t)*(t - dataIn[i].t)*(2*(dataIn[i+1].t - t) + (dataIn[i+1].t - dataIn[i].t));
                dir_B = -6*t*t + 2*(dataIn[i+1].t - dataIn[i].t + 2*dataIn[i+1].t + 4*dataIn[i].t)*t - (4*dataIn[i+1].t*dataIn[i].t + 2*dataIn[i].t*(dataIn[i+1].t - dataIn[i].t) + 2*dataIn[i].t*dataIn[i].t);
                dir2_B = -12*t + 2*((dataIn[i+1].t - dataIn[i].t) + 2*dataIn[i+1].t + 4*dataIn[i].t);
                phi_S = (dataIn[i].phi * A/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))) + (dataIn[i+1].phi * B/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)));
                phi_dirS = (dataIn[i].phi * dir_A/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))) + (dataIn[i+1].phi * dir_B/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)));
                phi_dir2S = (dataIn[i].phi * dir2_A/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))) + (dataIn[i+1].phi * dir2_B/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)));
                lambda_S = (dataIn[i].lambda * A/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))) + (dataIn[i+1].lambda * B/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)));
                lambda_dirS = (dataIn[i].lambda * dir_A/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))) + (dataIn[i+1].lambda * dir_B/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)));
                lambda_dir2S = (dataIn[i].lambda * dir2_A/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))) + (dataIn[i+1].lambda * dir2_B/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)));
                h_S = (dataIn[i].h * A/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))) + (dataIn[i+1].h * B/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)));
                h_dirS = (dataIn[i].h * dir_A/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))) + (dataIn[i+1].h * dir_B/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)));
                h_dir2S = (dataIn[i].h * dir2_A/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))) + (dataIn[i+1].h * dir2_B/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)));
                H = (dataIn[i].h*A)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))+(dataIn[i+1].h*B)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t));
                dir_H = (dataIn[i].h*dir_A)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))+(dataIn[i+1].h*dir_B)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t));
                dir2_H = (dataIn[i].h*dir2_A)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))+(dataIn[i+1].h*dir2_B)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t));
                fprintf(f_phi, "%lf %lf %lf %lf %lf %lf\n",phi_S, phi_dirS, phi_dir2S, t, dataIn[i].t, dataIn[i+1].t);
                fprintf(f_lambda, "%lf %lf %lf %lf %lf %lf\n",lambda_S, lambda_dirS, lambda_dir2S, t, dataIn[i].t, dataIn[i+1].t);
                fprintf(f_height, "%lf %lf %lf %lf %lf %lf\n",h_S, h_dirS, h_dir2S, t, dataIn[i].t, dataIn[i+1].t);
//                printf("%lf %lf %lf %lf %lf %lf\n",phi_S, phi_dirS, phi_dir2S, t, dataIn[i].t, dataIn[i+1].t);
//                printf("%lf %lf %lf %lf %lf %lf\n",lambda_S, lambda_dirS, lambda_dir2S, t, dataIn[i].t, dataIn[i+1].t);
//                printf("%lf %lf %lf %lf %lf %lf\n",h_S, h_dirS, h_dir2S, t, dataIn[i].t, dataIn[i+1].t);
                W = (dataIn[i].W*A)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))+(dataIn[i+1].W*B)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t));
                dir_W = (dataIn[i].W*dir_A)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))+(dataIn[i+1].W*dir_B)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t));
                dir2_W = (dataIn[i].W*dir2_A)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))+(dataIn[i+1].W*dir2_B)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t));
                Psi = (dataIn[i].psi*A)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))+(dataIn[i+1].psi*B)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t));
                dir_Psi = (dataIn[i].psi*dir_A)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))+(dataIn[i+1].psi*dir_B)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t));
                dir2_Psi = (dataIn[i].psi*dir2_A)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))+(dataIn[i+1].psi*dir2_B)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t));
                phi = (dataIn[i].phi*A)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t))+(dataIn[i+1].phi*B)/((dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t)*(dataIn[i+1].t - dataIn[i].t));
                teta = asin(dir_H/W);
                dir_teta = (1/(sqrt(1-((dir_H/W)*(dir_H/W)))))*(((dir2_H*W)-(dir_H*dir_W))/(W*W));
                double g = 9.780318*(1 + 0.0053024*sin(phi)*sin(phi) - 0.59*sin(2*phi)*(1e-5) - 0.3086*H*(1e-5)); //сила тяжести
                gamma = atan(-W*dir_Psi/g);
                dir_gamma = (1/(1+(W*Psi/g)*(W*Psi/g)))*(1/g)*(dir_W*dir_Psi + W*dir2_Psi);
                dir_We1 = -dir_W*sin(Psi)*cos(teta);
                dir_Wn1 = dir_W*cos(Psi)*cos(teta);
                dir_Wh1 = dir_W*sin(teta);
                We = -W*sin(Psi)*cos(teta);
                Wn = W*cos(Psi)*cos(teta);
                Wh = W*sin(teta);
                w_wave_E = dir_teta*cos(Psi) - dir_gamma*sin(Psi)*cos(teta);
                w_wave_N = dir_teta*sin(Psi) + dir_gamma*cos(Psi)*cos(teta);
                w_wave_H = dir_Psi + dir_gamma*sin(teta);
                dir_We2 = w_wave_N*Wh - w_wave_H*Wn;
                dir_Wn2 = w_wave_H*We - w_wave_E*Wh;
                dir_Wh2 = w_wave_E*Wn - w_wave_N*We;
                dir_Ve = dir_We1 + dir_We2;
                dir_Vn = dir_Wn1 + dir_Wn2;
                dir_Vh = dir_Wh1 + dir_Wh2;
                Ve = We + w_wave_E;//в методичке вообще не написано как искать эту величину
                Vn = Wn + w_wave_N;//в методичке вообще не написано как искать эту величину
                Vh = Wh + w_wave_H;//в методичке вообще не написано как искать эту величину
                double Azim = atan2(sin(dataIn[i+1].lambda - dataIn[i].lambda)*cos(dataIn[i+1].phi),cos(dataIn[i].phi)*sin(dataIn[i+1].phi) - sin(dataIn[i].phi)*cos(dataIn[i+1].phi)*cos(dataIn[i+1].lambda - dataIn[i].lambda));//азимутальный угол Bearing по этой ссылке(http://www.movable-type.co.uk/scripts/latlong.html)
                dir_v_ksi = dir_Vn*sin(Azim) + dir_Ve*cos(Azim);
                dir_v_eta = dir_Vn*cos(Azim) - dir_Ve*sin(Azim);
                dir_v_zeta = dir_Vh;
                v_ksi = Vn*sin(Azim) + Ve*cos(Azim);
                v_eta = Vn*cos(Azim) - Ve*sin(Azim);
                v_zeta = Vh;
                double b13 = cos(phi)*sin(dataIn[i].psi), b23 = cos(phi)*cos(dataIn[i].psi), b33 = sin(phi); //элементы МНК (вроде так ищутся, но я не совсем уверен)
                w_ksi = (-v_eta/a0)*(1 + e*e*(-0.5*b33*b33 + b23*b23) - H/a0) - v_ksi*e*e*pow(a0,-1)*b13*b23;
                w_eta = (v_ksi/a0)*(1 + e*e*(-0.5*b33*b33 + b23*b23) - H/a0) - v_eta*e*e*pow(a0,-1)*b13*b23;
                double G = a0*(1 + 0.5*sin(phi)*sin(phi) + H/a0);// радиус кривизны первого вертикала
                w_zeta = Ve*tan(Psi)/G;
                omega_ksi = omega*b13;
                omega_eta = omega*b23;
                omega_zeta = omega*b33;
                w_wave_ksi = w_wave_N*sin(Azim) + w_wave_E*cos(Azim);
                w_wave_eta = w_wave_N*cos(Azim) - w_wave_E*sin(Azim);
                w_wave_zeta = w_wave_H;
                w_abs_ksi = w_wave_ksi + w_ksi + omega_ksi;
                w_abs_eta = w_wave_eta + w_eta + omega_eta;
                w_abs_zeta = w_wave_zeta + w_zeta + omega_zeta;
                b_ksi = -v_zeta*(w_eta + 2*omega_eta) + 2*omega_zeta*v_eta;
                b_eta = -v_zeta*(w_ksi + 2*omega_ksi) + 2*omega_zeta*v_ksi;
                b_zeta = -v_eta*(w_ksi + 2*omega_ksi) +(w_eta + 2*omega_eta*v_ksi)*v_ksi - g;
                a_ksi = dir_v_ksi - b_ksi;
                a_eta = dir_v_eta - b_eta;
                a_zeta = dir_v_zeta - b_zeta;
                double psi_g = Psi - Azim; //гироскопический курс БЧЭ
                MNK[0][0] = cos(psi_g)*cos(gamma) - sin(psi_g)*cos(teta)*sin(gamma);
                MNK[1][0] = sin(psi_g)*cos(gamma) + cos(psi_g)*cos(teta)*sin(gamma);
                MNK[2][0] = sin(teta)*sin(gamma);
                MNK[0][1] = -cos(psi_g)*sin(gamma) - sin(psi_g)*cos(teta)*cos(gamma);
                MNK[1][1] = -sin(psi_g)*sin(gamma) - cos(psi_g)*cos(teta)*cos(gamma);
                MNK[2][1] = sin(teta)*cos(gamma);
                MNK[0][2] = sin(psi_g)*sin(teta);
                MNK[1][2] = -cos(psi_g)*sin(teta);
                MNK[2][2] = cos(teta);
                w_x = MNK[0][0]*w_abs_ksi + MNK[0][1]*w_abs_eta + MNK[0][2]*w_abs_zeta;
                w_y = MNK[1][0]*w_abs_ksi + MNK[1][1]*w_abs_eta + MNK[1][2]*w_abs_zeta;
                w_z = MNK[2][0]*w_abs_ksi + MNK[2][1]*w_abs_eta + MNK[2][2]*w_abs_zeta;
                a_x = MNK[0][0]*a_ksi + MNK[0][1]*a_eta + MNK[0][2]*a_zeta;
                a_y = MNK[1][0]*a_ksi + MNK[1][1]*a_eta + MNK[1][2]*a_zeta;
                a_z = MNK[2][0]*a_ksi + MNK[2][1]*a_eta + MNK[2][2]*a_zeta;
                fprintf(f_accel, "%lf %lf %lf\n", a_x, a_y, a_z);
                fprintf(f_angul, "%lf %lf %lf\n", w_x, w_y, w_z);
                printf("Acceleration X:%lf, Y:%lf, Z:%lf\n",a_x, a_y, a_z);
                printf("Angular velocity X:%lf, Y:%lf, Z:%lf\n", w_x, w_y, w_z);
                if(i == 0 && g == 0)
                {
                    Sleep(t);
                }
                if(i > 0 || g == 1){
                    t_new = t - t_old;
                    Sleep(t_new);
                }
                t_old = t;
                t += hag;
                g = 1;
            }while(t < dataIn[i+1].t);
        }
        fclose(f_phi);
        fclose(f_lambda);
        fclose(f_height);
    }
    else {
        printf("ERROR. Invalid number of points");
    }
}
