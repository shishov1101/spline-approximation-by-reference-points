#include <stdio.h>
#include "spline.h"

int main()
{
    int number;//количество эталонных точек
    int i;//счетчик цикла
    int f;//частота
    double t;//начальное время
    double deg_Phi;
    double deg_Lambda;
    double deg_Psi;
    printf("Enter number of points:");
    scanf("%d", &number); //Здесь выбирается количество эталонных точек
    printf("Enter the frequency of parameters output (The required frequency is more than 10 Hertz):");
    scanf("%d", &f); //Здесь выбирается частота выдаваемых параметров
    printf("Enter the starting time:");
    scanf("%lf", &t);
    data_input data_in[number];
    for(i = 0; i < number; i++) //инициализация эталонных точек
    {
        printf("Enter time 't' in ms value of number %d:", i);
        scanf("%lf", &(data_in[i].t));
        printf("Enter latitude in degrees'phi' value of number %d:", i);
        scanf("%lf", &deg_Phi);
        data_in[i].phi = deg_Phi*PI/180;
        printf("Enter longitude in degrees 'lambda' value of number %d:", i);
        scanf("%lf", &deg_Lambda);
        data_in[i].lambda = deg_Lambda*PI/180;
        printf("Enter height in metres'h' value of number %d:", i);
        scanf("%lf", &(data_in[i].h));
        printf("Enter in degrees'psi' value of number %d:", i);
        scanf("%lf", &deg_Psi);
        data_in[i].psi = deg_Psi*PI/180;
        printf("Enter 'W' value of number %d:", i);
        scanf("%lf", &(data_in[i].W));
    }
    spline(number, t, data_in);
    return  0;
}

