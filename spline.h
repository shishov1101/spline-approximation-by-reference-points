#define PI 3.14159265358979
//структура, состоящая из входных данных
typedef struct{
    double t;//время
    double phi;//широта
    double lambda;//долгота
    double h;//высота
    double psi;//истинный курс
    double W;//траекторная скорость
}data_input;

void* spline(int num, double t, data_input *dataIn);
