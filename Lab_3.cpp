#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <locale>
#include <sstream> 

using namespace std;

// Функция истинного решения

double true_sol(double x, double y) 
{
    return 1.0 - x * x - y * y;
}

// Функция правой части уравнения

double f(double x, double y) 
{
    return 4.0; 
}

// Граничные условия

double mu_left(double y) 
{
    return true_sol(-1.0, y);
}

double mu_right(double y) 
{
    return true_sol(1.0, y);
}

double mu_down(double x) 
{
    return true_sol(x, -1.0);
}

double mu_up(double x) 
{
    return true_sol(x, 1.0);
}

string format_num(double n)
{
    char buffer[20];
   
    snprintf(buffer, sizeof(buffer), "%*.*f", 10, 7, n);

    return string(buffer);
}

int main()
{

    // Параметры задачи
    setlocale(LC_ALL, "Russian");

    double a = -1.0, b = 1.0, c = -1.0, d = 1.0;
    int N = 4, M = 4;
    double h = (b - a) / N;
    double k = (d - c) / M;

    // Инициализация сетки (N+1 x M+1)
    vector<vector<double>> v(N + 1, vector<double>(M + 1, 0.0));

    // Задание граничных условий
    for (int j = 0; j <= M; j++)
    {
        double yj = c + j * k;
        v[0][j] = mu_left(yj);
        v[N][j] = mu_right(yj);
    }

    for (int i = 0; i <= N; i++)
    {
        double xi = a + i * h;
        v[i][0] = mu_down(xi);
        v[i][M] = mu_up(xi);
    }

    // Коэффициенты для разностной схемы
    double coeff_h = 1.0 / (h * h);
    double coeff_k = 1.0 / (k * k);
    double A = -2.0 * (coeff_h + coeff_k);

    // Параметры итерационного процесса
    double eps = 1e-6;
    int Nmax = 1000;
    int S = 0;
    bool end = false;

    // Для хранения значений внутренних узлов на каждой итерации
    vector<vector<vector<double>>> v_history;
    v_history.push_back(v);

    bool endByEps = false;
    bool endByIter = false;
    double eps_max = 0.0;

    // Итерационный процесс (метод Зейделя)
    while (!end)
    {
        S++;
        eps_max = 0.0;

        // Проход по внутренним узлам
        for (int j = 1; j < M; j++)
        {
            for (int i = 1; i < N; i++)
            {
                double v_old = v[i][j];

                double rhs = -f(a + i * h, c + j * k);

                rhs -= coeff_h * (v[i - 1][j] + v[i + 1][j]);
                rhs -= coeff_k * (v[i][j - 1] + v[i][j + 1]);

                v[i][j] = rhs / A;

                // Контроль изменения 
                double curr_diff = fabs(v[i][j] - v_old);

                if (curr_diff > eps_max)
                {
                    eps_max = curr_diff;
                }
            }
        }
        // Запись текущего состояния
        v_history.push_back(v);

        // Критерии остановки
        if (eps_max < eps)
        {
            endByEps = true;
            end = true;
        }
        if (S >= Nmax)
        {
            endByIter = true;
            end = true;
        }
    }

    // ______ Вывод результатов ______
    if (endByEps)
    {
        cout << "\nЗакончили по критерию точности" << endl;
    }
    if (endByIter)
    {
        cout << "\nЗакончили по критерию количества итераций\n" << endl;
    }

    cout << "При решении разностной схемы методом Зейделя с критериями остановки Nmax = " << Nmax
        << " и eps = " << eps << " за N = " << S << " итераций\nдостигнута точность eps_max = " << eps_max
        << " и получено численное решение:\n";

    // Выводим численное решение
    int n = (M - 1) * (N - 1);

    vector<double> x(n);

    for (int j = 0; j < M - 1; j++)
    {
        for (int i = 0; i < N - 1; i++)
        {
            x[j * (N - 1) + i] = v[i + 1][j + 1];
        }
    }

    cout << "Численное решение\n\n";

    cout << "-------------------------------------------------------------------------------" << endl;

    for (int row = 0; row <= M; row++) 
    {
        int i = M - row;  

        cout << " " << format_num(c + i * k) << " |";

    
        for (int j = 0; j <= N; j++) 
        {
            cout << " " << format_num(v[i][j]) << " |";
        }
        cout << endl;

        cout << "-------------------------------------------------------------------------------" << endl;
    }

    string header = "   y \\ x    |";

    for (int i = 0; i <= N; i++)
    {

        header += " " + format_num(a + i * h) + " |";
    }
    cout << header << endl;

    cout << "-------------------------------------------------------------------------------" << endl;


    // Считаем невязку 
    double max_r = 0.0; // Норма Чебышева ||R||беск = max|R_ij|

    vector<double> residuals;

    for (int i = 1; i < N; i++)
    {
        for (int j = 1; j < M; j++) 
        {
            double u = (v[i - 1][j] - 2 * v[i][j] + v[i + 1][j]) / (h * h)
                + (v[i][j - 1] - 2 * v[i][j] + v[i][j + 1]) / (k * k);

            double r_ij = u + f(a + i * h, c + j * k);

            // Норма Чебышева

            if (fabs(r_ij) > max_r) 
            {
                max_r = fabs(r_ij);
            }
            residuals.push_back(r_ij);
        }
    }

    cout << "Норма невязки (Чебышева): ||R||беск = " << max_r << endl;
 
    //double res_norm = 0.0;

    //for (double r : residuals) 

    //{
    //    res_norm += r * r;
    //}

    //res_norm = sqrt(res_norm);

    //cout << "\nНорма невязки: " << res_norm << endl;
    //cout << "\nМаксимум невязки: " << max_r << endl;

    // Норма общей погрешности  
    vector<double> trueSol;  // Объявление вектора trueSol

    // Заполняем trueSol значениями точного решения для внутренних узлов
    for (int j = 1; j < M; j++) 
    {
        for (int i = 1; i < N; i++) 
        {
            trueSol.push_back(true_sol(a + i * h, c + j * k)); 
        }
    }
    // Разность между численнымм и точным решением в каждом внутреннем узле сетки
    vector<double> resTrue(x.size());

    for (size_t i = 0; i < x.size(); i++)
    {
         resTrue[i] = x[i] - trueSol[i];
    }

    // Вывод точного решения
    cout << "\nТочное решение" << endl;  
    cout << endl;

    cout << "-----------------------------------------------------------------------------\n";

    for (int j = M; j >= 0; j--)
    {
   
        string row = format_num(c + j * k) + " | ";

        for (int i = 0; i <= N; i++)
        {
            row += format_num(true_sol(a + i * h, c + j * k)) + " | ";
        }
        cout << row << endl;
    }
    cout << "-----------------------------------------------------------------------------\n";

    header = "   y \\ x    |";

    for (int i = 0; i <= N; i++)
    {
        header += format_num(a + i * h) + " | ";
    }
    cout << header << endl;

    cout << string(header.length(), '-') << endl;


    // Вычисляем погрешность
    double op = 0;

    for (int i = 0; i < N + 1; i++)
    {
        for (int j = 0; j < M + 1; j++) 
        {
            op = max(op, fabs(v[i][j] - true_sol(a + i * h, c + j * k)));
        }
    }
    cout << "Норма погрешности (Чебышева): ||Z||беск = " << op << endl;

    // Численное решение после первой итерации
    cout << "\nЧисленное решение после первой итерации" << endl;
    cout << endl;
    cout << "-----------------------------------------------------------------------------" << endl;

    auto v_iter1 = v_history[1];

    for (int row = 0; row <= M; row++)
    {
        int j = M - row;  

        string row_str = format_num(c + j * k) + " | ";

        for (int i = 0; i <= N; i++)
        {
            row_str += format_num(v_iter1[i][j]) + " | ";
        }
        cout << row_str << endl;
    }
    cout << "-----------------------------------------------------------------------------" << endl;

    header = "   y \\ x   |";

    for (int i = 0; i <= N; i++)
    {
        header += " " + format_num(a + i * h) + " |";
    }
    cout << header << endl;

    cout << string(header.length(), '-') << endl;


    // Численное решение после второй итерации  
    cout << "\nЧисленное решение после второй итерации" << endl;
    cout << endl;
    cout << "-----------------------------------------------------------------------------" << endl;

    auto v_iter2 = v_history[2];

    for (int row = 0; row <= M; row++)
    {
        int j = M - row; 

        string row_str = format_num(c + j * k) + " | ";

        for (int i = 0; i <= N; i++)
        {
            row_str += format_num(v_iter2[i][j]) + " | ";
        }
        cout << row_str << endl;
    }
    cout << "-----------------------------------------------------------------------------" << endl;

    header = "   y \\ x   |";

    for (int i = 0; i <= N; i++)

    {
        header += " " + format_num(a + i * h) + " |";
    }

    cout << header << endl;

    cout << string(header.length(), '-') << endl;

}

