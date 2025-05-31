#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
//#include <matplot/matplot.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>

using namespace std;

const double PI = 3.141592653589793238462643;

vector<double> gen_sin(double A, double w, double krok, double x0, double x1, bool czy_pochodna, bool czy_wypisac) {
	vector<double> y;

	if (czy_pochodna == 0){ //zwykla funkcja
		for (double i = x0; i < x1; i += krok) {
		y.push_back(A * sin(w * i));
		}
	}
	else { //pochodna
		for (double i = x0; i < x1; i += krok) {
			y.push_back(A * cos(w * i));
		}
	}

	if (czy_wypisac == 1) {
		for (int i = 0; i < y.size(); i++) {
			cout << y[i] << ' ';
		}
	}

	/*vector<double> x;
	for (int i = 0; i < y.size(); i++)
		x.push_back(i);
	matplot::plot(x, y);
	matplot::show(); */
	return y;
}

vector<double> gen_cos(double A, double w, double krok, double x0, double x1, bool czy_pochodna, bool czy_wypisac) {
	vector<double> y;

	if (czy_pochodna == 0) { //zwykla funckja
		for (double i = x0; i < x1; i += krok) {
			y.push_back(A * cos(w * i));
		}
	}
	else { //pochodna
		for (double i = x0; i < x1; i += krok) {
			y.push_back((-1) * A * sin(w * i));
		}
	}
	/*vector<double> x;
	for (int i = 0; i < y.size(); i++)
		x.push_back(i);
	matplot::plot(x, y);
	matplot::show(); */
	if (czy_wypisac == 1) {
		for (int i = 0; i < y.size(); i++) {
			cout << y[i] << ' ';
		}
	}
	return y;
}

vector<double> gen_piloksztaltny(double T, double a, double czesc_okresu, double x0, double x1, bool czy_pochodna, bool czy_wypisac) {
	vector<double> y;
	double krok = T * czesc_okresu;
	if (czy_pochodna == 0) { //zwykly przebieg
 		for (double i = x0; i < x1; i += T){
			if (i > x1)
				break;
			for (double j = x0; j < T; j += krok) {
				y.push_back(a * j);
			}
		}
	}
	else { //pochodna
		for (double i = x0; i < x1; i += T) {
			if (i > x1)
				break;
			for (double j = x0; j < T; j += krok) {
				y.push_back(a);
			}
		}
	}

	if (czy_wypisac == 1) {
		for (int i = 0; i < y.size(); i++) {
			cout << y[i] << ' ';
		}
	}
	return y;
}

vector<double> gen_prostokatny(double T, double A, double czesc_okresu, double x0, double x1, bool czy_pochodna, bool czy_wypisac) {
	vector<double> y;
	double krok = T * czesc_okresu;
	
	if (czy_pochodna == 0) { //zwykly przebieg
		for (double i = x0; i < x1; i += T) {
			if (i > x1)
				break;
			for (double j = x0; j < T / 2; j += krok) {
				y.push_back(A);
			}
			for (double j = T / 2; j < T; j += krok) {
				y.push_back(-A);
			}
		}
	}
	else { //pochodna
		for (double i = x0; i < x1; i += T) {
			if (i > x1)
				break;
			for (double j = x0; j < T; j += krok) {
				y.push_back(0);
			}
		}
	}

	if (czy_wypisac == 1) {
		for (int i = 0; i < y.size(); i++) {
			cout << y[i] << ' ';
		}
	}
	return y;
}

vector<complex<double>> DTF(vector<double> funkcja) {
	vector<complex<double>> transformata;
	double e = 2.71828182845904;
	double PI = 3.14159265359;

	complex<double> temp;

	for (int k = 0; k < funkcja.size(); k++) {
		temp.imag(0);
		temp.real(0);
		for (int n = 0; n < funkcja.size(); n++) {
			temp.real(temp.real() + funkcja[k] * cos((-1) * 2 * PI * k * n / funkcja.size()));
			temp.imag(temp.imag() + funkcja[k] * sin((-1) * 2 * PI * k * n / funkcja.size()));
		}
		transformata.push_back(temp);
	}


	for (int i = 0; i < transformata.size(); i++) {
		cout << transformata[i].real() << transformata[i].imag() << "*i" << ' ';
	}
	return transformata;
}

vector<double> odw_DTF(vector<complex<double>> odwrotna) {
	vector<double> funkcja;
	double e = 2.71828182845904;

	double temp;

	for (int k = 0; k < odwrotna.size(); k++) {
		temp = 0;
		for (int n = 0; n < odwrotna.size(); n++) {
			temp += (odwrotna[k].real() / odwrotna.size()) * cos(2 * PI * k * n / odwrotna.size());
		}
		funkcja[k] = temp;
	}

	for (int i = 0; i < funkcja.size(); i++) {
		cout << funkcja[i] << ' ';
	}
	return funkcja;
}

vector<double> filtruj_srednia(vector<double> sygnal) {
	vector<double> wynik;
	int n = sygnal.size();

	if (n == 0)
		return wynik;

	// pierwszy punkt: srednia z pierwszego i drugiego
	wynik.push_back((sygnal[0] + sygnal[1]) / 2.0);

	// punkty srodkowe: pelna srednia z trzech
	for (int i = 1; i < n - 1; i++) {
		double srednia = (sygnal[i - 1] + sygnal[i] + sygnal[i + 1]) / 3.0;
		wynik.push_back(srednia);
	}

	// ostatni punkt: srednia z ostatniego i przedostatniego
	wynik.push_back((sygnal[n - 2] + sygnal[n - 1]) / 2.0);

	// wynik
	for (int i = 0; i < wynik.size(); i++) {
		cout << wynik[i] << ' ';
	}
	cout << endl;

	return wynik;
}


int main() {
	cout << "sin" << endl;
	gen_sin(1, 1, PI / 4, 0, 10 * PI, 0, 1);
	cout << endl << endl << endl;
	gen_sin(1, 1, PI / 4, 0, 10 * PI, 1, 1);
	cout << endl << endl << endl;

	cout << "cos" << endl;
	gen_cos(1, 1, PI / 4, 0, 10 * PI, 0, 1);
	cout << endl << endl << endl;
	gen_cos(1, 1, PI / 4, 0, 10 * PI, 1, 1);
	cout << endl << endl << endl;

	cout << "piloksztaltny" << endl;
	gen_piloksztaltny(1, 1, 0.1, 0, 10, 0, 1);
	cout << endl << endl << endl;
	gen_piloksztaltny(1, 1, 0.1, 0, 10, 1, 1);
	cout << endl << endl << endl;

	cout << "prostokatny" << endl;
	gen_prostokatny(1, 1, 0.1, 0, 10, 0, 1);
	cout << endl << endl << endl;
	gen_prostokatny(1, 1, 0.1, 0, 10, 1, 1);
	cout << endl << endl << endl;

	cout << "DTF" << endl;
	DTF(gen_sin(1, 1, PI / 4, 0, 10 * PI, 0, 0));
	cout << endl << endl << endl;
	DTF(gen_cos(1, 1, PI / 4, 0, 10 * PI, 0, 0));
	cout << endl << endl << endl;
	DTF(gen_piloksztaltny(1, 1, 0.1, 0, 10, 0, 0));
	cout << endl << endl << endl;
	DTF(gen_prostokatny(1, 1, 0.1, 0, 10, 0, 0));
	cout << endl << endl << endl;


	cout << "Filtrowanie sinusa" << endl;
	filtruj_srednia(gen_sin(1, 1, PI / 4, 0, 10 * PI, 0, 0));
	cout << endl;

	cout << "Filtrowanie cosinusa" << endl;
	filtruj_srednia(gen_cos(1, 1, PI / 4, 0, 10 * PI, 0, 0));
	cout << endl;

	cout << "Filtrowanie prostokatnego" << endl;
	filtruj_srednia(gen_prostokatny(1, 1, 0.1, 0, 10, 0, 0));
	cout << endl;

	cout << "Filtrowanie piloksztaltnego" << endl;
	filtruj_srednia(gen_piloksztaltny(1, 1, 0.1, 0, 10, 0, 0));
	cout << endl;
	return 0;
}

PYBIND11_MODULE(sygnaly, m) {
	m.def("gen_sin", &gen_sin, "generuje sinus (i jego pochodna)");
	m.def("gen_cos", &gen_cos, "generuje cosinus (i jego pochodna)");
	m.def("gen_piloksztaltny", &gen_piloksztaltny, "generuje sygnal piloksztaltny (i jego pochodna)");
	m.def("gen_prostokatny", &gen_prostokatny, "generuje sygnal prostokatny (i jego pochodna)");
	m.def("DTF", &DTF, "generuje Dyskretna Transformate Fouriera");
	m.def("odw_DTF", &odw_DTF, "generuje odwrotnosc Dyskretnej Transformaty Fouriera");
	m.def("filtruj_srednia", &filtruj_srednia, "generuje filtracje");
}