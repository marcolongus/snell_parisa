//System: Not interacting Persistent Random Walkers. 

#include <bits/stdc++.h>
#include "./headers/agentes.h"

#define CHECK(val_1, val_2) if (val_1 % val_2 == 0)
#define TIME(val_1, val_2) static_cast<float>(val_1) * val_2

using namespace std;

int main(void) {
	
	/* DEFINICIÓN DE ARCHIVOS DE SALIDA DEL PROGRAMA */
	ofstream final_state("data/evolution.txt", ios_base::app);
	ofstream epidemic("data/epidemia.txt", ios_base::app);
	ofstream metrica("data/metrica.txt", ios_base::app);
	ofstream simulation_data("data/simulation_data.txt", ios_base::app);
	//ofstream anim("data/animacion.txt", ios_base::app);
	print_simulation_parameters(simulation_data);

	/* METRICA Y PERFORMANCE */
	float  updates = 0;  // Contador de updates.
	size_t start_s = clock();

	/* SIMULACION */
	gen.seed(seed);
	for (size_t n_simulaciones = 0; n_simulaciones < 20; n_simulaciones++) {
		print_header(n_simulaciones);
		
		string anim_file = "data/animacion_" + to_string(n_simulaciones) + ".txt";
		ofstream anim(anim_file, ios_base::app);
		
		/* DECLARACIÓN DE VARIABLES */
		vector<particle> system,
						 system_new;
		vector<bool>     inter;        	// Flag de interacción.
		vector<size_t>   state_vector;  // (S, I, R) vector.
		
		inter.resize(N, false);
		state_vector.resize(spin, 0);

		/* INICIALIZAMOS GRILLA */
		vector<vector<set<size_t>>> grid;
		size_t num_grid_x = floor(L_x);
		size_t num_grid_y = floor(L_y);

		grid.resize(num_grid_x);
		for (size_t i = 0; i < grid.size(); i++) grid[i].resize(num_grid_y);

		/* CONDICIÓN INICIAL */
		init_system(system, state_vector, grid);
		system_new.resize(system.size());
		print_state(state_vector);
		
		/* EVOLUCIÓN DEL SISTEMA */
		int time_step = 0;
		KIND i_max = 0, t_max = 0;
		while (state_vector[1] < N) {
			CHECK(time_step, (int)2e04) print_epidemic_tofile(epidemic, state_vector, time_step);
			CHECK(time_step, (int)2e04) printf("Time: %0.f\n", TIME(time_step, delta_time));
			time_step++;
			update_system(system, system_new, state_vector, grid, inter, time_step, anim);
			if (i_max < state_vector[1]) {
				i_max = state_vector[1];
				t_max = TIME(time_step, delta_time);
			}
		}  // WHILE
		updates += static_cast<float>(time_step);

		/* ESCRITURA DE RESULTADOS */
		print_finalstate_tofile(final_state, state_vector, i_max, t_max, time_step);
		print_result_header();
		print_state(state_vector);
		cout << endl;
		anim.close();
	}  // FOR SIMUL

	/* IMPRESION EN PANTALLA */
	int   stop_s       = clock();
	float cps          = static_cast<float>(CLOCKS_PER_SEC);
	float time_elapsed = static_cast<float>(clock() - start_s) / cps;
	float metric       = updates*static_cast<float>(N) / time_elapsed;

	cout << "Time [seg]   : " << time_elapsed  		<< endl;
	cout << "Time [min]   : " << time_elapsed / 60  << endl;
	cout << "Metrica [pps]: " << metric       		<< endl;
	
	metrica  << metric << endl;

	/* CERRRAMOS ARCHIVOS */
	final_state.close();
	epidemic.close();
	metrica.close();
	simulation_data.close();
	// anim.close();

	return 0;
}

/* Autor: Benjamín R. Marcolongo. FAMAF-UNC.
/* Email: benjaminmarcolongo@gmail.com / benjaminmarcolongo@unc.edu.ar 
 *-----------------------------------------------------------------------------------------------------------
 * Programa para evolucionar un sistema de N agentes en el tiempo.
 *-----------------------------------------------------------------------------------------------------------
 * Agentes:
 *		i.   Las velocidades de los N agentes pueden ser uniformes o tomadas aleatoriamente.
 *		      1. Las distribuciones por defecto son de tipo exponencial o de ley de potencias.
 *		ii.  Pueden interactuar a través de un potencial u otra forma (una red neuronal i.e).
 *			  1. El potencial por defecto es de esferas blandas.
 *		ii.  Están confinados a un toroide topológico-
 *			  1. El tamaño característico del toro se denota con L (lado del cuadrado).
 *		iii. Poseen un estado interno, caracterizado por un número entero.
 *			  1. Este estado puede o no interactaur con la dinámica espacial.
 *		iv.  El estado interno puede evolucionar en el tiempo.
 *			  1. Esta dinámica está regulada, genéricamente, con una distribución de poisson.
 *-----------------------------------------------------------------------------------------------------------
 * Red compleja: (falta programar, esta en versiones anteriores)
 *			i.  El programa trackea todo el sistema de interacciones y guarda la red compleja resultante.
 *			ii. La red compleja puede ser la asociada a la propagación del estado interno o la de contactos
 *-----------------------------------------------------------------------------------------------------------
 * Estuctura de datos para optimizar la búsqueda de interacciones entre agentes:
 *	   i.   Utiliza un red-and-black tree (BST) implementado en c++ como set.
 			Esta estructura es equivalente a un linked list.
 *	   ii.  Cada agente está indexado por un int.
 *	   iii. Se construye una grilla con cuadrículas de tamaño 1x1 y cada a una se le asigna un set.
 *	   iv.  Cada set contiene los agentes que están en cada cuadrícula proyectando su posción sobre la grilla.
 *----------------------------------------------------------------------------------------------------------- 
*/
