#include "classparticle.h"
#define forn(i,a,b) for(int i=a; i<b; i++)

using namespace std;

void
init_system( vector<particle>            &system       ,  
			 vector<size_t>              &state_vector , 
			 vector<vector<set<size_t>>> &grid            )
{
	for(int p = 0; p < N; p++) {
		particle Agent;
		bool accepted = false;
		while(!accepted){
			accepted = true;
			Agent = create_particle();
			int i_index = floor(Agent.x),
				j_index = floor(Agent.y);
			forn(l, -2, 3) {
				forn(m, -2, 3) {
					int i = b_condition(i_index + l, L_x),
						j = b_condition(j_index + m, L_y);
					
					if (!grid[i][j].empty()) {
						for (auto element: grid[i][j]) {
							if (interact(Agent, system[element])) accepted = false;
						} 
					}
				}
			}
			if (accepted) grid[i_index][j_index].insert(p); 
		} // while
		system.push_back(Agent);
		state_vector[Agent.get_state()]++;
	} // For N	
}


void 
update_system(  vector<particle>            &system      , 
				vector<particle>            &system_new  ,
				vector<size_t>              &state_vector, 
				vector<vector<set<size_t>>> &grid        ,
				vector<bool>                &inter       ,
				int                         &time_step   ,
				ofstream                    &anim          )
{
	size_t healthy=0, infected=0, refract=0;
	state_vector = {0UL,0UL,0UL};
	#pragma omp parallel for schedule(guided) reduction(+:healthy, infected, refract) num_threads(6)
	for (size_t p=0; p<N; p++) {
			vector<int> index;
			index.push_back(p);
			inter[p] = false;
			/* chequeamos interacciones */
			forn(l, -2, 3) forn(m, -2, 3) {
				size_t i_index = b_condition(floor(system[p].x) + l, L_x),
					   j_index = b_condition(floor(system[p].y) + m, L_y);
				if(!grid[i_index][j_index].empty()) {
					for(auto element: grid[i_index][j_index]) {
						if (element !=p && interact_walls(system[p], system[element])) {
							inter[p] = true;
							index.push_back(element);
						}
					}//for
				}//if not empty
			} //for m, l
			/*fin de chequeo de interacciones*/
			system_new[p] = evolution(system, index, inter[p]);
			switch(system_new[p].get_state()){
				case(0):
					healthy++;
					break;
				case(1):
					infected++;
					break;
				default:
				refract++;
					break;
			}
		}//for p

		state_vector = {healthy, infected, refract};
		
		//Animacion (records prev. system state):
		if (animation and (time_step - 1) % anim_step == 0) {
			for(size_t p=0; p < N; p++) {
				anim << system[p].x                   << " ";
				anim << system[p].y                   << " ";
				anim << (time_step - 1) * delta_time  << " ";
				anim << system[p].get_state() << endl;
			}
		} // If animacion
		
		/* Update set */
		for(size_t p=0; p<N; p++) {
			int i_new = floor(system_new[p].x),
				j_new = floor(system_new[p].y);
			
			int i_old = floor(system[p].x),
				j_old = floor(system[p].y);

			if (grid[i_new][j_new].find(p) == grid[i_new][j_new].end()){
				grid[i_old][j_old].erase(p);
				grid[i_new][j_new].insert(p);
			} 
		} // For p set.
		system = system_new;
}


//Print functiones
void print_header(int n_simulaciones)
{
	cout << "=========================================================="   << endl;
	cout << "[+] Simulacion: " << n_simulaciones << endl;
	cout << "=========================================================="   << endl;
	cout << endl;
}

void print_result_header(void)
{
	cout << endl;
	cout << "Experimento data:"    << endl;
}


void print_mem_info(void)
{
	int   system_memory  = 2 * ( 5 * sizeof(KIND)*N);
	float space_memory   = 48.f * L * L;
	cout << "Memoria del sistema [Bytes]     : " << system_memory                 << " Bytes" << endl;
	cout << "Memoria del sistema [Megasbytes]: " << (float)system_memory / 1000000. << " Mb"    << endl;  
	cout << "Memoria de una particula 		 : " << (float)system_memory / (float)N << " Bytes" << endl;
	cout << "Memoria del espacio      		 : " <<  space_memory / 1e06f     	  << " Mb"    << endl;
	cout << "Memoria de un grid        		 : " <<  48                           << " Bytes" << endl;  
}


void print_epidemic_tofile(ofstream &file, vector<size_t> &state_vector, int &time_step)
{
	file << state_vector[0] << " ";
	file << state_vector[1] << " ";
	file << state_vector[2] << " ";
	file << delta_time * (double)time_step << endl;
}


void print_state(vector<size_t> state_vector){
		cout << "Healthy   : " << state_vector[0] << endl;
		cout << "Infected  : " << state_vector[1] << endl;
		cout << "Refractary: " << state_vector[2] << endl; 
		cout << endl;
}


void print_finalstate_tofile(ofstream       &file, 
							 vector<size_t> &state_vector, 
							 KIND 			&i_max, 
							 KIND 			&t_max, 
							 int 			&time_step)
{
	file << state_vector[2] << " ";
	file << delta_time*(KIND)time_step << " ";
	file << i_max << " ";
	file << t_max << endl;
}

void print_simulation_parameters(ofstream &file)
{
	file << "L               = "   << L 		 			 << endl;
	file << "L_x             = "   << L_x 		 			 << endl;
	file << "L_y             = "   << L_y 		 			 << endl;
	file << "N               = "   << N 		 			 << endl;
	file << "dt              = "   << delta_time 			 << endl;
	file << "active_vel      = "   << active_velocity 		 << endl;
	file << "Vel Distribut   = "   << velocity_distribution  << endl;
	file << "tau_left rot    = "   << alpha_left 			 << endl;
	file << "tau_right rot   = "   << alpha_right 			 << endl;
	file << "tau_t           = "   << tau_t 				 << endl;
	file << "tau_i           = "   << tau_i 				 << endl;
	file << "tau_r           = "   << tau_r 				 << endl;
	file << "delta           = "   << delta 				 << endl;
	file << "anim step time  = "   << anim_step * delta_time << endl;
	file << "anim_step       = "   << anim_step 			 << endl;
	file << "density         = "   << N / pow(L, 2)			 << endl;
	file << "seed            = "   << seed 					 << endl;
}