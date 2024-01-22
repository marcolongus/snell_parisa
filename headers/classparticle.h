#include "parameters.h"    //Módulo con parametros para la simulación.
using namespace std;

//Generador de números aleatorios en (0,1).
mt19937::result_type seed = time(0);
mt19937 gen(seed);                             //Standard mersenne_twister_engine seeded time(0)

uniform_real_distribution<KIND> dis(0., 1.); // dis(gen), número aleatorio real entre 0 y 1.
normal_distribution<KIND> normal_dist(0.0, 1.0);
/*Definimos la clase partículas y sus métodos */
/*Memoria que ocupa cada miembro de la clase: 4 doubbles (32 Bytes) + 1 int (4 Bytes) = 36 Bytes*/


/***************************************************************************************/
class particle {
	private:
		int state;
	public:
		KIND x,y;
		KIND velocity;
		KIND angle;

	particle();
	particle(KIND x1, KIND y1, KIND vel, KIND ang);

	bool is_healthy()    {return (state <  1);}
	bool is_infected()   {return (state == 1);}
	bool is_refractary() {return (state == 2);}

	int get_state()      {return state;}

	void set_healthy()    {state = 0;}
	void set_infected()   {state = 1;}
	void set_refractary() {state = 2;}
};


/* Constructor de partícula genérica*/
particle::particle() {
	velocity = 0;
	angle    = 0;
	x        = 0;
	y        = 0;
}

/* Constructor of a particle in a given phase-state (x,p) of the system */
particle::particle(KIND x1, KIND y1, KIND vel, KIND ang) {
	velocity = vel;
	angle    = ang;
	x        = x1;
	y        = y1;
}


particle create_particle(void) {
	KIND x, y, velocity, angle;
	x     = dis(gen) * L_x;
	y     = dis(gen) * L_y;
	angle = dis(gen) * dos_Pi;

	switch(velocity_distribution){
		case(0):
			velocity = -active_velocity * log(1. - dis(gen));
			break;
		case(1):
			velocity = pow(dis(gen) * (pow(v_max, 1 - k_powerl) - pow(v_min, 1 - k_powerl)) + pow(v_min, 1 - k_powerl), 1. / (1. - k_powerl));
			break;
		default:
			velocity = active_velocity;
			break;
	}
	
	// Creación de la partícula:
	particle A(x, y, velocity, angle);
	
	//Setting del estado interno de la partícula:
	// if   (dis(gen) < p_init){ A.set_infected();} // p_init de partículas infectadas.
	
	// INIT CONDITION	
	if  ( (0 <= A.x) and (A.x < x_wave)) A.set_infected(); // wave lineal.
	// if    ((pow(A.x - (L/2), 2) + pow(A.y - (L/2), 2)) <= 100)  {A.set_infected();} //wave radial
	// if  (((L/5 - 20) < A.x) and (A.x < L/5)) A.set_infected(); // wave lineal.
	// if (dis(gen) < p_rinit and !A.is_infected()) {A.set_refractary();}	
	else A.set_healthy();
	return A;
}


/* FUNCIONES AUXILIARES PARA LA CLASE */
/* Real boundary condition  and integer boundary condition functions */
KIND b_condition(KIND a, KIND L_mod) {
	return fmod((fmod(a, L_mod) + L_mod), L_mod);
}

int my_mod(int a, int b) {
	return ((a % b) + b) % b;
}


KIND distance(particle A, particle B) {
	KIND x1,x2,y1,y2,res;
	res = infinity;
	x2 = B.x; y2 = B.y;
	for(int i=-1; i<2; i++) for(int j=-1; j<2; j++) {
		x1 = A.x + i * L_x;
		y1 = A.y + j * L_y;
		res = min(res, (KIND)pow((x1 - x2), 2) + (KIND)pow((y1 - y2), 2));
	}
	return sqrt(res);
}

KIND distance_wall(particle A, particle B) {
	KIND x1, x2, y1, y2, res;
	x1 = A.x; y1 = A.y;
	x2 = B.x; y2 = B.y;
	res = (KIND)pow((x1 - x2), 2) + (KIND)pow((y1 - y2), 2);
	return sqrt(res);
}

KIND distance_x(particle A, particle B){
		KIND x1, x2, res;
		int j = 0;
		vector<KIND> dx;
		dx.resize(3,0);
		res = infinity;
		x2  = B.x;
		for(int i=-1; i<2; i++) {
			x1 = A.x + i * L_x;
			dx[i + 1] = x1 - x2;
			if (abs(dx[i + 1]) < res ) {
				res = abs(dx[i + 1]);
				j = i;
			} // if
		} // for
		return dx[j + 1];
}

KIND distance_y(particle A, particle B){
		KIND y1, y2, res;
		int j = 0;
		vector<KIND> dy;

		dy.resize(3,0);
		res = infinity;
		y2  = B.y;
		for(int i=-1; i<2; i++) {
			y1      = A.y + i * L_y;
			dy[i + 1] = y1 - y2;
			if (abs(dy[i + 1]) < res ) {
				res = abs(dy[i + 1]);
				j = i;
			} //if
		} // for
		return dy[j + 1];
}

KIND distance1(KIND dx, KIND dy) {
	return sqrt(pow(dx, 2) + pow(dy, 2));
}

bool interact(particle A, particle B) {
	return (distance(A, B) < diameter);
} 

bool interact_walls(particle A, particle B){
	return (distance_wall(A, B) < diameter);
} 



/* INTERACTION FUNCTIONS*/
/* Evolution time step function of the particle */
particle evolution(vector<particle> &system, vector<int> &index, bool inter){
	particle Agent = system[index[0]];
	
	bool right_regime = (Agent.x >= x_shift) and (Agent.y <=  L_y / delta * (Agent.x -  x_shift));
	
	/* DINÁMICA ESPACIAL DEL SISTEMA*/	
	if (right_regime) {
		Agent.angle += eta_right * normal_dist(gen) * sqrt_dt;  
	} else { 
		Agent.angle += eta_left * normal_dist(gen) * sqrt_dt;
	}

	/* PHYSICAL OVERLAP */
	inter = false;	
	if (inter) {
		vector<KIND> field, potencial;
		field.resize(2); potencial.resize(2,0); //inicia vector tamaño 2 en 0.

		for(size_t i=1; i < index.size(); i++) {
			KIND   dx_0i = distance_x(system[index[0]], system[index[i]]),
				   dy_0i = distance_y(system[index[0]], system[index[i]]),
				   d_0i  = distance1(dx_0i, dy_0i);
			potencial[0] += pow(d_0i,-3) * dx_0i;
			potencial[1] += pow(d_0i,-3) * dy_0i;

		} 
		
		for(size_t i=0; i<potencial.size(); i++) 
			potencial[i] = gamma_friction * potencial[i];

		field[0] = system[index[0]].velocity * cos(system[index[0]].angle) + potencial[0];
		field[1] = system[index[0]].velocity * sin(system[index[0]].angle) + potencial[1];

		Agent.x = b_condition(Agent.x + delta_time * field[0], L_x);
		Agent.y = b_condition(Agent.y + delta_time * field[1], L_y);
	} // if
	
	else {
		Agent.x = b_condition(Agent.x + Agent.velocity * cos(Agent.angle) * delta_time, L_x);
		Agent.y = b_condition(Agent.y + Agent.velocity * sin(Agent.angle) * delta_time, L_y);
	} // else


	// Reflective walls:
	if (Agent.x > L_x - 1 and cos(Agent.angle) > 0) Agent.angle = Pi - Agent.angle ;
	if (Agent.x < 1 and cos(Agent.angle) < 0) Agent.angle = Pi - Agent.angle;

	if (Agent.y > L_y - 1 and sin(Agent.angle) > 0) Agent.angle =  -Agent.angle;
	if (Agent.y < 1 and sin(Agent.angle) < 0) Agent.angle = - Agent.angle;


	/* DINÁMICA DE LA EPIDEMIA */
	bool flag = true; // Flag de infección.
	for (size_t i=1; i<index.size(); i++){
		if (Agent.is_healthy() && system[index[i]].is_infected()) {
			if (dis(gen) < p_transmision){
				Agent.set_infected();
				flag = false; // No puede volverse refractaria en esta instancia de evolución.
			}
		}
	} // for
	// if (Agent.is_refractary() && (dis(gen) < p_recfractary) ) Agent.set_healthy(); //SIRS
	// if (Agent.is_infected() && flag && (dis(gen) < p_infection) ) Agent.set_refractary();
	return Agent;
}
/***************************************************************************************/

